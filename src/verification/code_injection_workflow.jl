"""
High-level MMS verification workflow with code injection

This module provides the main user-facing workflow for MMS verification
with automatic code generation and injection.
"""

using Printf

"""
    MMSProblem

Complete MMS problem specification for code injection workflow.
"""
struct MMSProblem{P,M,V}
    pde::P                              # Symbolic PDE
    manufactured_solution::M             # Symbolic manufactured solution  
    variables::V                         # Symbolic variables (x, t, u, etc.)
    domain::Tuple                        # Spatial domain
    time_span::Union{Tuple,Nothing}      # Time span (if time-dependent)
    
    function MMSProblem(pde, mms, variables, domain; time_span=nothing)
        new{typeof(pde),typeof(mms),typeof(variables)}(
            pde, mms, variables, domain, time_span
        )
    end
end

"""
    SolverSetup

Configuration for user's solver that will be modified.
"""
struct SolverSetup
    source_file::String                  # Path to original solver code
    language::Symbol                     # :julia, :c, :python, :fortran
    output_file::String                  # Where to write modified code
    result_file::Union{String,Nothing}   # Where solver writes results
    run_command::Union{String,Nothing}   # How to run the solver
    
    function SolverSetup(source_file::String, language::Symbol;
                        output_file::String = source_file * ".mms",
                        result_file::Union{String,Nothing} = nothing,
                        run_command::Union{String,Nothing} = nothing)
        new(source_file, language, output_file, result_file, run_command)
    end
end

"""
    verify_with_code_injection(
        problem::MMSProblem,
        solver_setup::SolverSetup;
        resolutions::Vector = [0.1, 0.05, 0.025],
        result_reader::Union{Function,Nothing} = nothing,
        verbose::Bool = true
    )

Complete MMS verification workflow with automatic code generation and injection.

# Workflow
1. Compute forcing term symbolically
2. Generate forcing function code in target language
3. Inject code into solver file
4. Run solver at multiple resolutions
5. Analyze convergence

# Arguments
- `problem::MMSProblem`: Complete problem specification
- `solver_setup::SolverSetup`: Solver configuration
- `resolutions::Vector`: Grid resolutions to test
- `result_reader::Function`: Function to read results (path -> solution_array)
- `verbose::Bool`: Print progress

# Example
```julia
# Define problem
@variables x u(x)
Dx = Differential(x)
problem = MMSProblem(
    -Dx(Dx(u)),                    # PDE
    sin(π*x),                      # Manufactured solution
    (x, u),                        # Variables
    (0.0, 1.0)                     # Domain
)

# Setup solver
solver = SolverSetup(
    "my_poisson_solver.c",         # Original C code
    :c,                            # Language
    output_file = "my_poisson_solver_mms.c",
    run_command = "gcc my_poisson_solver_mms.c -o solver && ./solver"
)

# Run verification
results = verify_with_code_injection(problem, solver, resolutions=[0.1, 0.05])
```
"""
function verify_with_code_injection(
    problem::MMSProblem,
    solver_setup::SolverSetup;
    resolutions::Vector = [0.1, 0.05, 0.025],
    result_reader::Union{Function,Nothing} = nothing,
    verbose::Bool = true,
    additional_params::Dict = Dict()
)
    verbose && @info "MMS Verification with Code Injection" language=solver_setup.language
    
    # Step 1: Compute forcing term symbolically
    verbose && @info "Step 1: Computing forcing term symbolically"
    
    forcing_term = if problem.time_span === nothing
        # Spatial only (elliptic PDE)
        compute_forcing_term_1d(
            problem.pde, 
            problem.manufactured_solution,
            problem.variables[1],  # x
            problem.variables[2]   # u
        )
    else
        # Time-dependent
        compute_forcing_term_time_dependent(
            problem.pde,
            problem.manufactured_solution,
            problem.variables[1],  # x
            problem.variables[2],  # t
            problem.variables[3]   # u
        )
    end
    
    if verbose
        println("  Forcing term: ", forcing_term.symbolic_form)
    end
    
    # Step 2: Generate code
    verbose && @info "Step 2: Generating $(solver_setup.language) code"
    
    arg_names = if problem.time_span === nothing
        [:x]
    else
        [:x, :t]
    end
    
    forcing_code = generate_forcing_code(forcing_term, solver_setup.language, arg_names)
    
    if verbose
        println("Generated code:")
        println(forcing_code)
    end
    
    # Step 3: Inject into solver
    verbose && @info "Step 3: Injecting code into solver"
    
    marker_found = inject_forcing_code(
        solver_setup.source_file,
        forcing_code,
        solver_setup.output_file
    )
    
    if marker_found
        verbose && @info "  Code injected successfully at marker location"
    else
        @warn "  Marker not found - code appended to end of file"
    end
    
    verbose && @info "  Modified solver written to: $(solver_setup.output_file)"
    
    # Step 4: Run convergence study
    verbose && @info "Step 4: Running convergence study"
    
    convergence_data = ConvergenceData{Float64}[]
    
    for (idx, h) in enumerate(resolutions)
        verbose && @info "  Resolution $idx/$(length(resolutions)): h = $h"
        
        # Create grid
        grid = Grid1D(problem.domain..., h)
        verbose && @info "    Grid: $(grid.nx) points"
        
        # Modify solver parameters if needed
        modified_file = replace_parameters(
            solver_setup.output_file,
            solver_setup.output_file * ".run",
            h,
            additional_params
        )
        
        # Run solver
        if solver_setup.run_command !== nothing
            verbose && @info "    Running solver..."
            t_start = time()
            
            run_solver_command(solver_setup.run_command, modified_file)
            
            t_elapsed = time() - t_start
            verbose && @info "    Completed in $(round(t_elapsed, digits=3))s"
        else
            @warn "No run_command specified - skipping solver execution"
            t_elapsed = 0.0
        end
        
        # Read results
        if result_reader !== nothing && solver_setup.result_file !== nothing
            numerical_solution = result_reader(solver_setup.result_file)
            
            # Evaluate manufactured solution
            analytical_solution = if problem.time_span === nothing
                [evaluate_manufactured(problem.manufactured_solution, problem.variables[1] => x_val) 
                 for x_val in grid.x]
            else
                # For time-dependent, evaluate at final time
                t_final = problem.time_span[2]
                [evaluate_manufactured(problem.manufactured_solution, 
                                      problem.variables[1] => x_val,
                                      problem.variables[2] => t_final) 
                 for x_val in grid.x]
            end
            
            # Compute errors
            errors = compute_all_errors(numerical_solution, analytical_solution, grid)
            
            if verbose
                for (norm_name, error_value) in errors
                    @printf "    %s error: %.6e\n" norm_name error_value
                end
            end
            
            push!(convergence_data, ConvergenceData(h, grid.nx, errors, t_elapsed))
        else
            @warn "No result_reader specified - skipping error analysis for this resolution"
        end
    end
    
    # Step 5: Analyze convergence
    if !isempty(convergence_data)
        verbose && @info "Step 5: Analyzing convergence"
        
        convergence_rates = analyze_convergence(convergence_data)
        
        if verbose
            @info "Convergence Analysis:"
            for (key, value) in convergence_rates
                if !endswith(string(key), "_R2")
                    R2_key = Symbol(string(key) * "_R2")
                    R2 = haskey(convergence_rates, R2_key) ? convergence_rates[R2_key] : NaN
                    @printf "  %s: rate = %.3f (R² = %.4f)\n" key value R2
                end
            end
        end
        
        # Create result
        problem_info = Dict{Symbol,Any}(
            :pde => string(problem.pde),
            :manufactured_solution => string(problem.manufactured_solution),
            :domain => problem.domain
        )
        
        solver_info = Dict{Symbol,Any}(
            :source_file => solver_setup.source_file,
            :language => solver_setup.language,
            :output_file => solver_setup.output_file
        )
        
        return VerificationResult(
            problem_info,
            solver_info,
            convergence_data,
            convergence_rates
        )
    else
        @warn "No convergence data collected - cannot perform analysis"
        return nothing
    end
end

"""
Helper function to evaluate manufactured solution at a point
"""
function evaluate_manufactured(expr, pairs...)
    # This needs Symbolics substitution
    # For now, simplified version
    error("evaluate_manufactured needs implementation with Symbolics.substitute")
end

"""
Helper function to replace parameters in solver file
"""
function replace_parameters(input_file::String, output_file::String, h::Float64, params::Dict)
    content = read(input_file, String)
    
    # Common parameter replacements
    content = replace(content, r"MMS_H_VALUE\s*=\s*[\d.e+-]+" => "MMS_H_VALUE = $h")
    
    # Additional custom parameters
    for (key, value) in params
        pattern = Regex("$(key)\\s*=\\s*[\\d.e+-]+")
        content = replace(content, pattern => "$(key) = $value")
    end
    
    write(output_file, content)
    return output_file
end

"""
Helper function to run solver command
"""
function run_solver_command(command::String, modified_file::String)
    # Replace placeholder if present
    actual_command = replace(command, "SOLVER_FILE" => modified_file)
    run(`sh -c $actual_command`)
end

export MMSProblem, SolverSetup, verify_with_code_injection
