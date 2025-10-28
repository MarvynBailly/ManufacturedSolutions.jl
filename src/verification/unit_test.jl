"""
Unit testing interface for MMS verification

This module provides a simple unit test workflow:
1. Inject forcing term into user's code
2. Run the modified code
3. Verify solution against manufactured solution
"""

using DelimitedFiles
using Printf

"""
    MMSUnitTest

Configuration for an MMS unit test.
"""
struct MMSUnitTest
    test_name::String
    source_file::String              # User's Julia code file to test
    pde::Any                         # Symbolic PDE
    manufactured_solution::Any       # Symbolic manufactured solution
    variables::Tuple                 # (x, u) or (x, t, u)
    domain::Tuple{Float64,Float64}   # Spatial domain
    resolution::Float64              # Grid spacing for test
    output_file::String              # Where solver writes results
    marker::String                   # Injection marker in source
    tolerance::Float64               # Error tolerance for pass/fail
    
    function MMSUnitTest(
        test_name::String,
        source_file::String,
        pde,
        manufactured_solution,
        variables::Tuple,
        domain::Tuple;
        resolution::Float64 = 0.01,
        output_file::String = "mms_test_output.txt",
        marker::String = "MMS_FORCING",
        tolerance::Float64 = 1e-3
    )
        new(test_name, source_file, pde, manufactured_solution, variables,
            domain, resolution, output_file, marker, tolerance)
    end
end

"""
    run_mms_unit_test(test::MMSUnitTest; verbose::Bool=true)

Run a single MMS unit test.

# Workflow
1. Compute forcing term symbolically
2. Generate Julia code for forcing function
3. Inject into source file (creates modified version)
4. Run modified code
5. Read output
6. Compare against manufactured solution
7. Return pass/fail with error metrics

# Returns
- `Dict` with keys: `:passed`, `:L2_error`, `:Linf_error`, `:message`

# Example
```julia
@variables x u(x)
Dx = Differential(x)

test = MMSUnitTest(
    "Poisson Solver Test",
    "my_solver.jl",
    -Dx(Dx(u)),
    sin(π*x),
    (x, u),
    (0.0, 1.0),
    resolution = 0.01,
    tolerance = 1e-4
)

result = run_mms_unit_test(test)
if result[:passed]
    println("✓ Test PASSED")
else
    println("✗ Test FAILED: ", result[:message])
end
```
"""
function run_mms_unit_test(test::MMSUnitTest; verbose::Bool=true)
    verbose && println("="^70)
    verbose && println("MMS Unit Test: $(test.test_name)")
    verbose && println("="^70)
    
    # Step 1: Compute forcing term
    verbose && println("\n[1/6] Computing forcing term symbolically...")
    
    forcing = if length(test.variables) == 2
        # Spatial only
        x_var, u_var = test.variables
        compute_forcing_term_1d(test.pde, test.manufactured_solution, x_var, u_var)
    else
        # Time-dependent
        x_var, t_var, u_var = test.variables
        compute_forcing_term_time_dependent(test.pde, test.manufactured_solution, 
                                           x_var, t_var, u_var)
    end
    
    if verbose
        println("  Forcing: ", forcing.symbolic_form)
    end
    
    # Step 2: Generate Julia code
    verbose && println("\n[2/6] Generating Julia forcing function...")
    
    arg_names = length(test.variables) == 2 ? [:x] : [:x, :t]
    forcing_code = generate_forcing_code(forcing, :julia, arg_names)
    
    if verbose
        println("  Generated function:")
        for line in split(forcing_code, '\n')[1:min(3, length(split(forcing_code, '\n')))]
            println("    ", line)
        end
        println("    ...")
    end
    
    # Step 3: Inject into source file
    verbose && println("\n[3/6] Injecting forcing into source file...")
    
    modified_file = replace(test.source_file, r"\.jl$" => "_mms_test.jl")
    marker_found = inject_forcing_code(
        test.source_file,
        forcing_code,
        modified_file,
        marker = test.marker
    )
    
    if !marker_found
        warning_msg = "Warning: Marker '$(test.marker)' not found in $(test.source_file)"
        verbose && println("  ⚠ $warning_msg")
    else
        verbose && println("  ✓ Code injected successfully")
    end
    
    verbose && println("  Modified file: $modified_file")
    
    # Step 4: Run the modified code
    verbose && println("\n[4/6] Running solver...")
    
    try
        # Set environment variable for resolution
        ENV["MMS_TEST_RESOLUTION"] = string(test.resolution)
        ENV["MMS_TEST_OUTPUT"] = test.output_file
        
        # Run the modified solver
        cmd = `julia --project=. $modified_file`
        output = read(cmd, String)
        
        if verbose
            println("  Solver output:")
            for line in split(output, '\n')[1:min(5, length(split(output, '\n')))]
                println("    ", line)
            end
            if length(split(output, '\n')) > 5
                println("    ...")
            end
        end
        
    catch e
        error_msg = "Failed to run solver: $e"
        verbose && println("  ✗ $error_msg")
        return Dict(
            :passed => false,
            :L2_error => NaN,
            :Linf_error => NaN,
            :message => error_msg
        )
    end
    
    # Step 5: Read output
    verbose && println("\n[5/6] Reading solver output...")
    
    if !isfile(test.output_file)
        error_msg = "Output file not found: $(test.output_file)"
        verbose && println("  ✗ $error_msg")
        return Dict(
            :passed => false,
            :L2_error => NaN,
            :Linf_error => NaN,
            :message => error_msg
        )
    end
    
    # Read the output (assuming format: x u_numerical)
    data = try
        readdlm(test.output_file, skipstart=1)
    catch e
        error_msg = "Failed to read output file: $e"
        verbose && println("  ✗ $error_msg")
        return Dict(
            :passed => false,
            :L2_error => NaN,
            :Linf_error => NaN,
            :message => error_msg
        )
    end
    
    x_values = data[:, 1]
    u_numerical = data[:, 2]
    
    verbose && println("  Read $(length(x_values)) solution points")
    
    # Step 6: Compute analytical solution and compare
    verbose && println("\n[6/6] Verifying against manufactured solution...")
    
    # Evaluate manufactured solution
    u_mms_func = if length(test.variables) == 2
        x_var = test.variables[1]
        build_function(test.manufactured_solution, x_var, expression=Val{false})
    else
        x_var, t_var = test.variables[1:2]
        # For time-dependent, we need final time - assume it's in ENV or default
        t_final = get(ENV, "MMS_TEST_FINAL_TIME", "1.0")
        t_final = parse(Float64, t_final)
        mms_with_t = substitute(test.manufactured_solution, t_var => t_final)
        build_function(mms_with_t, x_var, expression=Val{false})
    end
    
    u_analytical = [u_mms_func(x) for x in x_values]
    
    # Compute errors
    h = x_values[2] - x_values[1]
    grid_dummy = Grid1D(test.domain..., h)
    
    L2_error = compute_L2_error(u_numerical, u_analytical, grid_dummy)
    Linf_error = compute_Linf_error(u_numerical, u_analytical)
    
    # Determine pass/fail
    passed = L2_error < test.tolerance
    
    if verbose
        @printf "  L2 error:   %.6e\n" L2_error
        @printf "  L∞ error:   %.6e\n" Linf_error
        @printf "  Tolerance:  %.6e\n" test.tolerance
        println()
        if passed
            println("  ✓ TEST PASSED")
        else
            println("  ✗ TEST FAILED")
        end
    end
    
    message = if passed
        "L2 error ($(round(L2_error, sigdigits=4))) < tolerance ($(test.tolerance))"
    else
        "L2 error ($(round(L2_error, sigdigits=4))) > tolerance ($(test.tolerance))"
    end
    
    return Dict(
        :passed => passed,
        :L2_error => L2_error,
        :Linf_error => Linf_error,
        :message => message,
        :num_points => length(x_values),
        :resolution => h
    )
end

"""
    @mms_test name solver_file pde u_mms variables domain

Macro for convenient MMS test definition.

# Example
```julia
@mms_test "My Poisson Test" "solver.jl" begin
    @variables x u(x)
    Dx = Differential(x)
    pde = -Dx(Dx(u))
    u_mms = sin(π*x)
    variables = (x, u)
    domain = (0.0, 1.0)
end
```
"""
macro mms_test(name, solver_file, block)
    quote
        # Execute the block to get variables
        $(esc(block))
        
        # Create and run test
        test = MMSUnitTest(
            $(esc(name)),
            $(esc(solver_file)),
            $(esc(:pde)),
            $(esc(:u_mms)),
            $(esc(:variables)),
            $(esc(:domain))
        )
        
        result = run_mms_unit_test(test)
        
        # Integration with Test framework
        if isdefined(Main, :Test)
            Main.Test.@test result[:passed]
        end
        
        result
    end
end

export MMSUnitTest, run_mms_unit_test, @mms_test
