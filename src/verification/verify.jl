"""
Main verification workflow for user solvers
"""

using Printf

"""
    verify_solver(
        user_solver::UserSolver,
        pde,
        manufactured_solution;
        resolutions::Vector = [0.1, 0.05, 0.025],
        error_norms::Vector{Symbol} = [:L1, :L2, :Linf],
        verbose::Bool = true,
        solver_kwargs...
    )

Verify a user's numerical solver using the Method of Manufactured Solutions.

# Arguments
- `user_solver::UserSolver`: User's solver wrapped in UserSolver
- `pde`: Symbolic PDE expression (optional - can be `nothing` if forcing is provided manually)
- `manufactured_solution`: ManufacturedSolution or function
- `resolutions::Vector`: Grid resolutions to test (default: [0.1, 0.05, 0.025])
- `error_norms::Vector{Symbol}`: Error norms to compute (default: [:L1, :L2, :Linf])
- `verbose::Bool`: Print progress information (default: true)
- `solver_kwargs...`: Additional arguments passed to user's solver

# Returns
- `VerificationResult`: Complete verification results

# Example
```julia
# With symbolic forcing (automatic)
@variables x u(x)
Dx = Differential(x)
pde_lhs = -Dx(Dx(u))
u_mms_sym = sin(π * x)

# Or with manual forcing
function my_solver_wrapper(grid, forcing; kwargs...)
    # forcing is automatically computed if pde is provided
    return numerical_solution
end

solver = UserSolver(my_solver_wrapper, name="My Solver")
results = verify_solver(solver, pde_lhs, u_mms, resolutions=[0.1, 0.05])
```
"""
function verify_solver(
    user_solver::UserSolver,
    pde,
    manufactured_solution;
    resolutions::Vector = [0.1, 0.05, 0.025],
    error_norms::Vector{Symbol} = [:L1, :L2, :Linf],
    verbose::Bool = true,
    solver_kwargs...
)
    verbose && @info "Starting solver verification" solver=user_solver.metadata[:name]
    
    # Store convergence data for each resolution
    convergence_data = ConvergenceData{Float64}[]
    
    for (idx, h) in enumerate(resolutions)
        verbose && @info "Resolution $idx/$(length(resolutions)): h = $h"
        
        # Create grid
        domain = manufactured_solution.domain
        grid = Grid1D(domain..., h)
        
        verbose && @info "  Grid: $(grid.nx) points, dx = $(grid.dx)"
        
        # Time the solver
        t_start = time()
        
        # Call user's solver with forcing term
        # For now, we'll need to pass the forcing function
        # This will be improved when we implement symbolic forcing computation
        numerical_solution = user_solver.solve_function(grid, manufactured_solution; solver_kwargs...)
        
        t_elapsed = time() - t_start
        
        verbose && @info "  Solver completed in $(round(t_elapsed, digits=4)) seconds"
        
        # Evaluate analytical solution on grid
        analytical_solution = manufactured_solution.(grid.x)
        
        # Compute errors
        errors = compute_all_errors(numerical_solution, analytical_solution, grid)
        
        if verbose
            for (norm_name, error_value) in errors
                @printf "  %s error: %.6e\n" norm_name error_value
            end
        end
        
        # Store data
        push!(convergence_data, ConvergenceData(
            h,
            grid.nx,
            errors,
            t_elapsed
        ))
    end
    
    # Analyze convergence
    verbose && @info "Analyzing convergence rates..."
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
    
    # Create result object
    problem_info = Dict{Symbol,Any}(
        :domain => manufactured_solution.domain,
        :manufactured_solution => get(manufactured_solution.metadata, :name, "Unknown")
    )
    
    result = VerificationResult(
        problem_info,
        user_solver.metadata,
        convergence_data,
        convergence_rates
    )
    
    verbose && @info "Verification complete!"
    
    return result
end

"""
    verify_solver(
        solve_function::Function,
        manufactured_solution::ManufacturedSolution;
        kwargs...
    )

Convenience method that automatically wraps the solver function.
"""
function verify_solver(
    solve_function::Function,
    manufactured_solution::ManufacturedSolution;
    solver_name::String = "User Solver",
    kwargs...
)
    solver = UserSolver(solve_function, name=solver_name)
    return verify_solver(solver, nothing, manufactured_solution; kwargs...)
end

export verify_solver
