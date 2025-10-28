"""
Simple unit testing interface using function wrappers (Julia only)

This provides a clean, straightforward approach for testing Julia solvers
without any code injection - just pass your solver function.
"""

"""
    SimpleMMS Test

Simple MMS test configuration for Julia solvers.
"""
struct SimpleMMSTest
    test_name::String
    solver_function::Function        # User's solver: (grid, forcing) -> solution
    pde::Any                        # Symbolic PDE
    manufactured_solution::Any      # Symbolic manufactured solution
    variables::Tuple                # (x, u) or (x, t, u)
    domain::Tuple{Float64,Float64} # Spatial domain
    boundary_conditions::Symbol     # :dirichlet, :periodic, :neumann
    resolution::Float64             # Grid spacing
    tolerance::Float64              # Error tolerance
    
    function SimpleMMSTest(
        test_name::String,
        solver_function::Function,
        pde,
        manufactured_solution,
        variables::Tuple,
        domain::Tuple;
        boundary_conditions::Symbol = :dirichlet,
        resolution::Float64 = 0.01,
        tolerance::Float64 = 1e-3
    )
        new(test_name, solver_function, pde, manufactured_solution,
            variables, domain, boundary_conditions, resolution, tolerance)
    end
end

"""
    run_simple_mms_test(test::SimpleMMSTest; verbose::Bool=true)

Run a simple MMS unit test using function wrapper approach.

# Workflow
1. Compute forcing term symbolically
2. Create grid
3. Evaluate forcing and manufactured solution
4. Call user's solver
5. Compute errors
6. Return pass/fail

# Example
```julia
@variables x u(x)
Dx = Differential(x)

function my_solver(grid, forcing)
    # ... your implementation ...
    return u_numerical
end

test = SimpleMMSTest(
    "My Test",
    my_solver,
    -Dx(Dx(u)),
    sin(π*x),
    (x, u),
    (0.0, 1.0)
)

result = run_simple_mms_test(test)
```
"""
function run_simple_mms_test(test::SimpleMMSTest; verbose::Bool=true)
    verbose && println("="^70)
    verbose && println("MMS Test: $(test.test_name)")
    verbose && println("="^70)
    
    # Step 1: Compute forcing term symbolically
    verbose && println("\n[1/5] Computing forcing term symbolically...")
    
    forcing_term = if length(test.variables) == 2
        # Spatial only
        x_var, u_var = test.variables
        compute_forcing_term_1d(test.pde, test.manufactured_solution, x_var, u_var)
    else
        # Time-dependent
        x_var, t_var, u_var = test.variables
        compute_forcing_term_time_dependent(test.pde, test.manufactured_solution,
                                          x_var, t_var, u_var)
    end
    
    verbose && println("  Forcing: ", forcing_term.symbolic_form)
    
    # Step 2: Create grid
    verbose && println("\n[2/5] Creating grid...")
    grid = Grid1D(test.domain..., test.resolution)
    verbose && println("  Domain: [$(test.domain[1]), $(test.domain[2])]")
    verbose && println("  Grid points: $(grid.nx)")
    verbose && println("  Spacing: h = $(grid.dx)")
    
    # Step 3: Evaluate manufactured solution
    verbose && println("\n[3/5] Evaluating manufactured solution...")
    
    x_var = test.variables[1]
    u_mms_func = build_function(test.manufactured_solution, x_var, expression=Val{false})
    u_analytical = [u_mms_func(xi) for xi in grid.x]
    
    verbose && println("  Solution range: [$(minimum(u_analytical)), $(maximum(u_analytical))]")
    
    # Step 4: Run solver
    verbose && println("\n[4/5] Running solver...")
    
    try
        u_numerical = test.solver_function(grid, forcing_term)
        
        if length(u_numerical) != grid.nx
            error_msg = "Solver returned $(length(u_numerical)) points, expected $(grid.nx)"
            verbose && println("  ✗ $error_msg")
            return Dict(
                :passed => false,
                :L2_error => NaN,
                :Linf_error => NaN,
                :message => error_msg
            )
        end
        
        verbose && println("  Numerical solution range: [$(minimum(u_numerical)), $(maximum(u_numerical))]")
        
        # Step 5: Compute errors
        verbose && println("\n[5/5] Computing errors...")
        
        L2_error = compute_L2_error(u_numerical, u_analytical, grid)
        Linf_error = compute_Linf_error(u_numerical, u_analytical)
        L1_error = compute_L1_error(u_numerical, u_analytical, grid)
        
        passed = L2_error < test.tolerance
        
        if verbose
            @printf "  L1 error:   %.6e\n" L1_error
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
            :L1_error => L1_error,
            :L2_error => L2_error,
            :Linf_error => Linf_error,
            :message => message,
            :num_points => grid.nx,
            :resolution => grid.dx,
            :u_numerical => u_numerical,
            :u_analytical => u_analytical,
            :grid => grid
        )
        
    catch e
        error_msg = "Solver failed: $e"
        verbose && println("  ✗ $error_msg")
        if verbose
            println("\nStacktrace:")
            showerror(stdout, e, catch_backtrace())
            println()
        end
        return Dict(
            :passed => false,
            :L2_error => NaN,
            :Linf_error => NaN,
            :message => error_msg
        )
    end
end

"""
    @test_solver solver_func pde u_mms variables domain [options]

Macro for convenient test creation and execution.
"""
macro test_solver(name, solver, pde_expr, u_mms_expr, vars, dom, options...)
    quote
        test = SimpleMMSTest(
            $(esc(name)),
            $(esc(solver)),
            $(esc(pde_expr)),
            $(esc(u_mms_expr)),
            $(esc(vars)),
            $(esc(dom)),
            $(esc(options))...
        )
        run_simple_mms_test(test)
    end
end

export SimpleMMSTest, run_simple_mms_test, @test_solver
