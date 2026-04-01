"""
Simple unit testing interface for verifying numerical solvers

This module provides single-resolution verification tests that check if a solver
produces solutions within a specified error tolerance. Unlike convergence tests
(which verify order of accuracy), these tests simply verify that the solver
produces the correct solution to within acceptable error bounds.

# When to use this vs. ConvergenceTest
- Use SimpleMMSTest: Quick verification that solver works correctly
- Use ConvergenceTest: Verify that discretization order is correct

Both are valuable! Simple tests catch bugs faster, convergence tests verify accuracy.
"""

"""
    SimpleMMSTest

Configuration for a single-resolution MMS verification test.

This test runs your solver once at a specified resolution and checks if the
error is below a tolerance. Think of it as a unit test for your solver.

# Example
```julia
test = SimpleMMSTest(
    "Poisson solver test",
    my_solver,
    -Dx(Dx(u)),          # PDE: -u'' = f
    sin(π*x),            # Manufactured solution
    (x, u),              # Variables
    (0.0, 1.0),          # Domain
    resolution = 0.01,   # Grid spacing
    tolerance = 1e-3     # Error must be < 0.001
)
result = run_simple_mms_test(test)
@assert result[:passed]  # Test passed!
```
"""
struct SimpleMMSTest
    test_name::String
    solver_function::Function        # Your solver: (grid, forcing) -> numerical_solution
    pde::Any                        # Symbolic PDE (e.g., -Dx(Dx(u)) for Poisson)
    manufactured_solution::Any      # Symbolic solution (e.g., sin(π*x))
    variables::Tuple                # Symbolic variables: (x, u) or (x, t, u)
    domain::Tuple{Float64,Float64} # Spatial domain (a, b)
    boundary_conditions::Symbol     # :dirichlet, :periodic, or :neumann
    resolution::Float64             # Grid spacing h
    tolerance::Float64              # Maximum acceptable L2 error
    
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

Run a single-resolution MMS verification test.

This function performs the complete MMS workflow:
1. Computes the forcing term needed to produce the manufactured solution
2. Creates a grid at the specified resolution
3. Runs your solver with the computed forcing
4. Compares numerical solution to analytical (manufactured) solution
5. Reports pass/fail based on whether error < tolerance

# Arguments
- `test::SimpleMMSTest`: Test configuration
- `verbose::Bool=true`: Print detailed progress information

# Returns
Dictionary with keys:
- `:passed`: Bool - did the test pass?
- `:L1_error`, `:L2_error`, `:Linf_error`: Error norms
- `:message`: Human-readable result message
- `:num_points`: Number of grid points used
- `:resolution`: Grid spacing used
- `:u_numerical`, `:u_analytical`, `:grid`: Detailed results

# Example
```julia
result = run_simple_mms_test(test)
if result[:passed]
    println("✓ Solver works correctly!")
    println("  L2 error: ", result[:L2_error])
else
    println("✗ Solver failed")
    println("  ", result[:message])
end
```
"""
function run_simple_mms_test(test::SimpleMMSTest; verbose::Bool=true)
    verbose && println("="^70)
    verbose && println("MMS Test: $(test.test_name)")
    verbose && println("="^70)
    
    # ==========================================================================
    # STEP 1: Compute forcing term symbolically
    # ==========================================================================
    # Given a PDE and manufactured solution, we need to compute what forcing
    # term would produce that solution. For example:
    #   PDE: -u'' = f
    #   Manufactured: u = sin(πx)
    #   Then: u'' = -π²sin(πx), so we need f = π²sin(πx)
    
    verbose && println("\n[1/5] Computing forcing term symbolically...")
    
    forcing_term = if length(test.variables) == 2
        # Spatial only (e.g., steady-state Poisson equation)
        x_var, u_var = test.variables
        compute_forcing_term_1d(test.pde, test.manufactured_solution, x_var, u_var)
    else
        # Time-dependent (e.g., heat equation, wave equation)
        x_var, t_var, u_var = test.variables
        compute_forcing_term_time_dependent(test.pde, test.manufactured_solution,
                                          x_var, t_var, u_var)
    end
    
    verbose && println("  Forcing: ", forcing_term.symbolic_form)
    
    # ==========================================================================
    # STEP 2: Create computational grid
    # ==========================================================================
    # Create a uniform 1D grid with the specified resolution
    
    verbose && println("\n[2/5] Creating grid...")
    grid = Grid1D(test.domain..., test.resolution)
    verbose && println("  Domain: [$(test.domain[1]), $(test.domain[2])]")
    verbose && println("  Grid points: $(grid.nx)")
    verbose && println("  Spacing: h = $(grid.dx)")
    
    # ==========================================================================
    # STEP 3: Evaluate analytical (manufactured) solution on the grid
    # ==========================================================================
    # We need to know the exact solution at each grid point so we can compare
    # it to the numerical solution later
    
    verbose && println("\n[3/5] Evaluating manufactured solution...")
    
    x_var = test.variables[1]
    u_mms_func = build_function(test.manufactured_solution, x_var, expression=Val{false})
    u_analytical = [u_mms_func(xi) for xi in grid.x]
    
    verbose && println("  Solution range: [$(minimum(u_analytical)), $(maximum(u_analytical))]")
    
    # ==========================================================================
    # STEP 4: Run the numerical solver
    # ==========================================================================
    # Call the user's solver with the grid and computed forcing term.
    # The solver should return a vector of solution values at grid points.
    
    verbose && println("\n[4/5] Running solver...")
    
    try
        u_numerical = test.solver_function(grid, forcing_term)
        
        # Sanity check: solver must return the correct number of points
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
        
        # ==========================================================================
        # STEP 5: Compute errors and check tolerance
        # ==========================================================================
        # Compare numerical solution to analytical solution using standard error norms:
        #   - L1: Average absolute error
        #   - L2: Root-mean-square error (most common)
        #   - L∞: Maximum absolute error
        
        verbose && println("\n[5/5] Computing errors...")
        
        L2_error = compute_L2_error(u_numerical, u_analytical, grid)
        Linf_error = compute_Linf_error(u_numerical, u_analytical)
        L1_error = compute_L1_error(u_numerical, u_analytical, grid)
        
        # Test passes if L2 error is below tolerance
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
