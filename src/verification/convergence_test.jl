"""
Convergence-based Verification using Method of Manufactured Solutions

This module provides **order-of-accuracy verification**, which is more rigorous
than simple error tolerance testing. Instead of checking if error < tolerance
at one resolution, we verify that error decreases at the correct rate as we
refine the grid.

## Why Convergence Testing?

A second-order accurate method should satisfy:
    error(h) ≈ C × h²

This means:
- Halving grid spacing (h → h/2) should quarter the error
- Going h → h/4 should reduce error by factor of 16
- The "convergence order" p in error ~ h^p should equal 2

By testing on multiple resolutions (h, h/2, h/4, h/8), we can compute the
observed convergence order using log-log regression and verify it matches
the method's theoretical order.

## How It Works

1. **Setup**: Choose manufactured solution u(x), compute forcing f(x)
2. **Solve**: Run solver on progressively finer grids (h, h/2, h/4, ...)
3. **Measure**: Compute L2 error at each resolution
4. **Analyze**: Fit log(error) vs log(h) to find slope = convergence order
5. **Verify**: Check observed order ≈ expected order within tolerance

If observed order differs significantly from expected, the solver has a bug.

## Example
```julia
using Symbolics

# Define symbolic problem
@variables x u
pde = -Differential(x)^2(u)  # -u''(x) = f(x)
u_exact = sin(π*x)           # Manufactured solution

# Create convergence test for second-order solver
test = ConvergenceTest(
    "Poisson 1D - Sine Wave",
    my_solver,              # Your solver function
    pde,
    u_exact,
    (x, u),
    (0.0, 1.0),            # Domain
    expected_order = 2.0,  # Method is 2nd-order accurate
    base_resolution = 0.1,
    num_refinements = 4    # Test at h, h/2, h/4, h/8
)

# Run verification
result = run_convergence_test(test)

if result[:passed]
    println("✓ Solver achieves expected order ", result[:observed_order])
else
    println("✗ Solver order ", result[:observed_order], 
            " differs from expected ", result[:expected_order])
end
```

## When to Use This vs SimpleMMSTest

- **ConvergenceTest**: When you want to verify convergence order (most rigorous)
- **SimpleMMSTest**: When you just want quick pass/fail at fixed resolution

Convergence testing is the gold standard for verification because it checks
the asymptotic behavior, not just error magnitude at one resolution.
"""

using LinearAlgebra
using Printf

"""
    ConvergenceTest

MMS test that verifies the order of accuracy by running on multiple resolutions.

# Fields
- `test_name::String`: Descriptive name for the test
- `solver_function::Function`: User's solver (grid, forcing) -> solution
- `pde`: Symbolic PDE expression
- `manufactured_solution`: Symbolic manufactured solution
- `variables::Tuple`: (x, u) symbolic variables
- `domain::Tuple{Float64,Float64}`: Spatial domain (a, b)
- `expected_order::Float64`: Expected convergence order (e.g., 2.0 for second-order)
- `base_resolution::Float64`: Coarsest grid spacing
- `num_refinements::Int`: Number of grid refinements
- `tolerance::Float64`: Tolerance for order verification (default: 0.2)
"""
struct ConvergenceTest
    test_name::String
    solver_function::Function
    pde::Any
    manufactured_solution::Any
    variables::Tuple
    domain::Tuple{Float64,Float64}
    expected_order::Float64
    base_resolution::Float64
    num_refinements::Int
    tolerance::Float64
    
    function ConvergenceTest(
        test_name::String,
        solver_function::Function,
        pde,
        manufactured_solution,
        variables::Tuple,
        domain::Tuple;
        expected_order::Float64 = 2.0,
        base_resolution::Float64 = 0.1,
        num_refinements::Int = 4,
        tolerance::Float64 = 0.2
    )
        new(test_name, solver_function, pde, manufactured_solution,
            variables, domain, expected_order, base_resolution,
            num_refinements, tolerance)
    end
end

"""
    run_convergence_test(test::ConvergenceTest; verbose::Bool=true)

Run convergence test by computing errors on multiple grid resolutions.

# Algorithm
1. Compute forcing term symbolically
2. For each resolution (h, h/2, h/4, ...):
   - Create grid
   - Run solver
   - Compute L2 error
3. Compute observed convergence order using log-log regression
4. Compare observed order vs. expected order

# Returns
Dictionary with:
- `:passed`: Bool indicating if observed order matches expected
- `:expected_order`: Expected convergence order
- `:observed_order`: Measured convergence order
- `:resolutions`: Array of grid spacings tested
- `:errors`: Array of L2 errors
- `:order_tolerance`: Tolerance for order verification
"""
function run_convergence_test(test::ConvergenceTest; verbose::Bool=true)
    verbose && println("="^70)
    verbose && println("Convergence Test: $(test.test_name)")
    verbose && println("="^70)
    verbose && println("Domain: [$(test.domain[1]), $(test.domain[2])]")
    verbose && println("Expected Order: $(test.expected_order)")
    verbose && println("Testing $(test.num_refinements) refinement levels")
    verbose && println()
    
    # Step 1: Compute forcing term symbolically
    verbose && println("[1/3] Computing forcing term symbolically...")
    
    x_var, u_var = test.variables
    forcing_term = compute_forcing_term_1d(test.pde, test.manufactured_solution, x_var, u_var)
    
    verbose && println("  Forcing: ", forcing_term.symbolic_form)
    verbose && println()
    
    # Step 2: Run solver on multiple resolutions
    verbose && println("[2/3] Running solver on multiple resolutions...")
    verbose && println()
    
    resolutions = Float64[]
    errors = Float64[]
    num_points = Int[]
    
    for i in 0:(test.num_refinements - 1)
        h = test.base_resolution / (2.0^i)
        push!(resolutions, h)
        
        # Create grid
        grid = Grid1D(test.domain..., h)
        
        # Evaluate manufactured solution
        u_mms_func = build_function(test.manufactured_solution, x_var, expression=Val{false})
        u_analytical = [u_mms_func(xi) for xi in grid.x]
        
        # Run solver
        try
            u_numerical = test.solver_function(grid, forcing_term)
            
            if length(u_numerical) != grid.nx
                error("Solver returned $(length(u_numerical)) points, expected $(grid.nx)")
            end
            
            # Compute error
            L2_error = compute_L2_error(u_numerical, u_analytical, grid)
            push!(errors, L2_error)
            push!(num_points, grid.nx)
            
            if verbose
                @printf "  Level %d: h = %.6f, nx = %4d, L2 error = %.6e\n" i h grid.nx L2_error
            end
            
        catch e
            verbose && println("  ✗ Solver failed at refinement level $i: $e")
            return Dict(
                :passed => false,
                :message => "Solver failed at refinement level $i: $e",
                :expected_order => test.expected_order,
                :observed_order => NaN
            )
        end
    end
    
    verbose && println()
    
    # Step 3: Compute observed convergence order
    verbose && println("[3/3] Computing observed convergence order...")
    
    if length(errors) < 2
        verbose && println("  ✗ Not enough data points for convergence analysis")
        return Dict(
            :passed => false,
            :message => "Not enough refinement levels",
            :expected_order => test.expected_order,
            :observed_order => NaN
        )
    end
    
    # Use least-squares fit: log(error) = log(C) + p*log(h)
    # where p is the convergence order
    # Since error ~ C*h^p, we have log(error) = log(C) + p*log(h)
    # The slope of log(error) vs log(h) gives us p
    log_h = log.(resolutions)
    log_errors = log.(errors)
    
    # Linear regression: slope = p
    n = length(log_h)
    sum_x = sum(log_h)
    sum_y = sum(log_errors)
    sum_xy = sum(log_h .* log_errors)
    sum_xx = sum(log_h .^ 2)
    
    # Slope of log-log plot IS the convergence order
    observed_order = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x^2)
    
    # Compute error reduction ratios
    verbose && println()
    verbose && println("  Error Reduction Ratios:")
    for i in 2:length(errors)
        ratio = errors[i-1] / errors[i]
        theoretical_ratio = 2.0^test.expected_order
        @printf "    Level %d->%d: %.4f (expected: %.4f for order %.1f)\n" (i-2) (i-1) ratio theoretical_ratio test.expected_order
    end
    
    verbose && println()
    @printf "  Observed Order:  %.4f\n" observed_order
    @printf "  Expected Order:  %.4f\n" test.expected_order
    @printf "  Difference:      %.4f\n" abs(observed_order - test.expected_order)
    @printf "  Tolerance:       %.4f\n" test.tolerance
    verbose && println()
    
    # Check if observed order matches expected
    order_diff = abs(observed_order - test.expected_order)
    passed = order_diff < test.tolerance
    
    if verbose
        if passed
            println("  ✓ CONVERGENCE TEST PASSED")
            println("    Observed order $(round(observed_order, digits=3)) matches expected order $(test.expected_order)")
        else
            println("  ✗ CONVERGENCE TEST FAILED")
            println("    Observed order $(round(observed_order, digits=3)) differs from expected order $(test.expected_order)")
            println("    Difference $(round(order_diff, digits=3)) exceeds tolerance $(test.tolerance)")
        end
        println()
    end
    
    message = if passed
        "Observed order $(round(observed_order, digits=3)) ≈ expected order $(test.expected_order)"
    else
        "Observed order $(round(observed_order, digits=3)) ≠ expected order $(test.expected_order)"
    end
    
    return Dict(
        :passed => passed,
        :expected_order => test.expected_order,
        :observed_order => observed_order,
        :order_difference => order_diff,
        :order_tolerance => test.tolerance,
        :resolutions => resolutions,
        :errors => errors,
        :num_points => num_points,
        :message => message
    )
end

export ConvergenceTest, run_convergence_test
