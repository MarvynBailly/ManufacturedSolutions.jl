"""
Example: Verification of a simple 1D Poisson solver

This example demonstrates how to use ManufacturedSolutions.jl to verify
a user's Poisson solver using the Method of Manufactured Solutions.

Problem: -u'' = f(x) on [0, 1] with u(0) = u(1) = 0
"""

using ManufacturedSolutions
using Printf

# ============================================================================
# User's Poisson Solver (the code we want to verify)
# ============================================================================

"""
Simple 1D Poisson solver using finite differences.

Solves: -u'' = f(x) on domain [0, 1] with Dirichlet BCs u(0) = u(1) = 0

# Arguments
- `grid::Grid1D`: Spatial grid
- `forcing`: Function f(x) representing the forcing term
- `u_left::Float64`: Left boundary value (default: 0.0)
- `u_right::Float64`: Right boundary value (default: 0.0)

# Returns
- `u::Vector`: Numerical solution
"""
function poisson_solver_1d(grid::Grid1D, forcing; u_left=0.0, u_right=0.0)
    nx = grid.nx
    dx = grid.dx
    x = grid.x
    
    # Set up the linear system: A*u = b
    # Using centered finite differences: -u''(x) ≈ -(u[i-1] - 2u[i] + u[i+1])/dx²
    
    # Create tridiagonal matrix A
    A = zeros(nx, nx)
    b = zeros(nx)
    
    # Boundary conditions (first and last rows)
    A[1, 1] = 1.0
    b[1] = u_left
    
    A[nx, nx] = 1.0
    b[nx] = u_right
    
    # Interior points: -u'' = f
    # -(u[i-1] - 2u[i] + u[i+1])/dx² = f[i]
    # Rearrange: -u[i-1] + 2u[i] - u[i+1] = dx² * f[i]
    for i in 2:(nx-1)
        A[i, i-1] = -1.0
        A[i, i]   =  2.0
        A[i, i+1] = -1.0
        b[i] = dx^2 * forcing(x[i])
    end
    
    # Solve the linear system
    u = A \ b
    
    return u
end

# ============================================================================
# Manufactured Solution Setup
# ============================================================================

"""
Manufactured solution: u(x) = sin(π*x)

This satisfies:
- u(0) = sin(0) = 0 ✓
- u(1) = sin(π) = 0 ✓
- u''(x) = -π²*sin(π*x)

So the PDE becomes:
-u''(x) = f(x)
-(-π²*sin(π*x)) = f(x)
f(x) = π²*sin(π*x)
"""

# Define the manufactured solution
u_manufactured = ManufacturedSolution(
    x -> sin(π * x),
    (0.0, 1.0),
    name = "sin(πx)",
    description = "Smooth solution satisfying homogeneous Dirichlet BCs"
)

# Define the forcing term (computed analytically for this example)
# In general, this would be computed symbolically
forcing_function(x) = π^2 * sin(π * x)

# ============================================================================
# Create a wrapper that matches our interface
# ============================================================================

"""
Wrapper for our Poisson solver that matches the ManufacturedSolutions interface.

The framework expects: solver(grid, manufactured_solution, ...) -> numerical_solution
"""
function poisson_solver_wrapper(grid::Grid1D, manufactured_solution::ManufacturedSolution; kwargs...)
    # Use the known forcing function
    # (In future, this will be computed automatically from the PDE and manufactured solution)
    return poisson_solver_1d(grid, forcing_function; kwargs...)
end

# ============================================================================
# Run Verification
# ============================================================================

println("="^70)
println("1D Poisson Solver Verification using MMS")
println("="^70)
println()
println("Problem: -u''(x) = f(x) on [0,1] with u(0) = u(1) = 0")
println("Manufactured Solution: u(x) = sin(πx)")
println("Forcing Term: f(x) = π²sin(πx)")
println()

# Create the user solver
solver = UserSolver(
    poisson_solver_wrapper,
    name = "1D Finite Difference Poisson Solver",
    description = "2nd order centered differences"
)

# Run verification with multiple grid resolutions
println("Running convergence study...")
println()

results = verify_solver(
    solver,
    nothing,  # PDE (will be used when symbolic processing is implemented)
    u_manufactured,
    resolutions = [0.1, 0.05, 0.025, 0.0125],
    verbose = true
)

println()
println("="^70)
println("Verification Summary")
println("="^70)

# Display convergence rates
println("\nConvergence Rates (Expected: 2.0 for 2nd order method):")
println("-"^50)
for (norm, rate) in results.convergence_rates
    if !endswith(string(norm), "_R2")
        R2_key = Symbol(string(norm) * "_R2")
        R2 = results.convergence_rates[R2_key]
        @printf "%-10s  Rate: %.3f  (R² = %.4f)\n" norm rate R2
    end
end

# Check if convergence is as expected
println("\nVerification Assessment:")
println("-"^50)
l2_rate = results.convergence_rates[:L2]
if abs(l2_rate - 2.0) < 0.1
    println("✓ PASSED: L2 convergence rate ≈ 2.0 (2nd order method)")
else
    println("✗ WARNING: L2 convergence rate = $(round(l2_rate, digits=3)), expected ≈ 2.0")
end

println()
println("="^70)
