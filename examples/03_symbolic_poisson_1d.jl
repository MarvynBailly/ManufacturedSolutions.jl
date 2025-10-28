"""
Example: Automatic Forcing Term Computation using Symbolic Interface

This example demonstrates the symbolic interface that automatically computes
forcing terms from the PDE and manufactured solution.

Problem: -u''(x) = f(x) on [0, 1] with u(0) = u(1) = 0
"""

using ManufacturedSolutions
using Printf

# ============================================================================
# Step 1: Define the PDE symbolically
# ============================================================================

@variables x u(x)
Dx = Differential(x)

# PDE: -u'' = f(x)
# Left-hand side of the PDE (what operates on u)
pde_lhs = -Dx(Dx(u))

println("="^70)
println("Symbolic PDE Definition")
println("="^70)
println("PDE: -u''(x) = f(x)")
println("Symbolic form: ", pde_lhs, " = f(x)")
println()

# ============================================================================
# Step 2: Define manufactured solution symbolically
# ============================================================================

# Choose: u(x) = sin(πx)
u_manufactured_symbolic = sin(π * x)

println("Manufactured Solution: u(x) = ", u_manufactured_symbolic)
println()

# ============================================================================
# Step 3: Automatically compute forcing term
# ============================================================================

println("Computing forcing term symbolically...")
forcing_term = compute_forcing_term_1d(pde_lhs, u_manufactured_symbolic, x, u)

println("Forcing term f(x) = ", forcing_term.symbolic_form)
println()

# ============================================================================
# Step 4: Create callable manufactured solution
# ============================================================================

# Build function from symbolic expression
u_mms_func = build_function(u_manufactured_symbolic, x, expression=Val{false})

# Wrap in ManufacturedSolution type
u_mms = ManufacturedSolution(
    u_mms_func,
    (0.0, 1.0),
    name = "sin(πx)",
    description = "Automatically derived from symbolic expression"
)

# ============================================================================
# Step 5: User's Poisson solver (same as before)
# ============================================================================

function poisson_solver_1d(grid::Grid1D, forcing; u_left=0.0, u_right=0.0)
    nx = grid.nx
    dx = grid.dx
    x = grid.x
    
    # Set up the linear system: A*u = b
    A = zeros(nx, nx)
    b = zeros(nx)
    
    # Boundary conditions
    A[1, 1] = 1.0
    b[1] = u_left
    
    A[nx, nx] = 1.0
    b[nx] = u_right
    
    # Interior points: -u'' = f
    for i in 2:(nx-1)
        A[i, i-1] = -1.0
        A[i, i]   =  2.0
        A[i, i+1] = -1.0
        b[i] = dx^2 * forcing(x[i])
    end
    
    # Solve
    u = A \ b
    return u
end

# ============================================================================
# Step 6: Create solver wrapper that uses automatic forcing
# ============================================================================

function solver_wrapper(grid::Grid1D, mms::ManufacturedSolution; kwargs...)
    # The forcing term is passed via the forcing_term variable
    # (We need to access it from the outer scope for now)
    return poisson_solver_1d(grid, forcing_term; kwargs...)
end

# ============================================================================
# Step 7: Run verification
# ============================================================================

println("="^70)
println("Running Verification with Automatic Forcing")
println("="^70)
println()

solver = UserSolver(
    solver_wrapper,
    name = "1D Poisson Solver (Symbolic Forcing)",
    description = "Using automatically computed forcing term"
)

results = verify_solver(
    solver,
    pde_lhs,
    u_mms,
    resolutions = [0.1, 0.05, 0.025, 0.0125],
    verbose = true
)

println()
println("="^70)
println("Verification Summary")
println("="^70)

println("\nSymbolic Derivation:")
println("-"^50)
println("PDE:              ", pde_lhs, " = f")
println("Manufactured Sol: ", u_manufactured_symbolic)
println("Forcing Term:     ", forcing_term.symbolic_form)

println("\nConvergence Rates (Expected: 2.0 for 2nd order method):")
println("-"^50)
for (norm, rate) in results.convergence_rates
    if !endswith(string(norm), "_R2")
        R2_key = Symbol(string(norm) * "_R2")
        R2 = results.convergence_rates[R2_key]
        @printf "%-10s  Rate: %.3f  (R² = %.4f)\n" norm rate R2
    end
end

println("\nVerification Assessment:")
println("-"^50)
l2_rate = results.convergence_rates[:L2]
if abs(l2_rate - 2.0) < 0.1
    println("✓ PASSED: L2 convergence rate ≈ 2.0 (2nd order method)")
    println("\n🎉 Symbolic forcing term computation successful!")
else
    println("✗ WARNING: L2 convergence rate = $(round(l2_rate, digits=3)), expected ≈ 2.0")
end

println()
println("="^70)
println("\nKey Takeaway: The forcing term f(x) = $(forcing_term.symbolic_form)")
println("was computed AUTOMATICALLY from the PDE and manufactured solution!")
println("="^70)
