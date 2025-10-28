"""
Example: Heat Equation with Automatic Symbolic Forcing

This example shows how to use the symbolic interface for time-dependent PDEs.
The forcing term is computed automatically from the PDE and manufactured solution.

Problem: u_t = κ*u_xx + f(x,t) on [0, 1] × [0, T] with periodic BCs
"""

using ManufacturedSolutions
using Printf

# ============================================================================
# Step 1: Define the PDE symbolically
# ============================================================================

@variables x t u(x, t)
Dx = Differential(x)
Dt = Differential(t)

# Heat equation: u_t - κ*u_xx = f(x,t)
κ = 0.1
pde_lhs = Dt(u) - κ * Dx(Dx(u))

println("="^70)
println("Symbolic Heat Equation")
println("="^70)
println("PDE: u_t = κ*u_xx + f(x,t)")
println("Symbolic form: ", pde_lhs, " = f(x,t)")
println("Diffusivity κ = ", κ)
println()

# ============================================================================
# Step 2: Define manufactured solution symbolically
# ============================================================================

# Choose: u(x,t) = exp(-t) * sin(2πx)
u_manufactured_symbolic = exp(-t) * sin(2π * x)

println("Manufactured Solution: u(x,t) = ", u_manufactured_symbolic)
println()

# ============================================================================
# Step 3: Automatically compute forcing term
# ============================================================================

println("Computing forcing term symbolically...")
forcing_term = compute_forcing_term_time_dependent(pde_lhs, u_manufactured_symbolic, x, t, u)

println("Forcing term f(x,t) = ", forcing_term.symbolic_form)
println()

# Create callable manufactured solution
u_mms_func = build_function(u_manufactured_symbolic, [x, t], expression=Val{false})

# ============================================================================
# Step 4: User's heat equation solver (same as example 02)
# ============================================================================

function heat_solver_1d_periodic(grid::Grid1D, forcing_func, u_initial; 
                                  κ=0.1, t_final=0.1, cfl=0.25)
    nx = grid.nx
    dx = grid.dx
    x_vals = grid.x
    
    # Time step (CFL condition)
    dt = cfl * dx^2 / κ
    nt = ceil(Int, t_final / dt)
    dt = t_final / nt
    
    println("  Time integration: nt = $nt, dt = $dt")
    
    # Initialize
    u = u_initial.(x_vals)
    u_new = similar(u)
    
    # Time stepping
    for n in 1:nt
        t_current = n * dt
        
        for i in 1:nx
            # Periodic indices
            i_left = (i == 1) ? nx : i-1
            i_right = (i == nx) ? 1 : i+1
            
            # Centered difference for u_xx
            u_xx = (u[i_left] - 2*u[i] + u[i_right]) / dx^2
            
            # Forward Euler: u^(n+1) = u^n + dt*(κ*u_xx + f)
            u_new[i] = u[i] + dt * (κ * u_xx + forcing_func(x_vals[i], t_current))
        end
        
        u, u_new = u_new, u
    end
    
    return u
end

# ============================================================================
# Step 5: Create solver wrapper
# ============================================================================

t_final = 0.1

# Initial condition from manufactured solution at t=0
function u_initial(x_val)
    return u_mms_func([x_val, 0.0])
end

function solver_wrapper(grid::Grid1D, mms; kwargs...)
    return heat_solver_1d_periodic(
        grid, 
        forcing_term,  # Automatically computed!
        u_initial,
        κ = κ,
        t_final = t_final,
        cfl = 0.25
    )
end

# Manufactured solution at final time
struct FinalTimeSolution
    func::Function
    t_final::Float64
end

(fts::FinalTimeSolution)(x_val) = fts.func([x_val, fts.t_final])

u_mms = ManufacturedSolution(
    FinalTimeSolution(u_mms_func, t_final),
    (0.0, 1.0),
    name = "exp(-t)*sin(2πx) at t=$t_final",
    description = "Time-dependent solution with automatic forcing"
)

# ============================================================================
# Step 6: Run verification
# ============================================================================

println("="^70)
println("Running Verification with Symbolic Time-Dependent Forcing")
println("="^70)
println()

solver = UserSolver(
    solver_wrapper,
    name = "1D Heat Solver (Symbolic Forcing)",
    description = "Automatically computed forcing for u_t = κu_xx + f"
)

results = verify_solver(
    solver,
    pde_lhs,
    u_mms,
    resolutions = [0.05, 0.025, 0.0125, 0.00625],
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

println("\nConvergence Rates:")
println("-"^50)
println("Expected: ~1.0 (Forward Euler is 1st order in time)")
println()

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
if l2_rate > 0.9 && l2_rate < 1.1
    println("✓ PASSED: L2 convergence rate ≈ 1.0 (1st order time integration)")
    println("\n🎉 Symbolic time-dependent forcing term computation successful!")
else
    println("⚠ Note: L2 convergence rate = $(round(l2_rate, digits=3))")
end

println()
println("="^70)
println("\nKey Achievement: The forcing term")
println("  f(x,t) = ", forcing_term.symbolic_form)
println("\nwas computed AUTOMATICALLY by taking derivatives of")
println("  u(x,t) = ", u_manufactured_symbolic)
println("\nand substituting into")
println("  ", pde_lhs, " = f(x,t)")
println("="^70)
