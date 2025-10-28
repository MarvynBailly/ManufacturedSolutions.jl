"""
Example: Verification of a 1D Heat Equation Solver

This example demonstrates verification of a time-dependent PDE solver
using the Method of Manufactured Solutions.

Problem: u_t = κ*u_xx on [0, 1] × [0, T] with periodic boundary conditions
"""

using ManufacturedSolutions
using Printf

# ============================================================================
# User's Heat Equation Solver (the code we want to verify)
# ============================================================================

"""
Simple 1D heat equation solver using forward Euler in time and 
centered differences in space with periodic boundary conditions.

Solves: u_t = κ*u_xx on [0, 1] with periodic BCs

# Arguments
- `grid::Grid1D`: Spatial grid
- `forcing`: Function f(x,t) representing the forcing term (if any)
- `u_initial`: Initial condition function u(x,0)
- `κ::Float64`: Thermal diffusivity (default: 0.1)
- `t_final::Float64`: Final time (default: 0.1)
- `cfl::Float64`: CFL number for stability (default: 0.25)

# Returns
- `u::Vector`: Numerical solution at t_final
"""
function heat_solver_1d_periodic(grid::Grid1D, forcing, u_initial; 
                                  κ=0.1, t_final=0.1, cfl=0.25)
    nx = grid.nx
    dx = grid.dx
    x = grid.x
    
    # Time step (CFL condition for forward Euler: dt ≤ dx²/(2κ))
    dt = cfl * dx^2 / κ
    nt = ceil(Int, t_final / dt)
    dt = t_final / nt  # Adjust to hit t_final exactly
    
    println("  Time integration: nt = $nt, dt = $dt")
    
    # Initialize
    u = u_initial.(x)
    u_new = similar(u)
    
    # Time stepping
    for n in 1:nt
        t = n * dt
        
        # Interior points with periodic BCs
        for i in 1:nx
            # Periodic indices
            i_left = (i == 1) ? nx : i-1
            i_right = (i == nx) ? 1 : i+1
            
            # Centered difference for u_xx
            u_xx = (u[i_left] - 2*u[i] + u[i_right]) / dx^2
            
            # Forward Euler: u^(n+1) = u^n + dt*(κ*u_xx + f)
            u_new[i] = u[i] + dt * (κ * u_xx + forcing(x[i], t))
        end
        
        # Update
        u, u_new = u_new, u
    end
    
    return u
end

# ============================================================================
# Manufactured Solution Setup
# ============================================================================

"""
Manufactured solution: u(x,t) = exp(-t) * sin(2π*x)

This is periodic: u(0,t) = u(1,t) = 0 for all t

Derivatives:
- u_t = -exp(-t) * sin(2π*x)
- u_x = 2π * exp(-t) * cos(2π*x)
- u_xx = -4π² * exp(-t) * sin(2π*x)

PDE: u_t = κ*u_xx + f
-exp(-t)*sin(2π*x) = κ*(-4π²)*exp(-t)*sin(2π*x) + f(x,t)

Forcing term:
f(x,t) = -exp(-t)*sin(2π*x) + 4π²*κ*exp(-t)*sin(2π*x)
       = exp(-t)*sin(2π*x) * (4π²*κ - 1)
"""

κ = 0.1
t_final = 0.1

# Manufactured solution
u_manufactured = ManufacturedSolution(
    (x, t) -> exp(-t) * sin(2π * x),
    (0.0, 1.0),
    name = "exp(-t)*sin(2πx)",
    description = "Time-decaying periodic solution"
)

# Initial condition
u_initial(x) = u_manufactured(x, 0.0)

# Forcing function (derived analytically)
forcing_function(x, t) = exp(-t) * sin(2π * x) * (4 * π^2 * κ - 1)

# ============================================================================
# Create solver wrapper
# ============================================================================

function heat_solver_wrapper(grid::Grid1D, manufactured_solution::ManufacturedSolution; 
                             κ=0.1, t_final=0.1, cfl=0.25)
    return heat_solver_1d_periodic(
        grid, 
        forcing_function, 
        u_initial,
        κ=κ, 
        t_final=t_final, 
        cfl=cfl
    )
end

# Wrapper for evaluation at final time
struct TimeDependentSolution
    mms::ManufacturedSolution
    t_final::Float64
end

function (tds::TimeDependentSolution)(x)
    return tds.mms(x, tds.t_final)
end

# Create time-dependent manufactured solution
u_mms_at_tfinal = ManufacturedSolution(
    TimeDependentSolution(u_manufactured, t_final),
    (0.0, 1.0),
    name = "exp(-$t_final)*sin(2πx)",
    description = "Solution at t = $t_final"
)

# ============================================================================
# Run Verification
# ============================================================================

println("="^70)
println("1D Heat Equation Solver Verification using MMS")
println("="^70)
println()
println("Problem: u_t = κ*u_xx + f on [0,1] with periodic BCs")
println("Thermal diffusivity: κ = $κ")
println("Final time: t = $t_final")
println("Manufactured Solution: u(x,t) = exp(-t)*sin(2πx)")
println()

solver = UserSolver(
    heat_solver_wrapper,
    name = "1D Forward Euler Heat Solver",
    description = "1st order in time, 2nd order in space"
)

println("Running convergence study...")
println()

results = verify_solver(
    solver,
    nothing,
    u_mms_at_tfinal,
    resolutions = [0.05, 0.025, 0.0125, 0.00625],
    verbose = true,
    κ = κ,
    t_final = t_final,
    cfl = 0.25
)

println()
println("="^70)
println("Verification Summary")
println("="^70)

println("\nConvergence Rates:")
println("-"^50)
println("Expected: ~2.0 in space (2nd order spatial discretization)")
println("Note: Forward Euler is 1st order in time, which may affect overall rate")
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
if l2_rate > 1.5  # Spatial convergence dominates for small time steps
    println("✓ PASSED: Spatial convergence is being observed")
    println("  L2 rate = $(round(l2_rate, digits=3))")
else
    println("⚠ Note: Convergence rate = $(round(l2_rate, digits=3))")
    println("  This may be influenced by temporal discretization error")
end

println()
println("="^70)
