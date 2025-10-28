"""
1D Reaction-Diffusion ODE Solvers for Testing

Solves: -u'' + k*u = f(x) on [a,b]

This is a second-order ODE combining diffusion (-u'') and reaction (+k*u).
"""

using LinearAlgebra

"""
    reaction_diffusion_solver_correct(grid::Grid1D, forcing::ForcingTerm; k=1.0)

CORRECT implementation of 1D reaction-diffusion solver with Dirichlet BCs.

Solves: -u'' + k*u = f on [a,b] with u(a) = u(b) = 0

Parameters:
- k: Reaction coefficient (default: 1.0)

Uses centered finite differences with second-order accuracy.
"""
function reaction_diffusion_solver_correct(grid::Grid1D, forcing::ForcingTerm; k=1.0)
    nx = grid.nx
    h = grid.dx
    
    # Allocate
    u = zeros(nx)
    
    # Build tridiagonal system
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left BC: u(a) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        elseif i == nx
            # Right BC: u(b) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        else
            # Interior: -u'' + k*u = f
            # Discretization: -(u[i-1] - 2u[i] + u[i+1])/h^2 + k*u[i] = f[i]
            # Rearranged: -u[i-1] + (2 + k*h^2)*u[i] - u[i+1] = h^2*f[i]
            A[i, i-1] = -1.0
            A[i, i] = 2.0 + k * h^2  # ✓ CORRECT: includes reaction term
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve linear system
    u = A \ rhs
    return u
end

"""
    reaction_diffusion_solver_broken(grid::Grid1D, forcing::ForcingTerm; k=1.0)

INTENTIONALLY BROKEN implementation to demonstrate verification failure.

This solver FORGETS the reaction term k*u in the discretization!
It will fail MMS tests and show higher errors.
"""
function reaction_diffusion_solver_broken(grid::Grid1D, forcing::ForcingTerm; k=1.0)
    nx = grid.nx
    h = grid.dx
    
    # Allocate
    u = zeros(nx)
    
    # Build tridiagonal system
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left BC: u(a) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        elseif i == nx
            # Right BC: u(b) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        else
            # Interior: -u'' + k*u = f
            # ❌ BUG: Missing the k*u term!
            A[i, i-1] = -1.0
            A[i, i] = 2.0  # ❌ Should be 2.0 + k*h^2
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve linear system
    u = A \ rhs
    return u
end

"""
    reaction_diffusion_solver_very_broken(grid::Grid1D, forcing::ForcingTerm; k=1.0)

VERY BROKEN implementation with wrong sign in discretization.

This solver has the WRONG SIGN for the reaction term.
It solves -u'' - k*u = f instead of -u'' + k*u = f
"""
function reaction_diffusion_solver_very_broken(grid::Grid1D, forcing::ForcingTerm; k=1.0)
    nx = grid.nx
    h = grid.dx
    
    # Allocate
    u = zeros(nx)
    
    # Build tridiagonal system
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left BC: u(a) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        elseif i == nx
            # Right BC: u(b) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        else
            # Interior: -u'' + k*u = f
            # ❌ BUG: Wrong sign! Should be + k*h^2, not - k*h^2
            A[i, i-1] = -1.0
            A[i, i] = 2.0 - k * h^2  # ❌ Wrong sign!
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve linear system
    u = A \ rhs
    return u
end

export reaction_diffusion_solver_correct, reaction_diffusion_solver_broken, 
       reaction_diffusion_solver_very_broken
