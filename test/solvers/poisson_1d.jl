"""
1D Poisson Solvers for Testing

This file contains reference implementations of 1D Poisson solvers
that can be verified using MMS unit tests.

Solves: -u'' = f(x) on [a,b]
"""

using LinearAlgebra

"""
    poisson_solver_dirichlet(grid::Grid1D, forcing::ForcingTerm)

Simple 1D Poisson solver with Dirichlet boundary conditions.

Solves: -u'' = f on [a,b] with u(a) = u(b) = 0

Uses centered finite differences with second-order accuracy.
"""
function poisson_solver_dirichlet(grid::Grid1D, forcing::ForcingTerm)
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
            # Interior: -u'' = f
            # Discretization: -(u[i-1] - 2u[i] + u[i+1])/h^2 = f[i]
            A[i, i-1] = -1.0
            A[i, i] = 2.0
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve linear system
    u = A \ rhs
    return u
end

"""
    poisson_solver_periodic(grid::Grid1D, forcing::ForcingTerm)

1D Poisson solver with periodic boundary conditions.

Solves: -u'' = f with u(a) = u(b), u'(a) = u'(b)

Note: For this to have a solution, the forcing must satisfy
the compatibility condition: ∫f dx = 0
"""
function poisson_solver_periodic(grid::Grid1D, forcing::ForcingTerm)
    nx = grid.nx
    h = grid.dx
    
    # For periodic BCs, we have nx-1 unknowns (last point = first point)
    n = nx - 1
    u = zeros(n)
    
    # Build system with periodic connections
    A = zeros(n, n)
    rhs = zeros(n)
    
    for i in 1:n
        # Periodic wrapping
        im1 = i == 1 ? n : i - 1
        ip1 = i == n ? 1 : i + 1
        
        A[i, im1] = -1.0
        A[i, i] = 2.0
        A[i, ip1] = -1.0
        rhs[i] = h^2 * forcing(grid.x[i])
    end
    
    # Solve
    u = A \ rhs
    
    # Add periodic point
    return vcat(u, u[1])
end

"""
    poisson_solver_neumann(grid::Grid1D, forcing::ForcingTerm; left_flux=0.0, right_flux=0.0)

1D Poisson solver with Neumann (flux) boundary conditions.

Solves: -u'' = f with -u'(a) = left_flux, -u'(b) = right_flux

Note: For this to have a solution, the forcing and fluxes must satisfy
the compatibility condition: ∫f dx + right_flux - left_flux = 0
"""
function poisson_solver_neumann(grid::Grid1D, forcing::ForcingTerm; 
                               left_flux=0.0, right_flux=0.0)
    nx = grid.nx
    h = grid.dx
    
    # Build system
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left Neumann BC: -u'(a) = left_flux
            # Using one-sided difference: -(u[2] - u[1])/h = left_flux
            A[i, i] = -1.0/h
            A[i, i+1] = 1.0/h
            rhs[i] = left_flux
        elseif i == nx
            # Right Neumann BC: -u'(b) = right_flux
            # Using one-sided difference: -(u[nx] - u[nx-1])/h = right_flux
            A[i, i-1] = -1.0/h
            A[i, i] = 1.0/h
            rhs[i] = right_flux
        else
            # Interior: -u'' = f
            A[i, i-1] = -1.0
            A[i, i] = 2.0
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve (may be singular if compatibility not satisfied)
    u = A \ rhs
    return u
end

export poisson_solver_dirichlet, poisson_solver_periodic, poisson_solver_neumann
