"""
Reference 1D Poisson Solvers for MMS Testing

This file contains **educational reference implementations** of 1D Poisson
solvers with different boundary conditions. These are used as test cases to
demonstrate the Method of Manufactured Solutions framework.

## The Poisson Equation

We solve the 1D Poisson equation:
    -u''(x) = f(x)  for x ∈ [a, b]

with various boundary conditions:
- **Dirichlet**: Specify solution values at boundaries (u(a), u(b) given)
- **Periodic**: Solution wraps around (u(a) = u(b), u'(a) = u'(b))
- **Neumann**: Specify flux (derivative) at boundaries (u'(a), u'(b) given)

## Discretization

All solvers use **second-order centered finite differences**:
    u''(x) ≈ [u(x-h) - 2u(x) + u(x+h)] / h²

This leads to a **tridiagonal linear system** Au = f that we solve directly.

## Why These Solvers?

These are intentionally simple, well-tested implementations that:
1. Demonstrate correct numerical methods
2. Achieve expected second-order accuracy
3. Serve as reference solutions for MMS verification
4. Show how to handle different boundary conditions

If your solver fails MMS tests but these pass, the issue is in your implementation,
not in the verification framework!
"""

using LinearAlgebra

"""
    poisson_solver_dirichlet(grid::Grid1D, forcing::ForcingTerm)

Solve the 1D Poisson equation with **Dirichlet (fixed value) boundary conditions**.

## Problem
    -u''(x) = f(x)  for x ∈ [a, b]
    u(a) = 0,  u(b) = 0

## Method
Uses **centered finite differences** with second-order accuracy:
    -u''(x) ≈ -[u(x-h) - 2u(x) + u(x+h)]/h² = f(x)

This discretization leads to a tridiagonal system:
```
[  1    0    0  ...  0  ]   [ u₁ ]   [ 0     ]
[ -1    2   -1  ...  0  ]   [ u₂ ]   [ h²f₂  ]
[  0   -1    2  ...  0  ]   [ u₃ ] = [ h²f₃  ]
[ ...               ... ]   [ .. ]   [ ...   ]
[  0    0    0  ...  1  ]   [ uₙ ]   [ 0     ]
```

The first and last rows enforce boundary conditions u(a) = 0, u(b) = 0.

## Arguments
- `grid`: Spatial discretization (Grid1D)
- `forcing`: Right-hand side function f(x)

## Returns
- Solution vector u at grid points

## Convergence
This method is **second-order accurate**: error ~ O(h²)

## Example
```julia
grid = Grid1D(0.0, 1.0, 101)
forcing = ForcingTerm(nothing, x -> π^2 * sin(π*x))
u = poisson_solver_dirichlet(grid, forcing)
# Exact solution: u(x) = sin(πx)
```
"""
function poisson_solver_dirichlet(grid::Grid1D, forcing::ForcingTerm)
    nx = grid.nx
    h = grid.dx
    
    # Allocate solution vector
    u = zeros(nx)
    
    # Build tridiagonal system Au = rhs
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left boundary: u(a) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        elseif i == nx
            # Right boundary: u(b) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        else
            # Interior points: -u'' = f
            # Stencil: -(u[i-1] - 2u[i] + u[i+1])/h² = f[i]
            # Rearranging: -u[i-1] + 2u[i] - u[i+1] = h²f[i]
            A[i, i-1] = -1.0
            A[i, i] = 2.0
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing(grid.x[i])
        end
    end
    
    # Solve the linear system (direct solve for tridiagonal)
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
