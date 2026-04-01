"""
Automatic computation of forcing terms for the Method of Manufactured Solutions

This module uses symbolic mathematics (via Symbolics.jl) to automatically compute
the forcing term needed to produce a given manufactured solution. This eliminates
the need for manual differentiation and reduces errors.

# The Core Idea

For a PDE:  L[u] = f(x)

Where L is a differential operator (e.g., -d²/dx² for Poisson), we want to find
the forcing f(x) that produces a known solution u(x).

Solution: Substitute u(x) into L[u] and simplify to get f(x) = L[u(x)]

# Example

PDE: -u'' = f on [0,1]
Choose: u(x) = sin(πx)
Compute: u'(x) = π cos(πx)
         u''(x) = -π² sin(πx)
Therefore: f(x) = -u''(x) = π² sin(πx)

This module does all that computation automatically!
"""

using Symbolics

"""
    compute_forcing_term(pde_lhs, manufactured_solution_sym, variables)

Compute the forcing term by substituting the manufactured solution into the PDE.

For a PDE of the form: L[u] = f
where L is a differential operator, we substitute u = u_manufactured to get:
f = L[u_manufactured]

# Arguments
- `pde_lhs`: Left-hand side of PDE (symbolic expression with differential operators)
- `manufactured_solution_sym`: Symbolic manufactured solution
- `variables`: Tuple of symbolic variables (e.g., (x,) or (x, t))

# Returns
- `ForcingTerm`: Forcing term with both symbolic and callable forms

# Example
```julia
@variables x u(x)
Dx = Differential(x)

# PDE: -u'' = f
pde_lhs = -Dx(Dx(u))

# Manufactured solution: u(x) = sin(πx)
u_mms = sin(π * x)

# Compute forcing: f(x) = -(-π²sin(πx)) = π²sin(πx)
forcing = compute_forcing_term(pde_lhs, u_mms, (x,))
```
"""
function compute_forcing_term(pde_lhs, manufactured_solution_sym, variables)
    # The PDE left-hand side contains differential operators acting on u.
    # For example, for Poisson: pde_lhs = -Dx(Dx(u))
    # We need to substitute our manufactured solution for u.
    
    # Substitute: Replace u(x) with the manufactured solution (e.g., sin(πx))
    pde_with_mms = substitute(pde_lhs, Dict(variables[end] => manufactured_solution_sym))
    
    # Expand derivatives: Symbolics.jl needs to actually compute the derivatives
    # This turns -Dx(Dx(sin(πx))) into -(-π² sin(πx))
    forcing_symbolic = expand_derivatives(pde_with_mms)
    
    # Simplify: Clean up the expression (combine like terms, etc.)
    forcing_symbolic = simplify(forcing_symbolic)
    
    # Build callable function: Convert symbolic expression to actual Julia function
    # We need only the independent variables (x, t, etc.), not the dependent variable (u)
    independent_vars = variables[1:end-1]  # e.g., (x,) or (x, t)
    
    # This creates a function we can call like: forcing_func(0.5) to get f(0.5)
    forcing_func = build_function(forcing_symbolic, independent_vars..., expression=Val{false})
    
    # Return both the function and symbolic form (useful for debugging/display)
    return ForcingTerm(forcing_func, forcing_symbolic)
end

"""
    compute_forcing_term_1d(pde_lhs, manufactured_solution_sym, x_var, u_var)

Convenience function for 1D spatial problems.

# Arguments
- `pde_lhs`: Left-hand side of PDE
- `manufactured_solution_sym`: Symbolic manufactured solution expression
- `x_var`: Spatial variable symbol
- `u_var`: Dependent variable symbol

# Example
```julia
@variables x u(x)
Dx = Differential(x)

pde_lhs = -Dx(Dx(u))
u_mms = sin(π * x)

forcing = compute_forcing_term_1d(pde_lhs, u_mms, x, u)
```
"""
function compute_forcing_term_1d(pde_lhs, manufactured_solution_sym, x_var, u_var)
    return compute_forcing_term(pde_lhs, manufactured_solution_sym, (x_var, u_var))
end

"""
    compute_forcing_term_time_dependent(pde_lhs, manufactured_solution_sym, x_var, t_var, u_var)

Convenience function for time-dependent 1D problems.

# Arguments
- `pde_lhs`: Left-hand side of PDE
- `manufactured_solution_sym`: Symbolic manufactured solution expression
- `x_var`: Spatial variable symbol
- `t_var`: Time variable symbol  
- `u_var`: Dependent variable symbol

# Example
```julia
@variables x t u(x,t)
Dx = Differential(x)
Dt = Differential(t)

# Heat equation: u_t - κ*u_xx = f
κ = 0.1
pde_lhs = Dt(u) - κ*Dx(Dx(u))
u_mms = exp(-t) * sin(2π * x)

forcing = compute_forcing_term_time_dependent(pde_lhs, u_mms, x, t, u)
```
"""
function compute_forcing_term_time_dependent(pde_lhs, manufactured_solution_sym, x_var, t_var, u_var)
    # Substitute manufactured solution into PDE
    pde_with_mms = substitute(pde_lhs, Dict(u_var => manufactured_solution_sym))
    
    # Expand all derivatives
    forcing_symbolic = expand_derivatives(pde_with_mms)
    
    # Simplify if possible
    forcing_symbolic = simplify(forcing_symbolic)
    
    # Build a callable function from the symbolic expression
    forcing_func = build_function(forcing_symbolic, [x_var, t_var], expression=Val{false})
    
    return ForcingTerm(forcing_func, forcing_symbolic)
end

"""
    @pde_system(expr)

Macro for defining PDE systems with a clean syntax.

# Example
```julia
@pde_system begin
    @variables x u(x)
    @operator Dx = Differential(x)
    @pde -Dx(Dx(u)) = 0  # Laplace equation
end
```
"""
macro pde_system(expr)
    # This is a placeholder for future enhancement
    # For now, users can use the functions directly
    error("@pde_system macro not yet implemented. Use compute_forcing_term functions directly.")
end

export compute_forcing_term, compute_forcing_term_1d, compute_forcing_term_time_dependent
