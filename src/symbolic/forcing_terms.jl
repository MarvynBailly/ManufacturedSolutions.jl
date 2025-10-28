"""
Symbolic processing for PDE forcing term computation
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
    # Substitute manufactured solution into PDE
    # The PDE LHS contains operations on u (the unknown)
    # We need to replace u with the manufactured solution
    pde_with_mms = substitute(pde_lhs, Dict(variables[end] => manufactured_solution_sym))
    
    # Expand all derivatives
    forcing_symbolic = expand_derivatives(pde_with_mms)
    
    # Simplify if possible
    forcing_symbolic = simplify(forcing_symbolic)
    
    # Build a callable function from the symbolic expression
    # Remove the independent variables from the tuple (keep only spatial/temporal vars)
    independent_vars = variables[1:end-1]  # All but the last (which is the dependent variable)
    
    forcing_func = build_function(forcing_symbolic, independent_vars..., expression=Val{false})
    
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
