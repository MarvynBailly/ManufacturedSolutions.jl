"""
Error Norms for Numerical Verification

This module computes discrete error norms to measure how closely a numerical
solution matches the exact (manufactured) solution. We use three standard norms:

- **L1 norm**: Average absolute error (integral norm)
- **L2 norm**: Root mean square error (energy norm)  
- **L∞ norm**: Maximum pointwise error (supremum norm)

These norms measure different aspects of solution accuracy and are commonly
used in convergence studies to verify that errors decrease at the expected rate.

# Example
```julia
grid = Grid1D(0.0, 1.0, 101)
u_numerical = solve_my_pde(grid)
u_exact = sin.(π .* grid.x)

errors = compute_all_errors(u_numerical, u_exact, grid)
# errors[:L1]   = 0.001234  (average error)
# errors[:L2]   = 0.001567  (RMS error)
# errors[:Linf] = 0.002341  (max error)
```
"""

using LinearAlgebra

"""
    compute_L1_error(numerical, analytical, grid::Grid1D)

Compute the L1 (integral) norm of the error.

The L1 norm measures the **average absolute error** across the domain. It's
computed using the discrete trapezoidal approximation to the integral:

    L1 = ∫|u_numerical - u_exact| dx ≈ Δx × Σ|u_num[i] - u_exact[i]|

This norm gives equal weight to errors at all points and is robust to isolated
spikes (unlike L∞). For smooth solutions, L1 error typically converges one
order slower than L2 error.

# Arguments
- `numerical`: Numerical solution values at grid points
- `analytical`: Exact solution values at grid points  
- `grid`: Grid structure (provides spacing Δx)

# Returns
- L1 error norm (scalar, same units as solution)

# Example
```julia
grid = Grid1D(0.0, 1.0, 101)
u_num = [0.0, 0.1, 0.2, ...]
u_exact = sin.(π .* grid.x)
error_L1 = compute_L1_error(u_num, u_exact, grid)
```
"""
function compute_L1_error(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    @assert length(numerical) == length(analytical) == grid.nx "Array sizes must match grid size"
    
    # Sum absolute pointwise errors
    error_sum = sum(abs.(numerical .- analytical))
    
    # Scale by grid spacing to approximate integral
    return grid.dx * error_sum
end

"""
    compute_L2_error(numerical, analytical, grid::Grid1D)

Compute the L2 (energy) norm of the error.

The L2 norm measures the **root mean square error**, which relates to the
"energy" of the error in many physical systems:

    L2 = √(∫|u_numerical - u_exact|² dx) ≈ √(Δx × Σ|u_num[i] - u_exact[i]|²)

This is the most commonly used norm in convergence studies because it:
- Penalizes larger errors more heavily (due to squaring)
- Connects to energy estimates in stability analysis
- Often shows cleaner convergence behavior than L1

For a second-order accurate method, L2 error typically scales like h².

# Arguments
- `numerical`: Numerical solution values at grid points
- `analytical`: Exact solution values at grid points
- `grid`: Grid structure (provides spacing Δx)

# Returns
- L2 error norm (scalar, same units as solution)

# Example
```julia
grid = Grid1D(0.0, 1.0, 101)
error_L2 = compute_L2_error(u_numerical, u_exact, grid)
println("L2 error: ", error_L2)
```
"""
function compute_L2_error(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    @assert length(numerical) == length(analytical) == grid.nx "Array sizes must match grid size"
    
    # Sum squared pointwise errors
    error_squared_sum = sum(abs2.(numerical .- analytical))
    
    # Scale by grid spacing and take square root
    return sqrt(grid.dx * error_squared_sum)
end

"""
    compute_Linf_error(numerical, analytical)

Compute the L∞ (maximum) norm of the error.

The L∞ norm measures the **worst-case pointwise error** anywhere in the domain:

    L∞ = max|u_numerical - u_exact|

This norm answers the question: "What's the largest error at any single point?"
It's useful for:
- Safety-critical applications (worst-case guarantees)
- Debugging (finding where errors are largest)
- Detecting local issues (boundary errors, discontinuities)

Note: L∞ doesn't require grid spacing since it's just a maximum over points.
It can be more sensitive to isolated bad points than L1/L2 norms.

# Arguments
- `numerical`: Numerical solution values at grid points
- `analytical`: Exact solution values at grid points

# Returns
- L∞ error norm (scalar, same units as solution)

# Example
```julia
error_Linf = compute_Linf_error(u_numerical, u_exact)
println("Maximum error: ", error_Linf)
```
"""
function compute_Linf_error(numerical::AbstractVector, analytical::AbstractVector)
    @assert length(numerical) == length(analytical) "Array sizes must match"
    
    # Find maximum absolute pointwise error
    return maximum(abs.(numerical .- analytical))
end

"""
    compute_all_errors(numerical, analytical, grid::Grid1D)

Compute all three standard error norms (L1, L2, L∞) and return as a dictionary.

This convenience function computes all error norms at once, which is useful for:
- Convergence studies (check multiple norms simultaneously)
- Comprehensive error reporting
- Comparing different aspects of solution quality

The L2 norm is typically most important for convergence order verification,
while L∞ helps identify localized issues.

# Arguments
- `numerical`: Numerical solution values at grid points
- `analytical`: Exact solution values at grid points
- `grid`: Grid structure

# Returns
- Dictionary with keys `:L1`, `:L2`, `:Linf` mapping to error values

# Example
```julia
grid = Grid1D(0.0, 1.0, 101)
errors = compute_all_errors(u_numerical, u_exact, grid)

println("L1 error:   ", errors[:L1])
println("L2 error:   ", errors[:L2])
println("Linf error: ", errors[:Linf])
```
"""
function compute_all_errors(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    return Dict{Symbol,Float64}(
        :L1 => compute_L1_error(numerical, analytical, grid),
        :L2 => compute_L2_error(numerical, analytical, grid),
        :Linf => compute_Linf_error(numerical, analytical)
    )
end

export compute_L1_error, compute_L2_error, compute_Linf_error, compute_all_errors

