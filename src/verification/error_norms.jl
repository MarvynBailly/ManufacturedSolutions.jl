"""
Error computation and norms for verification
"""

using LinearAlgebra

"""
    compute_L1_error(numerical, analytical, grid::Grid1D)

Compute the L1 (integral) error norm.

L1 = ∫|u_numerical - u_analytical| dx ≈ Δx * Σ|u_num[i] - u_analytical[i]|
"""
function compute_L1_error(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    @assert length(numerical) == length(analytical) == grid.nx "Array sizes must match grid size"
    
    error_sum = sum(abs.(numerical .- analytical))
    return grid.dx * error_sum
end

"""
    compute_L2_error(numerical, analytical, grid::Grid1D)

Compute the L2 (root mean square) error norm.

L2 = √(∫|u_numerical - u_analytical|² dx) ≈ √(Δx * Σ|u_num[i] - u_analytical[i]|²)
"""
function compute_L2_error(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    @assert length(numerical) == length(analytical) == grid.nx "Array sizes must match grid size"
    
    error_squared_sum = sum(abs2.(numerical .- analytical))
    return sqrt(grid.dx * error_squared_sum)
end

"""
    compute_Linf_error(numerical, analytical)

Compute the L∞ (maximum) error norm.

L∞ = max|u_numerical - u_analytical|
"""
function compute_Linf_error(numerical::AbstractVector, analytical::AbstractVector)
    @assert length(numerical) == length(analytical) "Array sizes must match"
    
    return maximum(abs.(numerical .- analytical))
end

"""
    compute_all_errors(numerical, analytical, grid::Grid1D)

Compute all standard error norms and return as a dictionary.
"""
function compute_all_errors(numerical::AbstractVector, analytical::AbstractVector, grid::Grid1D)
    return Dict{Symbol,Float64}(
        :L1 => compute_L1_error(numerical, analytical, grid),
        :L2 => compute_L2_error(numerical, analytical, grid),
        :Linf => compute_Linf_error(numerical, analytical)
    )
end

export compute_L1_error, compute_L2_error, compute_Linf_error, compute_all_errors
