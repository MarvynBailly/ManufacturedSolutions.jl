"""
Convergence analysis and order of accuracy estimation
"""

using Statistics

"""
    estimate_convergence_rate(errors::Vector{T}, resolutions::Vector{T}) where T

Estimate the convergence rate (order of accuracy) from error data.

Uses log-log linear regression: log(error) ≈ p*log(h) + c
where p is the convergence rate.

# Arguments
- `errors`: Vector of errors at different resolutions (from fine to coarse or vice versa)
- `resolutions`: Vector of grid spacings (same order as errors)

# Returns
- `p`: Estimated convergence rate
- `R²`: Coefficient of determination (goodness of fit)
"""
function estimate_convergence_rate(errors::Vector{T}, resolutions::Vector{T}) where T
    @assert length(errors) == length(resolutions) >= 2 "Need at least 2 data points"
    @assert all(errors .> 0) "Errors must be positive"
    @assert all(resolutions .> 0) "Resolutions must be positive"
    
    # Log-log transformation
    log_h = log.(resolutions)
    log_err = log.(errors)
    
    # Linear regression: log_err = p * log_h + c
    n = length(log_h)
    mean_log_h = mean(log_h)
    mean_log_err = mean(log_err)
    
    # Slope (convergence rate)
    numerator = sum((log_h .- mean_log_h) .* (log_err .- mean_log_err))
    denominator = sum((log_h .- mean_log_h).^2)
    p = numerator / denominator
    
    # Compute R² (coefficient of determination)
    predictions = p .* log_h .+ (mean_log_err - p * mean_log_h)
    ss_res = sum((log_err .- predictions).^2)
    ss_tot = sum((log_err .- mean_log_err).^2)
    R² = 1 - ss_res / ss_tot
    
    return p, R²
end

"""
    estimate_convergence_rate_pairwise(errors::Vector{T}, resolutions::Vector{T}) where T

Estimate convergence rate using consecutive pairs of data points.

For each consecutive pair, computes: p = log(e₁/e₂) / log(h₁/h₂)

# Returns
- Vector of pairwise convergence rates
"""
function estimate_convergence_rate_pairwise(errors::Vector{T}, resolutions::Vector{T}) where T
    @assert length(errors) == length(resolutions) >= 2 "Need at least 2 data points"
    
    n = length(errors)
    rates = zeros(n-1)
    
    for i in 1:(n-1)
        rates[i] = log(errors[i] / errors[i+1]) / log(resolutions[i] / resolutions[i+1])
    end
    
    return rates
end

"""
    analyze_convergence(convergence_data::Vector{ConvergenceData{T}}) where T

Analyze convergence data and compute convergence rates for all error norms.

# Returns
- Dictionary mapping error norm names to convergence rates
"""
function analyze_convergence(convergence_data::Vector{ConvergenceData{T}}) where T
    @assert length(convergence_data) >= 2 "Need at least 2 resolutions for convergence analysis"
    
    # Extract resolutions
    resolutions = [data.resolution for data in convergence_data]
    
    # Get all error norm types from first data point
    error_types = keys(convergence_data[1].errors)
    
    convergence_rates = Dict{Symbol,T}()
    
    for error_type in error_types
        # Extract errors for this norm
        errors = [data.errors[error_type] for data in convergence_data]
        
        # Estimate convergence rate
        rate, R² = estimate_convergence_rate(errors, resolutions)
        
        # Store the rate
        convergence_rates[error_type] = rate
        
        # Also store R² for quality assessment
        convergence_rates[Symbol(string(error_type) * "_R2")] = R²
    end
    
    return convergence_rates
end

export estimate_convergence_rate, estimate_convergence_rate_pairwise, analyze_convergence
