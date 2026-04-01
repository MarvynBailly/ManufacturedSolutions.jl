"""
Main test runner for ManufacturedSolutions.jl

This file runs all test suites when you execute `Pkg.test()`.
"""

using Test
using ManufacturedSolutions

@testset "ManufacturedSolutions.jl Tests" begin
    
    @testset "Simple MMS Tests - Poisson 1D" begin
        include("verification_1d_poisson.jl")
    end
    
    @testset "Simple MMS Tests - Reaction-Diffusion" begin
        include("verification_reaction_diffusion.jl")
    end
    
    @testset "Convergence Tests - Poisson 1D" begin
        include("convergence_poisson_1d.jl")
    end
    
    @testset "Convergence Tests - Reaction-Diffusion" begin
        include("convergence_reaction_diffusion.jl")
    end
    
end

println("\n" * "="^70)
println("All tests completed successfully! ✓")
println("="^70)
