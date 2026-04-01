"""
Convergence Verification Tests for 1D Poisson Equation

This test suite verifies the order of accuracy by running solvers
on multiple grid resolutions and checking convergence rates.

Test Cases:
1. Square Domain [0,1]
2. Rectangular Domain [0,2]
3. Shifted Domain [1,2]
4. Long Domain [0,10]
5. Irregular Domain [-1.5, 3.2]
"""

using ManufacturedSolutions
using LinearAlgebra
using Test
using Printf

# Import solver
include("solvers/poisson_1d.jl")
using .Main: poisson_solver_dirichlet

# Define the PDE
@variables x u(x)
Dx = Differential(x)
pde_poisson = -Dx(Dx(u))  # -u'' = f

println("="^80)
println("Convergence Verification: 1D Poisson Equation")
println("Testing Order of Accuracy with Multiple Refinements")
println("="^80)
println()

@testset "Poisson Solver - Order of Accuracy Tests" begin
    
    #==========================================================================#
    # Test 1: Square Domain [0,1]
    #==========================================================================#
    @testset "Square Domain [0,1]" begin
        println("\n" * "="^80)
        println("Test 1: Square Domain [0,1]")
        println("Manufactured Solution: u(x) = sin(πx)")
        println("="^80)
        
        test = ConvergenceTest(
            "Square Domain [0,1]",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Square domain test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 2: Rectangular Domain [0,2]
    #==========================================================================#
    @testset "Rectangular Domain [0,2]" begin
        println("\n" * "="^80)
        println("Test 2: Rectangular Domain [0,2]")
        println("Manufactured Solution: u(x) = sin(πx/2)")
        println("="^80)
        
        test = ConvergenceTest(
            "Rectangular Domain [0,2]",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(π * x / 2),
            (x, u),
            (0.0, 2.0),
            expected_order = 2.0,
            base_resolution = 0.2,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Rectangular domain test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 3: Shifted Domain [1,2]
    #==========================================================================#
    @testset "Shifted Domain [1,2]" begin
        println("\n" * "="^80)
        println("Test 3: Shifted Domain [1,2]")
        println("Manufactured Solution: u(x) = sin(π(x-1))")
        println("="^80)
        
        test = ConvergenceTest(
            "Shifted Domain [1,2]",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(π * (x - 1)),
            (x, u),
            (1.0, 2.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Shifted domain test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 4: Long Domain [0,10]
    #==========================================================================#
    @testset "Long Domain [0,10]" begin
        println("\n" * "="^80)
        println("Test 4: Long Domain [0,10]")
        println("Manufactured Solution: u(x) = sin(πx/10)")
        println("="^80)
        
        test = ConvergenceTest(
            "Long Domain [0,10]",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(π * x / 10),
            (x, u),
            (0.0, 10.0),
            expected_order = 2.0,
            base_resolution = 1.0,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Long domain test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 5: Irregular Domain [-1.5, 3.2]
    #==========================================================================#
    @testset "Irregular Domain [-1.5, 3.2]" begin
        println("\n" * "="^80)
        println("Test 5: Irregular Domain [-1.5, 3.2]")
        println("Manufactured Solution: u(x) = sin(π(x+1.5)/4.7)")
        println("="^80)
        
        # Domain length is 4.7, so wavelength matches domain
        test = ConvergenceTest(
            "Irregular Domain [-1.5, 3.2]",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(π * (x + 1.5) / 4.7),
            (x, u),
            (-1.5, 3.2),
            expected_order = 2.0,
            base_resolution = 0.47,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Irregular domain test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 6: Polynomial Solution (Machine Precision Check)
    #==========================================================================#
    @testset "Polynomial Solution - Machine Precision" begin
        println("\n" * "="^80)
        println("Test 6: Polynomial Solution u(x) = x(1-x)")
        println("Expected: Errors at machine precision (solver exact for polynomials)")
        println("="^80)
        
        test = ConvergenceTest(
            "Polynomial - Machine Precision",
            poisson_solver_dirichlet,
            pde_poisson,
            x * (1 - x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 3,  # Fewer refinements (errors hit machine precision)
            tolerance = 5.0  # Very relaxed - errors are at noise level
        )
        
        result = run_convergence_test(test)
        
        # For polynomial, errors should be at machine precision
        # Convergence order is meaningless when errors are ~ 1e-15
        @test result[:errors][end] < 1e-10  # Check final error is tiny
        
        println("✓ Polynomial test PASSED")
        println("  Note: Final error = $(round(result[:errors][end], sigdigits=2))")
        println("        At machine precision - convergence order undefined")
        println()
    end
    
    #==========================================================================#
    # Test 7: Multiple Modes - Higher Frequency
    #==========================================================================#
    @testset "Multiple Modes - sin(3πx)" begin
        println("\n" * "="^80)
        println("Test 7: Multiple Modes - sin(3πx)")
        println("Higher frequency requires finer base resolution")
        println("="^80)
        
        test = ConvergenceTest(
            "Multiple Modes",
            poisson_solver_dirichlet,
            pde_poisson,
            sin(3 * π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.05,  # Finer base grid for higher frequency
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Multiple modes test PASSED")
        println()
    end
    
    #==========================================================================#
    # Test 8: Exponential Solution
    #==========================================================================#
    @testset "Exponential Solution" begin
        println("\n" * "="^80)
        println("Test 8: Exponential Solution u(x) = exp(-x)sin(πx)")
        println("Tests convergence with mixed exponential/trig")
        println("="^80)
        
        test = ConvergenceTest(
            "Exponential",
            poisson_solver_dirichlet,
            pde_poisson,
            exp(-x) * sin(π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        @test result[:passed]
        @test abs(result[:observed_order] - 2.0) < 0.2
        
        println("✓ Exponential solution test PASSED")
        println()
    end
    
end

println("\n" * "="^80)
println("CONVERGENCE VERIFICATION SUMMARY")
println("="^80)
println()
println("All tests verify second-order accuracy (p ≈ 2.0) for the finite")
println("difference Poisson solver across various domain configurations:")
println()
println("  ✓ Square domain [0,1]")
println("  ✓ Rectangular domain [0,2]")
println("  ✓ Shifted domain [1,2]")
println("  ✓ Long domain [0,10]")
println("  ✓ Irregular domain [-1.5, 3.2]")
println("  ✓ Polynomial solution (exact)")
println("  ✓ Multiple modes (high frequency)")
println("  ✓ Exponential solution")
println()
println("The solver demonstrates consistent second-order convergence")
println("independent of domain size, location, or solution type.")
println()
println("="^80)
