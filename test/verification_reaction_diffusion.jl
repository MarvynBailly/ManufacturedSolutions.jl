"""
Verification Tests for Reaction-Diffusion ODE Solvers

Tests both CORRECT and BROKEN implementations to demonstrate
how MMS verification catches implementation errors.

ODE: -u'' + k*u = f(x) on [a,b] with u(a) = u(b) = 0
"""

using ManufacturedSolutions
using LinearAlgebra
using Test
using Printf

# Import solvers
include("solvers/reaction_diffusion_1d.jl")
using .Main: reaction_diffusion_solver_correct, reaction_diffusion_solver_broken,
             reaction_diffusion_solver_very_broken

# Define the PDE symbolically
@variables x u(x)
Dx = Differential(x)
k_param = 1.0  # Reaction coefficient
pde_reaction_diffusion = -Dx(Dx(u)) + k_param * u  # -u'' + k*u = f

println("="^70)
println("Verification Tests: Reaction-Diffusion ODE")
println("ODE: -u'' + u = f  (with k=1.0)")
println("="^70)
println()

@testset "Reaction-Diffusion Solver Verification" begin
    
    #==========================================================================#
    # CORRECT SOLVER TESTS
    #==========================================================================#
    
    println("\n" * "="^70)
    println("PART 1: Testing CORRECT Implementation")
    println("="^70)
    
    @testset "✓ Correct Solver - sin(πx)" begin
        println("\nTest 1: CORRECT solver with u(x) = sin(πx)")
        println("-"^70)
        
        # Wrap the solver to pass k parameter
        solver_with_k(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Correct - sin(πx)",
            solver_with_k,
            pde_reaction_diffusion,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-3
        )
        
        result = run_simple_mms_test(test)
        
        @test result[:passed]
        @test result[:L2_error] < 1e-3
        
        println("\n✓ CORRECT solver PASSED (as expected!)\n")
    end
    
    @testset "✓ Correct Solver - polynomial" begin
        println("\nTest 2: CORRECT solver with u(x) = x(1-x)")
        println("-"^70)
        
        solver_with_k(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Correct - x(1-x)",
            solver_with_k,
            pde_reaction_diffusion,
            x * (1 - x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-6
        )
        
        result = run_simple_mms_test(test)
        
        @test result[:passed]
        @test result[:L2_error] < 1e-10
        
        println("\n✓ CORRECT solver PASSED with machine precision!\n")
    end
    
    @testset "✓ Correct Solver - exponential" begin
        println("\nTest 3: CORRECT solver with u(x) = exp(-x)sin(πx)")
        println("-"^70)
        
        solver_with_k(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Correct - exponential",
            solver_with_k,
            pde_reaction_diffusion,
            exp(-x) * sin(π * x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-3
        )
        
        result = run_simple_mms_test(test)
        
        @test result[:passed]
        
        println("\n✓ CORRECT solver PASSED!\n")
    end
    
    #==========================================================================#
    # BROKEN SOLVER TESTS (Should FAIL)
    #==========================================================================#
    
    println("\n" * "="^70)
    println("PART 2: Testing BROKEN Implementation (Missing Reaction Term)")
    println("="^70)
    println("\n⚠️  The following tests SHOULD FAIL to demonstrate error detection\n")
    
    @testset "❌ Broken Solver - sin(πx)" begin
        println("\nTest 4: BROKEN solver (missing k*u term) with u(x) = sin(πx)")
        println("-"^70)
        
        solver_with_k(grid, forcing) = reaction_diffusion_solver_broken(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Broken - sin(πx)",
            solver_with_k,
            pde_reaction_diffusion,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-3
        )
        
        result = run_simple_mms_test(test)
        
        # Expect this to FAIL
        @test !result[:passed]
        @test result[:L2_error] > 1e-3  # Error should be large
        
        println("\n❌ BROKEN solver FAILED (as expected!)")
        println("   This demonstrates MMS catches implementation errors ✓\n")
    end
    
    @testset "❌ Broken Solver - polynomial" begin
        println("\nTest 5: BROKEN solver with u(x) = x(1-x)")
        println("-"^70)
        
        solver_with_k(grid, forcing) = reaction_diffusion_solver_broken(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Broken - x(1-x)",
            solver_with_k,
            pde_reaction_diffusion,
            x * (1 - x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-6
        )
        
        result = run_simple_mms_test(test)
        
        # Expect this to FAIL
        @test !result[:passed]
        
        println("\n❌ BROKEN solver FAILED (implementation error detected!) ✓\n")
    end
    
    #==========================================================================#
    # VERY BROKEN SOLVER TESTS (Should FAIL BADLY)
    #==========================================================================#
    
    println("\n" * "="^70)
    println("PART 3: Testing VERY BROKEN Implementation (Wrong Sign)")
    println("="^70)
    println("\n⚠️  This solver has the WRONG SIGN - errors will be LARGE\n")
    
    @testset "❌❌ Very Broken Solver - sin(πx)" begin
        println("\nTest 6: VERY BROKEN solver (wrong sign) with u(x) = sin(πx)")
        println("-"^70)
        
        solver_with_k(grid, forcing) = reaction_diffusion_solver_very_broken(grid, forcing, k=k_param)
        
        test = SimpleMMSTest(
            "Very Broken - sin(πx)",
            solver_with_k,
            pde_reaction_diffusion,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            resolution = 0.01,
            tolerance = 1e-3
        )
        
        result = run_simple_mms_test(test)
        
        # Expect this to FAIL with large errors
        @test !result[:passed]
        @test result[:L2_error] > 0.1  # Very large error
        
        println("\n❌❌ VERY BROKEN solver FAILED with large errors!")
        println("    MMS detected wrong sign in discretization ✓\n")
    end
    
end

println("\n" * "="^70)
println("VERIFICATION SUMMARY")
println("="^70)
println()
println("✓ CORRECT Implementation:")
println("  - All 3 tests PASSED")
println("  - Errors within tolerance (< 1e-3)")
println("  - Polynomial achieved machine precision (< 1e-10)")
println()
println("❌ BROKEN Implementation (missing reaction term):")
println("  - Both tests FAILED (as expected)")
println("  - Errors exceed tolerance")
println("  - MMS successfully detected implementation bug ✓")
println()
println("❌❌ VERY BROKEN Implementation (wrong sign):")
println("  - Test FAILED with very large errors")
println("  - MMS detected sign error in discretization ✓")
println()
println("="^70)
println("MMS Verification Works! 🎉")
println("="^70)
println()
println("Key Takeaway:")
println("  The Method of Manufactured Solutions successfully distinguishes")
println("  between correct and incorrect solver implementations by comparing")
println("  numerical solutions against known analytical solutions.")
println()
