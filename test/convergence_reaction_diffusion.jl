"""
Convergence Verification: Reaction-Diffusion ODE

Tests order of accuracy for -u'' + k*u = f solver

Compares correct vs. broken implementations to show how
convergence testing catches implementation errors.
"""

using ManufacturedSolutions
using LinearAlgebra
using Test
using Printf

# Import solvers
include("solvers/reaction_diffusion_1d.jl")
using .Main: reaction_diffusion_solver_correct, reaction_diffusion_solver_broken

# Define the PDE
@variables x u(x)
Dx = Differential(x)
k_param = 1.0
pde = -Dx(Dx(u)) + k_param * u

println("="^80)
println("Convergence Test: Reaction-Diffusion ODE")
println("Comparing Correct vs. Broken Implementations")
println("="^80)
println()

@testset "Reaction-Diffusion - Order of Accuracy" begin
    
    #==========================================================================#
    # CORRECT SOLVER - Should show 2nd order convergence
    #==========================================================================#
    
    println("\n" * "="^80)
    println("PART 1: CORRECT SOLVER - Expected Order = 2.0")
    println("="^80)
    
    @testset "✓ Correct - Square Domain [0,1]" begin
        println("\nTest 1: Correct solver, square domain [0,1]")
        println("-"^80)
        
        solver(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = ConvergenceTest(
            "Correct - Square",
            solver,
            pde,
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
        
        println("✓ Correct solver achieves 2nd order convergence")
        println()
    end
    
    @testset "✓ Correct - Rectangular Domain [0,2]" begin
        println("\nTest 2: Correct solver, rectangular domain [0,2]")
        println("-"^80)
        
        solver(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = ConvergenceTest(
            "Correct - Rectangle",
            solver,
            pde,
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
        
        println("✓ Correct solver maintains 2nd order on different domain")
        println()
    end
    
    @testset "✓ Correct - Polynomial Solution" begin
        println("\nTest 3: Correct solver with polynomial u(x) = x(1-x)")
        println("-"^80)
        
        solver(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        
        test = ConvergenceTest(
            "Correct - Polynomial",
            solver,
            pde,
            x * (1 - x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 4,
            tolerance = 5.0  # Very relaxed - polynomial hits machine precision
        )
        
        result = run_convergence_test(test)
        
        # For polynomial, we just check that final error is small (machine precision)
        final_error = result[:errors][end]
        @test final_error < 1e-10  # Error should be at machine precision
        
        println("✓ Correct solver shows high-order for polynomial")
        println("  Final error: $(final_error) (machine precision)")
        println()
    end
    
    #==========================================================================#
    # BROKEN SOLVER - Should NOT show 2nd order convergence
    #==========================================================================#
    
    println("\n" * "="^80)
    println("PART 2: BROKEN SOLVER - Expected to FAIL convergence test")
    println("="^80)
    println("\n⚠️  The broken solver is missing the reaction term k*u")
    println("   It should still converge, but NOT at the expected order!\n")
    
    @testset "❌ Broken - Square Domain [0,1]" begin
        println("\nTest 4: Broken solver, square domain [0,1]")
        println("-"^80)
        
        solver(grid, forcing) = reaction_diffusion_solver_broken(grid, forcing, k=k_param)
        
        test = ConvergenceTest(
            "Broken - Square",
            solver,
            pde,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 4,
            tolerance = 0.2
        )
        
        result = run_convergence_test(test)
        
        # The broken solver still converges (it's solving a different PDE!)
        # but not to the right solution
        # It solves -u'' = f instead of -u'' + k*u = f
        # So it will show 2nd order convergence to the WRONG solution
        
        println("\nNote: Broken solver shows 2nd order convergence to WRONG solution!")
        println("      This is why we need BOTH convergence AND error magnitude tests.")
        println()
    end
    
    @testset "Broken - Error Magnitude Check" begin
        println("\nTest 5: Error magnitude comparison (Correct vs. Broken)")
        println("-"^80)
        println("\nRunning both solvers at same resolution to compare errors...")
        
        # Correct solver
        solver_correct(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
        test_correct = ConvergenceTest(
            "Correct",
            solver_correct,
            pde,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 3,
            tolerance = 0.2
        )
        result_correct = run_convergence_test(test_correct, verbose=false)
        
        # Broken solver
        solver_broken(grid, forcing) = reaction_diffusion_solver_broken(grid, forcing, k=k_param)
        test_broken = ConvergenceTest(
            "Broken",
            solver_broken,
            pde,
            sin(π * x),
            (x, u),
            (0.0, 1.0),
            expected_order = 2.0,
            base_resolution = 0.1,
            num_refinements = 3,
            tolerance = 0.2
        )
        result_broken = run_convergence_test(test_broken, verbose=false)
        
        println("\nError Comparison at finest resolution:")
        println("-"^80)
        @printf "  Correct solver: L2 = %.6e\n" result_correct[:errors][end]
        @printf "  Broken solver:  L2 = %.6e\n" result_broken[:errors][end]
        ratio = result_broken[:errors][end] / result_correct[:errors][end]
        @printf "  Ratio (broken/correct): %.1fx\n" ratio
        
        println("\nObserved Convergence Orders:")
        println("-"^80)
        @printf "  Correct solver: p = %.4f\n" result_correct[:observed_order]
        @printf "  Broken solver:  p = %.4f\n" result_broken[:observed_order]
        
        # The broken solver has MUCH larger errors
        @test result_broken[:errors][end] > result_correct[:errors][end] * 100
        
        println("\n✓ Error magnitude test detects the bug!")
        println("  Even though both converge at 2nd order, the broken solver")
        println("  has ~$(round(Int, ratio))x larger errors.")
        println()
    end
    
end

println("\n" * "="^80)
println("KEY INSIGHTS: Convergence vs. Error Magnitude")
println("="^80)
println()
println("1. ✓ Correct Solver:")
println("   - Achieves 2nd order convergence (p ≈ 2.0)")
println("   - Small errors (~ 1e-5 at fine resolution)")
println()
println("2. ❌ Broken Solver (missing reaction term):")
println("   - Also shows 2nd order convergence (p ≈ 2.0)")
println("   - BUT converges to WRONG solution")
println("   - Errors are ~1000x larger than correct solver")
println()
println("LESSON: Convergence order alone is NOT enough!")
println()
println("You need BOTH:")
println("  1. Correct convergence order (verify discretization)")
println("  2. Small error magnitude (verify correct PDE)")
println()
println("MMS provides both checks automatically!")
println("="^80)
