"""
Quick demonstration: Run all three solvers side-by-side

This script runs the correct, broken, and very broken solvers
on the same test case to show the error differences.
"""

using ManufacturedSolutions
using Printf

include("solvers/reaction_diffusion_1d.jl")
using .Main: reaction_diffusion_solver_correct, reaction_diffusion_solver_broken,
             reaction_diffusion_solver_very_broken

# Define the PDE
@variables x u(x)
Dx = Differential(x)
k_param = 1.0
pde = -Dx(Dx(u)) + k_param * u

# Manufactured solution
manufactured_sol = sin(π * x)

println("="^80)
println("Side-by-Side Comparison: Correct vs. Broken Solvers")
println("="^80)
println("\nODE: -u'' + u = f  with u(0) = u(1) = 0")
println("Manufactured Solution: u(x) = sin(πx)")
println("Domain: [0,1] with h = 0.01")
println()

# Test parameters
domain = (0.0, 1.0)
res = 0.01
tol = 1e-3

println("="^80)
println("TEST 1: ✅ CORRECT SOLVER")
println("="^80)

solver1(grid, forcing) = reaction_diffusion_solver_correct(grid, forcing, k=k_param)
test1 = SimpleMMSTest("Correct", solver1, pde, manufactured_sol, (x, u), domain, 
                      resolution=res, tolerance=tol)
result1 = run_simple_mms_test(test1)

println("\n" * "="^80)
println("TEST 2: ❌ BROKEN SOLVER (Missing Reaction Term)")
println("="^80)

solver2(grid, forcing) = reaction_diffusion_solver_broken(grid, forcing, k=k_param)
test2 = SimpleMMSTest("Broken", solver2, pde, manufactured_sol, (x, u), domain, 
                      resolution=res, tolerance=tol)
result2 = run_simple_mms_test(test2)

println("\n" * "="^80)
println("TEST 3: ❌❌ VERY BROKEN SOLVER (Wrong Sign)")
println("="^80)

solver3(grid, forcing) = reaction_diffusion_solver_very_broken(grid, forcing, k=k_param)
test3 = SimpleMMSTest("Very Broken", solver3, pde, manufactured_sol, (x, u), domain, 
                      resolution=res, tolerance=tol)
result3 = run_simple_mms_test(test3)

# Summary comparison
println("\n" * "="^80)
println("SUMMARY COMPARISON")
println("="^80)
println()

function print_comparison_row(label, correct_val, broken_val, very_broken_val)
    ratio1 = broken_val / correct_val
    ratio2 = very_broken_val / correct_val
    @printf("%-20s | %12.2e | %12.2e (%5.0fx) | %12.2e (%5.0fx)\n", 
            label, correct_val, broken_val, ratio1, very_broken_val, ratio2)
end

println("Error Metric         | Correct      | Broken (ratio)      | Very Broken (ratio)")
println("-"^80)
print_comparison_row("L1 Error", result1[:L1_error], result2[:L1_error], result3[:L1_error])
print_comparison_row("L2 Error", result1[:L2_error], result2[:L2_error], result3[:L2_error])
print_comparison_row("L∞ Error", result1[:Linf_error], result2[:Linf_error], result3[:Linf_error])
println("-"^80)
println()

println("Test Status:")
println("  ✅ Correct:      $(result1[:passed] ? "PASS" : "FAIL")")
println("  ❌ Broken:       $(result2[:passed] ? "PASS" : "FAIL") (missing k*u term)")
println("  ❌❌ Very Broken: $(result3[:passed] ? "PASS" : "FAIL") (wrong sign)")
println()

println("Key Observations:")
println("  1. Broken solver has ~$(round(Int, result2[:L2_error]/result1[:L2_error]))x larger error than correct solver")
println("  2. Very broken solver has ~$(round(Int, result3[:L2_error]/result1[:L2_error]))x larger error than correct solver")
println("  3. MMS verification clearly distinguishes correct from incorrect implementations")
println()

println("="^80)
println("Conclusion: MMS is an effective automatic verification tool! ✓")
println("="^80)
