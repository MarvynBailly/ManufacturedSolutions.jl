"""
Example: MMS Unit Test for Poisson Solver

This demonstrates the simplified unit testing workflow:
1. Define PDE and manufactured solution
2. Run MMS unit test
3. Get pass/fail result

No manual forcing computation, no convergence studies - just verify correctness!
"""

using ManufacturedSolutions

println("="^70)
println("MMS Unit Test: Poisson Solver Verification")
println("="^70)
println()

# Define the test symbolically
@variables x u(x)
Dx = Differential(x)

# Create the test
test = MMSUnitTest(
    "Poisson 1D Test",                          # Test name
    "examples/solvers/poisson_solver_1d.jl",   # Solver to test
    -Dx(Dx(u)),                                 # PDE: -u'' = f
    sin(π * x),                                 # Manufactured solution
    (x, u),                                     # Variables
    (0.0, 1.0),                                 # Domain
    resolution = 0.01,                          # Grid spacing
    output_file = "test_poisson_output.txt",   # Where solver writes results
    tolerance = 1e-3                            # Pass/fail threshold
)

println("Test Configuration:")
println("  Solver: $(test.source_file)")
println("  PDE: -u'' = f")
println("  Manufactured solution: sin(πx)")
println("  Resolution: h = $(test.resolution)")
println("  Tolerance: $(test.tolerance)")
println()

# Run the test
result = run_mms_unit_test(test, verbose=true)

# Print summary
println()
println("="^70)
println("Test Summary")
println("="^70)
println("Test name:    $(test.test_name)")
println("Status:       ", result[:passed] ? "✓ PASSED" : "✗ FAILED")
println("L2 error:     ", result[:L2_error])
println("L∞ error:     ", result[:Linf_error])
println("Grid points:  ", result[:num_points])
println("Message:      ", result[:message])
println()

if result[:passed]
    println("🎉 Solver verified successfully!")
else
    println("⚠️  Solver failed verification")
    exit(1)
end
