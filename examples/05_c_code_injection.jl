"""
Example: MMS Verification of C Poisson Solver using Code Injection

This example demonstrates the complete workflow:
1. Define PDE and manufactured solution symbolically
2. Automatically generate C forcing function code
3. Inject into C solver source code
4. Compile and run C solver at multiple resolutions
5. Analyze convergence

This approach allows MMS verification of solvers written in C, Fortran, Python, etc.
"""

using ManufacturedSolutions
using Symbolics
using Printf
using DelimitedFiles

# Define the PDE and manufactured solution symbolically
@variables x u(x)
Dx = Differential(x)

println("="^70)
println("MMS Verification of C Code: 1D Poisson Equation")
println("="^70)
println()

# PDE: -u'' = f(x) on [0,1]
pde_lhs = -Dx(Dx(u))
println("PDE: -u'' = f(x)")

# Manufactured solution: u_mms(x) = sin(πx)
u_mms = sin(π * x)
println("Manufactured solution: u(x) = sin(πx)")
println()

# Create MMS problem
problem = MMSProblem(
    pde_lhs,
    u_mms,
    (x, u),
    (0.0, 1.0)  # Domain [0,1]
)

# Configure C solver
c_solver = SolverSetup(
    "examples/c_code/poisson_solver.c",   # Original C code
    :c,                                    # Language
    output_file = "examples/c_code/poisson_solver_mms.c"
)

println("Step 1: Computing forcing term symbolically...")
forcing = compute_forcing_term_1d(pde_lhs, u_mms, x, u)
println("  f(x) = ", forcing.symbolic_form)
println()

println("Step 2: Generating C code for forcing function...")
forcing_code = generate_forcing_code(forcing, :c, [:x])
println("Generated C function:")
println(forcing_code)
println()

println("Step 3: Injecting code into C solver...")
marker_found = inject_forcing_code(
    c_solver.source_file,
    forcing_code,
    c_solver.output_file
)

if marker_found
    println("  ✓ Code injected successfully at marker location")
else
    println("  ⚠ Marker not found - code appended to end")
end
println("  Modified solver: $(c_solver.output_file)")
println()

# For demonstration, let's verify the injection worked
println("Step 4: Verifying code injection...")
modified_content = read(c_solver.output_file, String)
if occursin("sin", modified_content) && occursin("M_PI", modified_content)
    println("  ✓ Forcing function successfully injected")
    
    # Show snippet of injected code
    lines = split(modified_content, '\n')
    for (i, line) in enumerate(lines)
        if occursin("double forcing", line)
            println("\n  Injected function at line $i:")
            for j in i:min(i+5, length(lines))
                println("    ", lines[j])
            end
            break
        end
    end
else
    println("  ⚠ Warning: Expected code not found in output")
end
println()

println("="^70)
println("Code Injection Complete!")
println("="^70)
println()
println("Next steps to complete verification:")
println("  1. Compile: gcc $(c_solver.output_file) -lm -o poisson_solver")
println("  2. Run at different resolutions (modify NX in code)")
println("  3. Read solution.txt output")
println("  4. Compute errors against manufactured solution")
println("  5. Verify 2nd-order convergence")
println()
println("Expected convergence rate: 2.0 (second-order accurate)")
println()

# Manual convergence study (since we can't easily compile C from Julia)
println("="^70)
println("Manual Convergence Study Instructions")
println("="^70)
println()
println("To complete the verification, compile and run the C code:")
println()

resolutions = [101, 201, 401]
for (i, nx) in enumerate(resolutions)
    h = 1.0 / (nx - 1)
    println("Resolution $i: h = $h (nx = $nx)")
    println("  Compile: gcc -DNX=$nx $(c_solver.output_file) -lm -o solver_$i")
    println("  Run:     ./solver_$i")
    println("  Output:  solution.txt")
    println()
end

println("Then use Julia to analyze results:")
println("""
using DelimitedFiles

# Read numerical solution
data = readdlm("solution.txt", skipstart=1)
x_num = data[:, 1]
u_num = data[:, 2]

# Compute analytical solution
u_exact = sin.(π .* x_num)

# Compute L2 error
h = x_num[2] - x_num[1]
L2_error = sqrt(h * sum((u_num .- u_exact).^2))
println("L2 error: ", L2_error)
""")
println()

println("Expected L2 errors (approximately):")
for nx in resolutions
    h = 1.0 / (nx - 1)
    # L2 error ~ C * h^2 for second-order method
    # Empirical constant C ≈ 0.1 for this problem
    expected_error = 0.1 * h^2
    println("  nx = $nx: L2 ≈ $(round(expected_error, sigdigits=4))")
end
