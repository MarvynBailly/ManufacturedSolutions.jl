"""
Example Poisson solver to be tested with MMS unit tests

Solves: -u'' = f(x) on [0,1] with Dirichlet BCs u(0)=u(1)=0

This file demonstrates the required structure for MMS unit testing:
1. Place a marker comment where the forcing function should be injected
2. Read resolution from ENV["MMS_TEST_RESOLUTION"]
3. Write output to ENV["MMS_TEST_OUTPUT"]
"""

using LinearAlgebra
using Printf

# MMS_FORCING_FUNCTION_HERE

"""
Solve -u'' = f(x) using finite differences
"""
function solve_poisson_1d(nx, forcing_func)
    # Domain
    x_left = 0.0
    x_right = 1.0
    h = (x_right - x_left) / (nx - 1)
    
    # Grid
    x = range(x_left, x_right, length=nx)
    
    # Allocate solution
    u = zeros(nx)
    
    # Build tridiagonal system: -u_{i-1} + 2u_i - u_{i+1} = h^2 * f_i
    A = zeros(nx, nx)
    rhs = zeros(nx)
    
    for i in 1:nx
        if i == 1
            # Left BC: u(0) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        elseif i == nx
            # Right BC: u(1) = 0
            A[i, i] = 1.0
            rhs[i] = 0.0
        else
            # Interior: -u'' = f
            A[i, i-1] = -1.0
            A[i, i] = 2.0
            A[i, i+1] = -1.0
            rhs[i] = h^2 * forcing_func(x[i])
        end
    end
    
    # Solve
    u = A \ rhs
    
    return collect(x), u
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # Get resolution from environment (set by MMS test)
    h = parse(Float64, get(ENV, "MMS_TEST_RESOLUTION", "0.01"))
    output_file = get(ENV, "MMS_TEST_OUTPUT", "mms_test_output.txt")
    
    # Compute grid size
    nx = Int(round(1.0 / h)) + 1
    
    println("Solving Poisson equation with h = $h (nx = $nx)")
    
    # Solve
    x, u = solve_poisson_1d(nx, forcing_term)
    
    # Write output
    open(output_file, "w") do io
        println(io, "# x u")
        for i in 1:length(x)
            @printf(io, "%.15e %.15e\n", x[i], u[i])
        end
    end
    
    println("Solution written to: $output_file")
    println("Max |u|: $(maximum(abs.(u)))")
end
