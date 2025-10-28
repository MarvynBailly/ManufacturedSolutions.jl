# MMS Unit Testing Workflow - Quick Guide

## Overview

The MMS unit testing workflow provides a simple way to verify numerical solvers with a single function call. No convergence studies, no manual forcing - just verify correctness!

## Quick Start

### 1. Prepare Your Solver

Add a marker comment where the forcing function should be injected:

```julia
# MMS_FORCING_FUNCTION_HERE

function solve_poisson_1d(nx, forcing_func)
    # Your solver implementation
    # Use: f_value = forcing_func(x)
end
```

**Important**: 
- Your solver should read resolution from `ENV["MMS_TEST_RESOLUTION"]`
- Your solver should write output to `ENV["MMS_TEST_OUTPUT"]`
- Output format: text file with columns `x u_numerical`

### 2. Create a Unit Test

```julia
using ManufacturedSolutions

# Define PDE and manufactured solution
@variables x u(x)
Dx = Differential(x)

test = MMSUnitTest(
    "My Solver Test",
    "my_solver.jl",
    -Dx(Dx(u)),           # PDE
    sin(π*x),             # Manufactured solution
    (x, u),               # Variables
    (0.0, 1.0),           # Domain
    resolution = 0.01,
    tolerance = 1e-3
)

# Run the test
result = run_mms_unit_test(test)

# Check result
if result[:passed]
    println("✓ PASSED")
else
    println("✗ FAILED")
end
```

### 3. What Happens Automatically

The tool will:
1. **Compute forcing term** symbolically: f(x) = π²sin(πx)
2. **Generate Julia code** for the forcing function
3. **Inject** the code into your solver (creates `my_solver_mms_test.jl`)
4. **Run** your solver
5. **Read** the output
6. **Compare** against manufactured solution
7. **Return** pass/fail with error metrics

## Solver Template

Your solver needs this structure:

```julia
#!/usr/bin/env julia

using LinearAlgebra

# MMS_FORCING_FUNCTION_HERE

function my_solver(...)
    # Your implementation
    # Use forcing_term(x) to evaluate forcing
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # Get parameters from environment
    h = parse(Float64, get(ENV, "MMS_TEST_RESOLUTION", "0.01"))
    output_file = get(ENV, "MMS_TEST_OUTPUT", "output.txt")
    
    # Solve
    x, u = my_solver(h, forcing_term)
    
    # Write output
    open(output_file, "w") do io
        println(io, "# x u")
        for i in 1:length(x)
            println(io, "$(x[i]) $(u[i])")
        end
    end
end
```

## Result Structure

The `run_mms_unit_test` function returns a dictionary:

```julia
result = Dict(
    :passed => true/false,
    :L2_error => 5.8e-5,
    :Linf_error => 8.2e-5,
    :message => "L2 error (5.816e-5) < tolerance (0.001)",
    :num_points => 101,
    :resolution => 0.01
)
```

## Integration with Julia Test Framework

```julia
using Test
using ManufacturedSolutions

@testset "My Solver Tests" begin
    @variables x u(x)
    Dx = Differential(x)
    
    test = MMSUnitTest(
        "Poisson Test",
        "my_solver.jl",
        -Dx(Dx(u)),
        sin(π*x),
        (x, u),
        (0.0, 1.0)
    )
    
    result = run_mms_unit_test(test, verbose=false)
    @test result[:passed]
    @test result[:L2_error] < 1e-3
end
```

## Example Output

```
======================================================================
MMS Unit Test: Poisson 1D Test
======================================================================

[1/6] Computing forcing term symbolically...
  Forcing: 9.869604401089358sin(πx)

[2/6] Generating Julia forcing function...
  Generated function:
    function forcing_term(x)
      return 9.869604401089358*sin(π*x)
    end

[3/6] Injecting forcing into source file...
  ✓ Code injected successfully

[4/6] Running solver...
  Solver output:
    Solving Poisson equation with h = 0.01 (nx = 101)
    Max |u|: 1.0000822507622098

[5/6] Reading solver output...
  Read 101 solution points

[6/6] Verifying against manufactured solution...
  L2 error:   5.816007e-05
  L∞ error:   8.225076e-05
  Tolerance:  1.000000e-03

  ✓ TEST PASSED

🎉 Solver verified successfully!
```

## Advantages Over Traditional MMS

**Traditional MMS workflow**:
- Manually compute forcing term
- Manually code forcing function
- Run solver at multiple resolutions
- Manually compute errors
- Manually fit convergence rates
- ~50-100 lines of test code

**MMS Unit Test workflow**:
- Everything automatic
- Single resolution test
- Pass/fail result
- ~10 lines of test code

## When to Use Unit Tests vs Convergence Studies

**Use Unit Tests When:**
- You want quick verification during development
- You're doing test-driven development
- You want CI/CD integration
- You just need to know if the solver is "correct enough"

**Use Convergence Studies When:**
- You need to verify the order of accuracy
- You're writing a paper and need convergence plots
- You're comparing multiple numerical methods
- You need detailed error analysis

Both workflows are available in ManufacturedSolutions.jl!

## Common Issues

### Issue: Test fails with large error

**Cause**: Solver might not be using the injected forcing function

**Solution**: Check that `forcing_term(x)` is called in your solver

### Issue: Output file not found

**Cause**: Solver not reading `ENV["MMS_TEST_OUTPUT"]`

**Solution**: Make sure your solver writes to the specified output file

### Issue: Marker not found

**Cause**: Wrong marker in source file

**Solution**: Use exactly `# MMS_FORCING_FUNCTION_HERE`

## See Also

- `examples/06_unit_test_poisson.jl` - Complete working example
- `examples/solvers/poisson_solver_1d.jl` - Example solver template
- Main README.md - Full package documentation
