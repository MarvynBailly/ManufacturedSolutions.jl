# ManufacturedSolutions.jl

**A Julia package for verifying numerical PDE solvers using the Method of Manufactured Solutions (MMS).**

[![Project Status](https://img.shields.io/badge/status-active%20development-yellow.svg)](https://github.com/MarvynBailly/ManufacturedSolutions.jl)

## Overview

ManufacturedSolutions.jl provides a clean, automated framework for verifying numerical PDE solvers **written in any programming language** (Julia, C, Fortran, Python, etc.). This package helps you **verify your existing numerical methods** by:

1. **Symbolically computing forcing terms** from your PDE and manufactured solution
2. **Generating code** in your target language (C, Fortran, Python, Julia)
3. **Injecting forcing functions** into your solver source files
4. **Analyzing convergence rates** to verify correctness

---

## Table of Contents

- [Quick Start](#quick-start)
  - [Julia Solver Example](#julia-solver-example)
  - [C Code Injection Example](#c-code-injection-example)
- [Features](#features)
- [Installation](#installation)
- [Complete Usage Guide](#complete-usage-guide)
  - [Method of Manufactured Solutions](#method-of-manufactured-solutions)
  - [Symbolic Interface](#symbolic-interface)
  - [Code Injection Workflow](#code-injection-workflow)
- [API Reference](#api-reference)
- [Examples](#examples)
- [Project Status](#project-status)
- [References](#references)
- [License](#license)

---

## Quick Start

### Julia Solver Example

Verify a Julia solver with automatic forcing term computation:

```julia
using ManufacturedSolutions

# Step 1: Define PDE symbolically
@variables x u(x)
Dx = Differential(x)
pde_lhs = -Dx(Dx(u))  # Poisson: -u'' = f

# Step 2: Choose manufactured solution
u_mms_sym = sin(π * x)

# Step 3: Automatically compute forcing term!
forcing = compute_forcing_term_1d(pde_lhs, u_mms_sym, x, u)
# Result: f(x) = π²sin(πx) ≈ 9.8696sin(πx)  ✓ Computed automatically!

# Step 4: Verify your solver
u_mms_func = build_function(u_mms_sym, x, expression=Val{false})
u_mms = ManufacturedSolution(u_mms_func, (0.0, 1.0))

solver_wrapper(grid, mms) = my_poisson_solver(grid, forcing)
solver = UserSolver(solver_wrapper, name="My Solver")

results = verify_solver(solver, pde_lhs, u_mms, 
                       resolutions=[0.1, 0.05, 0.025])

# Step 5: Check convergence
println("L2 convergence rate: ", results.convergence_rates[:L2])
# Output: L2 convergence rate: 2.002  ✓ (2nd order accurate)
```

### C Code Injection Example

Verify a C solver by automatically generating and injecting forcing code:

```julia
using ManufacturedSolutions

# Step 1: Define PDE and manufactured solution
@variables x u(x)
Dx = Differential(x)
pde_lhs = -Dx(Dx(u))
u_mms = sin(π * x)

# Step 2: Compute forcing symbolically
forcing = compute_forcing_term_1d(pde_lhs, u_mms, x, u)

# Step 3: Generate C code
c_code = generate_forcing_code(forcing, :c, [:x])
# Generates: double forcing(double x) { return 9.8696*sin(M_PI*x); }

# Step 4: Inject into your C solver file
inject_forcing_code("my_solver.c", c_code, "my_solver_mms.c")
# ✓ Code injected at // MMS_FORCING_FUNCTION_HERE marker

# Step 5: Compile and run your C code
# gcc my_solver_mms.c -lm -o solver && ./solver
```

**Supported languages**: Julia, C, Python, Fortran

---

## Features

### ✅ Currently Available

#### Core Verification
- **Symbolic Forcing Computation**: Automatic derivation using Symbolics.jl - no manual calculus!
- **Multi-Language Support**: Generate code for C, Fortran, Python, Julia
- **Code Injection**: Automatically inject forcing functions into source files
- **1D Problem Support**: Full support for 1D PDEs (elliptic, parabolic, hyperbolic)
- **Error Norms**: L1, L2, L∞ with proper discrete formulations
- **Convergence Analysis**: Automatic order of accuracy estimation via log-log regression
- **Time-Dependent PDEs**: Support for parabolic/hyperbolic problems

#### User Experience
- **Clean API**: Simple, intuitive function interface
- **Comprehensive Results**: Timing, error data, convergence rates, R² goodness-of-fit
- **Flexible Integration**: Works with your existing solvers
- **Non-Intrusive**: Original files unchanged (creates new files)
- **Detailed Logging**: Optional verbose output for debugging

### 🔄 Coming Soon

- 2D and 3D grid support
- Curvilinear/elliptic grid support (for grid generators)
- More boundary condition types
- Power operations in code generation (x^2 → pow(x,2))
- Automatic compilation and execution
- Visualization tools
- Symbolic BC extraction

---

## Installation

```julia
# Clone the repository
git clone https://github.com/MarvynBailly/ManufacturedSolutions.jl

# Activate and instantiate
using Pkg
Pkg.activate("path/to/ManufacturedSolutions.jl")
Pkg.instantiate()
```

### Dependencies

- Julia 1.0+
- Symbolics.jl (for automatic differentiation)
- LinearAlgebra (stdlib)
- Printf (stdlib)
- Statistics (stdlib)

---

## Complete Usage Guide

### Method of Manufactured Solutions

MMS is a powerful code verification technique that works by:

1. **Choose** a manufactured solution `u(x)` that you know analytically
2. **Substitute** it into your PDE to compute the required forcing term `f(x)`
3. **Solve** your PDE with this forcing term using your numerical method
4. **Compare** the numerical solution to the known manufactured solution
5. **Analyze** convergence rates to verify your solver's order of accuracy

#### Example: Poisson Equation

**Given PDE**: `-u''(x) = f(x)` on [0,1]

**Choose**: `u(x) = sin(πx)`

**Compute forcing**:
- u'(x) = π cos(πx)
- u''(x) = -π² sin(πx)
- **Forcing**: f(x) = -(-π² sin(πx)) = π² sin(πx)

**Verify**: Your solver should reproduce u(x) = sin(πx) when given f(x) = π² sin(πx).

**ManufacturedSolutions.jl automates steps 2-5!**

---

### Symbolic Interface

The symbolic interface automatically computes forcing terms - no manual calculus needed!

#### 1D Spatial Problems (Elliptic PDEs)

```julia
using ManufacturedSolutions

@variables x u(x)
Dx = Differential(x)

# Define your PDE
pde = -Dx(Dx(u))  # Poisson equation

# Choose manufactured solution
u_mms = sin(π*x)

# Compute forcing automatically
forcing = compute_forcing_term_1d(pde, u_mms, x, u)

# Use the forcing:
println(forcing.symbolic_form)  # Display: 9.8696sin(πx)
f_value = forcing(0.5)          # Evaluate at x=0.5
```

#### Time-Dependent Problems (Parabolic/Hyperbolic PDEs)

```julia
@variables x t u(x,t)
Dx = Differential(x)
Dt = Differential(t)

# Heat equation: u_t - κu_xx = f
κ = 0.1
pde = Dt(u) - κ*Dx(Dx(u))

# Manufactured solution
u_mms = exp(-t) * sin(2π*x)

# Compute forcing automatically
forcing = compute_forcing_term_time_dependent(pde, u_mms, x, t, u)

# Use in solver:
f_value = forcing(0.5, 0.1)  # Evaluate at x=0.5, t=0.1
```

#### Benefits

✅ **No manual derivation** - Symbolics.jl does the calculus  
✅ **Fewer errors** - No mistakes in derivative computation  
✅ **Complex solutions** - Try solutions you wouldn't want to differentiate by hand  
✅ **Rapid iteration** - Change manufactured solution instantly  
✅ **Self-documenting** - Symbolic form shows exactly what forcing is used

---

### Code Injection Workflow

Verify solvers written in **any programming language** by generating and injecting forcing code.

#### Step 1: Prepare Your Solver

Add a marker comment where the forcing function should go:

**C Example** (`my_poisson.c`):
```c
#include <math.h>

// MMS_FORCING_FUNCTION_HERE
// Default implementation (will be replaced)
double forcing(double x) {
    return 0.0;
}

void solve_poisson(int n, double *x, double *u) {
    // Your solver implementation...
    for (int i = 0; i < n; i++) {
        double f = forcing(x[i]);
        // ... use f in your discretization ...
    }
}
```

**Python Example** (`my_solver.py`):
```python
import numpy as np

# MMS_FORCING_FUNCTION_HERE
def forcing(x):
    return 0.0

def solve(x):
    f = forcing(x)
    # ... your solver ...
```

**Fortran Example** (`solver.f90`):
```fortran
! MMS_FORCING_FUNCTION_HERE
function forcing(x) result(f)
    real(8) :: x, f
    f = 0.0d0
end function
```

#### Step 2: Generate and Inject Code

```julia
using ManufacturedSolutions

# Define problem
@variables x u(x)
Dx = Differential(x)
pde = -Dx(Dx(u))
u_mms = sin(π*x)

# Compute forcing
forcing = compute_forcing_term_1d(pde, u_mms, x, u)

# Generate C code
c_code = generate_forcing_code(forcing, :c, [:x])
# Produces:
# double forcing(double x) {
#     return 9.869604401089358*sin(M_PI*x);
# }

# Inject into your C file
inject_forcing_code(
    "my_poisson.c",       # Original file
    c_code,               # Generated code
    "my_poisson_mms.c"    # Output file
)
# ✓ Code injected successfully!
```

#### Step 3: Compile and Run

```bash
# Compile the modified solver
gcc my_poisson_mms.c -lm -o solver

# Run at different resolutions
./solver  # Outputs solution to file
```

#### Step 4: Analyze Results

```julia
using DelimitedFiles

# Read numerical solution
data = readdlm("solution.txt", skipstart=1)
x_num = data[:, 1]
u_num = data[:, 2]

# Compute analytical solution
u_exact = sin.(π .* x_num)

# Compute error
h = x_num[2] - x_num[1]
L2_error = sqrt(h * sum((u_num .- u_exact).^2))
println("L2 error: ", L2_error)

# Repeat for multiple resolutions and check convergence rate
```

#### Supported Languages

| Language | Template | Math Functions | Constants |
|----------|----------|----------------|-----------|
| **C** | `double forcing(double x)` | `sin`, `cos`, `exp` | `M_PI` |
| **Python** | `def forcing(x):` | `math.sin`, `math.cos` | `math.pi` |
| **Fortran** | `function forcing(x)` | `dsin`, `dcos`, `dexp` | `4.0d0*datan(1.0d0)` |
| **Julia** | `forcing(x)` | `sin`, `cos`, `exp` | `π` |

---

## API Reference

### Core Types

#### `Grid1D`
```julia
Grid1D(left, right, nx)   # Create with number of points
Grid1D(left, right, dx)   # Create with spacing
```
Uniform 1D grid with specified domain and resolution.

#### `ManufacturedSolution`
```julia
ManufacturedSolution(expression, domain; name="", description="")
```
Wrapper for analytical solution function.

#### `ForcingTerm`
```julia
forcing.symbolic_form  # Symbolic expression
forcing(x)            # Evaluate at point
forcing(x, t)         # Time-dependent evaluation
```
Callable forcing function with symbolic representation.

#### `UserSolver`
```julia
UserSolver(solve_function; name="", description="")
```
Wrapper for your solver function. Function signature:
```julia
function my_solver(grid::Grid1D, mms::ManufacturedSolution; kwargs...)
    # ... your implementation ...
    return numerical_solution::Vector{Float64}
end
```

#### `VerificationResult`
```julia
results.convergence_data   # Per-resolution error data
results.convergence_rates  # Dict{Symbol,Float64} with rates
results.problem_info       # Problem metadata
results.solver_info        # Solver metadata
```

### Main Functions

#### `verify_solver`
```julia
verify_solver(
    solver::UserSolver,
    pde,                          # Symbolic PDE (or nothing)
    manufactured_solution,
    resolutions::Vector;
    verbose::Bool = true,
    kwargs...
)
```
Main verification workflow. Returns `VerificationResult`.

#### Symbolic Forcing Computation

```julia
# 1D spatial problems (elliptic)
forcing = compute_forcing_term_1d(pde_lhs, u_mms_sym, x_var, u_var)

# Time-dependent problems (parabolic/hyperbolic)
forcing = compute_forcing_term_time_dependent(pde_lhs, u_mms_sym, x_var, t_var, u_var)
```

#### Code Generation

```julia
# Generate code in target language
code = generate_forcing_code(forcing, :c, [:x])        # C
code = generate_forcing_code(forcing, :python, [:x])   # Python
code = generate_forcing_code(forcing, :fortran, [:x])  # Fortran
code = generate_forcing_code(forcing, :julia, [:x])    # Julia

# For time-dependent:
code = generate_forcing_code(forcing, :c, [:x, :t])
```

#### Code Injection

```julia
marker_found = inject_forcing_code(
    source_file,      # Original file path
    forcing_code,     # Generated code string
    output_file;      # Where to write modified file
    marker = "MMS_FORCING"  # Optional: custom marker
)
```

#### Error Norms

```julia
L1 = compute_L1_error(numerical, analytical, grid)
L2 = compute_L2_error(numerical, analytical, grid)
Linf = compute_Linf_error(numerical, analytical)
errors = compute_all_errors(numerical, analytical, grid)
```

#### Convergence Analysis

```julia
rate, R2 = estimate_convergence_rate(errors, resolutions)
rates = analyze_convergence(convergence_data)
```

---

## Examples

The `examples/` directory contains complete working examples:

### Julia Solvers

```bash
# Basic Poisson solver (manual forcing)
julia --project=. examples/01_simple_poisson_1d.jl
# Expected output: L2 rate ≈ 2.0 (2nd order convergence) ✓

# Heat equation solver (manual forcing)
julia --project=. examples/02_heat_equation_1d.jl
# Expected output: L2 rate ≈ 1.0 (Forward Euler, 1st order in time) ✓

# Poisson with symbolic forcing
julia --project=. examples/03_symbolic_poisson_1d.jl
# Demonstrates automatic forcing computation ✓

# Heat equation with symbolic forcing
julia --project=. examples/04_symbolic_heat_1d.jl
# Demonstrates time-dependent symbolic forcing ✓
```

### C Code Injection

```bash
# Code injection demonstration
julia --project=. examples/05_c_code_injection.jl
# Shows complete workflow: symbolic → generate C → inject ✓
```

### Expected Results

All examples should show convergence rates matching the theoretical order:

**Poisson Solver** (2nd order spatial):
```
L2    Rate: 2.002  (R² = 1.0000)  ✓
L1    Rate: 1.998  (R² = 1.0000)  ✓
Linf  Rate: 2.002  (R² = 1.0000)  ✓
```

**Heat Equation** (1st order in time, Forward Euler):
```
L2    Rate: 0.993  (R² = 1.0000)  ✓
L1    Rate: 0.998  (R² = 1.0000)  ✓
Linf  Rate: 0.937  (R² = 0.9995)  ✓
```

---

## Project Status

### ✅ Phase 1 Complete: 1D Foundation
- Core type system
- Error norms and convergence analysis
- Working Julia solver verification
- Comprehensive documentation

### ✅ Phase 2 Complete: Symbolic Interface
- Automatic forcing term computation
- Symbolics.jl integration
- 1D spatial and time-dependent support
- Working symbolic examples

### 🔄 Phase 3 In Progress: Code Injection
- Multi-language code generation (C, Python, Fortran, Julia)
- File injection system with markers
- C Poisson solver example
- **Known issues**: Duplicate function definitions (minor), power operations not yet supported

### ⏳ Upcoming: 2D Support and Elliptic Grids
- 2D grid types
- 2D symbolic forcing
- Curvilinear coordinate support
- Elliptic grid generator verification (primary goal!)

See `DEVELOPMENT_ROADMAP.md` for the complete 10-phase plan.

---

## Use Cases

This package is designed for:

- ✅ **Verifying elliptic grid generators** (primary goal)
- ✅ Verifying custom PDE solvers in any language
- ✅ Testing numerical method implementations
- ✅ Convergence studies for research papers
- ✅ Code verification for publications
- ✅ Teaching numerical methods (students can verify their solvers)

---

## Contributing

Contributions welcome! This is an active research project.

**Current priorities:**
1. Fix duplicate function issue in code injection
2. Add power operation support (x^2 → pow(x,2))
3. 2D grid support
4. Elliptic grid examples

---

## References

### Method of Manufactured Solutions
- Roache, P.J. (2002): "Code Verification by the Method of Manufactured Solutions", *Journal of Fluids Engineering*
- Salari, K. & Knupp, P. (2000): "Code Verification by the Method of Manufactured Solutions", Sandia National Labs Report SAND2000-1444

### Symbolic Computation
- Gowda, S. et al. (2021): "High-performance symbolic-numerics via multiple dispatch", *arXiv:2105.03949* (Symbolics.jl)

---

## License

See LICENSE file.

---

## Acknowledgments

Developed as part of research on elliptic grid generation and PDE solver verification. Based on prototypes in `example/MMS/`.

**Author**: MarvynBailly  
**Status**: Active Development  
**Version**: 0.3.0-dev (Code Injection Alpha)
