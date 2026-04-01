# ManufacturedSolutions.jl# ManufacturedSolutions.jl



> **A Julia package for verifying numerical PDE solvers using the Method of Manufactured Solutions (MMS)**> **A Julia package for verifying numerical PDE solvers using the Method of Manufactured Solutions (MMS)**



[![Project Status](https://img.shields.io/badge/status-active%20development-yellow.svg)](https://github.com/MarvynBailly/ManufacturedSolutions.jl)[![Project Status](https://img.shields.io/badge/status-active%20development-yellow.svg)](https://github.com/MarvynBailly/ManufacturedSolutions.jl)



## What is This?## What is This?



This package helps you verify that your numerical PDE solver is implemented correctly by:This package helps you verify that your numerical PDE solver is implemented correctly by:



1. **You provide**: A PDE and a manufactured solution

2. **We compute**: The forcing term needed (automatically, using symbolic math)

3. **You run**: Your solver with this forcing1. **You provide**: A PDE and a manufactured solutionManufacturedSolutions.jl provides a clean, automated framework for verifying numerical PDE solvers **written in any programming language** (Julia, C, Fortran, Python, etc.). This package helps you **verify your existing numerical methods** by:

4. **We check**: Whether your solver produces the correct answer and convergence rate

2. **We compute**: The forcing term needed (automatically, using symbolic math)

**No manual calculus required!** The package does all the symbolic differentiation for you.

3. **You run**: Your solver with this forcing1. **Symbolically computing forcing terms** from your PDE and manufactured solution

## Quick Example

4. **We check**: Whether your solver produces the correct answer and convergence rate2. **Generating code** in your target language (C, Fortran, Python, Julia)

```julia

using ManufacturedSolutions3. **Injecting forcing functions** into your solver source files



# 1. Define your PDE symbolically (Poisson equation: -u'' = f)**No manual calculus required!** The package does all the symbolic differentiation for you.4. **Analyzing convergence rates** to verify correctness

@variables x u(x)

Dx = Differential(x)

pde = -Dx(Dx(u))

## Quick Example---

# 2. Choose a manufactured solution

u_manufactured = sin(π * x)



# 3. Automatically compute what forcing term is needed```julia## Table of Contents

forcing = compute_forcing_term_1d(pde, u_manufactured, x, u)

# Result: f(x) = π² sin(πx) ≈ 9.87 sin(πx)  ✓using ManufacturedSolutions



# 4. Test your solver at multiple resolutions- [Quick Start](#quick-start)

test = ConvergenceTest(

    "My Poisson Solver Test",# 1. Define your PDE symbolically (Poisson equation: -u'' = f)  - [Julia Solver Example](#julia-solver-example)

    my_solver_function,

    pde,@variables x u(x)  - [C Code Injection Example](#c-code-injection-example)

    u_manufactured,

    (x, u),Dx = Differential(x)- [Features](#features)

    (0.0, 1.0),              # domain [0,1]

    expected_order = 2.0,     # expect 2nd order accuracypde = -Dx(Dx(u))- [Installation](#installation)

    base_resolution = 0.1,    # coarsest grid spacing

    num_refinements = 4       # test 4 different resolutions- [Complete Usage Guide](#complete-usage-guide)

)

# 2. Choose a manufactured solution  - [Method of Manufactured Solutions](#method-of-manufactured-solutions)

result = run_convergence_test(test)

# Output: ✓ Observed order 2.002 matches expected order 2.0u_manufactured = sin(π * x)  - [Symbolic Interface](#symbolic-interface)

```

  - [Code Injection Workflow](#code-injection-workflow)

## Installation

# 3. Automatically compute what forcing term is needed- [API Reference](#api-reference)

```julia

using Pkgforcing = compute_forcing_term_1d(pde, u_manufactured, x, u)- [Examples](#examples)

Pkg.add(url="https://github.com/MarvynBailly/ManufacturedSolutions.jl")

```# Result: f(x) = π² sin(πx) ≈ 9.87 sin(πx)  ✓- [Project Status](#project-status)



Or clone and develop locally:- [References](#references)



```bash# 4. Test your solver at multiple resolutions- [License](#license)

git clone https://github.com/MarvynBailly/ManufacturedSolutions.jl

cd ManufacturedSolutions.jltest = ConvergenceTest(

julia --project=. -e 'using Pkg; Pkg.instantiate()'

```    "My Poisson Solver Test",---



## What is the Method of Manufactured Solutions?    my_solver_function,



The Method of Manufactured Solutions (MMS) is a powerful technique for verifying that your numerical code is implemented correctly:    pde,## Quick Start



**The Problem**: How do you know your PDE solver works if you don't know the exact solution?    u_manufactured,



**The Solution**: Work backwards!    (x, u),### Julia Solver Example



1. Choose a solution you **do** know (the "manufactured solution")    (0.0, 1.0),              # domain [0,1]

2. Plug it into your PDE to find what forcing term would produce it

3. Solve your PDE with that forcing term    expected_order = 2.0,     # expect 2nd order accuracyVerify a Julia solver with automatic forcing term computation:

4. Check if you get back the solution you started with

    base_resolution = 0.1,    # coarsest grid spacing

**Example**: For the Poisson equation `-u'' = f`:

- Choose: `u(x) = sin(πx)` (we know this exactly!)    num_refinements = 4       # test 4 different resolutions```julia

- Compute: `u''(x) = -π² sin(πx)`, so we need `f(x) = π² sin(πx)`

- Solve: Run your solver with `f(x) = π² sin(πx)`)using ManufacturedSolutions

- Verify: Does your numerical solution match `sin(πx)`?



If your solver is correct, the numerical and analytical solutions should match (up to discretization error), and the error should decrease at the expected rate as you refine the grid.

result = run_convergence_test(test)# Step 1: Define PDE symbolically

## Features

# Output: ✓ Observed order 2.002 matches expected order 2.0@variables x u(x)

### Verification Methods

```Dx = Differential(x)

**Simple Unit Tests** - Test single cases:

```juliapde_lhs = -Dx(Dx(u))  # Poisson: -u'' = f

test = SimpleMMSTest(

    "Test name",## Installation

    solver_function,

    pde,# Step 2: Choose manufactured solution

    manufactured_solution,

    (x, u),```juliau_mms_sym = sin(π * x)

    (0.0, 1.0),

    resolution = 0.01,using Pkg

    tolerance = 1e-3

)Pkg.add(url="https://github.com/MarvynBailly/ManufacturedSolutions.jl")# Step 3: Automatically compute forcing term!

result = run_simple_mms_test(test)  # Pass/fail based on error tolerance

``````forcing = compute_forcing_term_1d(pde_lhs, u_mms_sym, x, u)



**Convergence Tests** - Verify order of accuracy:# Result: f(x) = π²sin(πx) ≈ 9.8696sin(πx)  ✓ Computed automatically!

```julia

test = ConvergenceTest(Or clone and develop locally:

    "Test name",

    solver_function,# Step 4: Verify your solver

    pde,

    manufactured_solution,```bashu_mms_func = build_function(u_mms_sym, x, expression=Val{false})

    (x, u),

    (0.0, 1.0),git clone https://github.com/MarvynBailly/ManufacturedSolutions.jlu_mms = ManufacturedSolution(u_mms_func, (0.0, 1.0))

    expected_order = 2.0,       # e.g., 2nd order method

    base_resolution = 0.1,cd ManufacturedSolutions.jl

    num_refinements = 4          # test at h, h/2, h/4, h/8

)julia --project=. -e 'using Pkg; Pkg.instantiate()'solver_wrapper(grid, mms) = my_poisson_solver(grid, forcing)

result = run_convergence_test(test)  # Measures actual convergence rate

``````solver = UserSolver(solver_wrapper, name="My Solver")



### Automatic Symbolic Computation



The package automatically computes forcing terms using [Symbolics.jl](https://symbolics.juliasymbolics.org/):## What is the Method of Manufactured Solutions?results = verify_solver(solver, pde_lhs, u_mms, 



**Spatial PDEs**:                       resolutions=[0.1, 0.05, 0.025])

```julia

@variables x u(x)The Method of Manufactured Solutions (MMS) is a powerful technique for verifying that your numerical code is implemented correctly:

Dx = Differential(x)

# Step 5: Check convergence

# Any PDE you can write symbolically

pde = -Dx(Dx(u)) + k*u  # Reaction-diffusion**The Problem**: How do you know your PDE solver works if you don't know the exact solution?println("L2 convergence rate: ", results.convergence_rates[:L2])

forcing = compute_forcing_term_1d(pde, sin(π*x), x, u)

```# Output: L2 convergence rate: 2.002  ✓ (2nd order accurate)



**Time-dependent PDEs**:**The Solution**: Work backwards!```

```julia

@variables x t u(x,t)1. Choose a solution you **do** know (the "manufactured solution")

Dx = Differential(x)

Dt = Differential(t)2. Plug it into your PDE to find what forcing term would produce it### C Code Injection Example



# Heat equation3. Solve your PDE with that forcing term

pde = Dt(u) - κ*Dx(Dx(u))

forcing = compute_forcing_term_time_dependent(pde, exp(-t)*sin(π*x), x, t, u)4. Check if you get back the solution you started withVerify a C solver by automatically generating and injecting forcing code:

```



### Test Organization

**Example**: For the Poisson equation `-u'' = f`:```julia

Tests are organized by domain type for comprehensive coverage:

- Choose: `u(x) = sin(πx)` (we know this exactly!)using ManufacturedSolutions

```julia

# Different domain sizes and locations- Compute: `u''(x) = -π² sin(πx)`, so we need `f(x) = π² sin(πx)`

test_suite = [

    ("Square domain [0,1]",       (0.0, 1.0)),- Solve: Run your solver with `f(x) = π² sin(πx)`# Step 1: Define PDE and manufactured solution

    ("Rectangle [0,2]",           (0.0, 2.0)),

    ("Shifted domain [1,2]",      (1.0, 2.0)),- Verify: Does your numerical solution match `sin(πx)`?@variables x u(x)

    ("Long domain [0,10]",        (0.0, 10.0)),

    ("Irregular [-1.5, 3.2]",     (-1.5, 3.2))Dx = Differential(x)

]

If your solver is correct, the numerical and analytical solutions should match (up to discretization error), and the error should decrease at the expected rate as you refine the grid.pde_lhs = -Dx(Dx(u))

# Different manufactured solutions

solutions = [u_mms = sin(π * x)

    sin(π*x),           # Smooth trigonometric

    x*(1-x),            # Polynomial## Features

    exp(-x)*sin(π*x),   # Exponential decay

    sin(3*π*x)          # Higher frequency# Step 2: Compute forcing symbolically

]

```### Verification Methodsforcing = compute_forcing_term_1d(pde_lhs, u_mms, x, u)



## Example Test Results



When you run a convergence test, you get detailed output showing the convergence behavior:**Simple Unit Tests** - Test single cases:# Step 3: Generate C code



``````juliac_code = generate_forcing_code(forcing, :c, [:x])

======================================================================

Convergence Test: Square Domain [0,1]test = SimpleMMSTest(# Generates: double forcing(double x) { return 9.8696*sin(M_PI*x); }

======================================================================

Domain: [0.0, 1.0]    "Test name",

Expected Order: 2.0

Testing 4 refinement levels    solver_function,# Step 4: Inject into your C solver file



[2/3] Running solver on multiple resolutions...    pde,inject_forcing_code("my_solver.c", c_code, "my_solver_mms.c")



  Level 0: h = 0.100000, nx =   11, L2 error = 5.845e-03    manufactured_solution,# ✓ Code injected at // MMS_FORCING_FUNCTION_HERE marker

  Level 1: h = 0.050000, nx =   21, L2 error = 1.456e-03

  Level 2: h = 0.025000, nx =   41, L2 error = 3.636e-04    (x, u),

  Level 3: h = 0.012500, nx =   81, L2 error = 9.088e-05

    (0.0, 1.0),# Step 5: Compile and run your C code

[3/3] Computing observed convergence order...

    resolution = 0.01,# gcc my_solver_mms.c -lm -o solver && ./solver

  Error Reduction Ratios:

    Level 0→1: 4.015 (expected: 4.000 for order 2.0) ✓    tolerance = 1e-3```

    Level 1→2: 4.004 (expected: 4.000 for order 2.0) ✓

    Level 2→3: 4.001 (expected: 4.000 for order 2.0) ✓)



  Observed Order:  2.002result = run_simple_mms_test(test)  # Pass/fail based on error tolerance**Supported languages**: Julia, C, Python, Fortran

  Expected Order:  2.000

  Difference:      0.002 < tolerance (0.2)```



  ✓ CONVERGENCE TEST PASSED---

```

**Convergence Tests** - Verify order of accuracy:

The key insight: For a 2nd order method, when you halve the grid spacing (h → h/2), the error should decrease by a factor of 4 (since 2² = 4). The test verifies this automatically.

```julia## Features

## Running the Tests

test = ConvergenceTest(

The package includes comprehensive test suites demonstrating correct and broken implementations:

    "Test name",### ✅ Currently Available

### Convergence Tests

```bash    solver_function,

# Verify Poisson solver across 8 different domains

julia --project=. test/convergence_poisson_1d.jl    pde,#### Core Verification



# Compare correct vs. broken reaction-diffusion solver    manufactured_solution,- **Symbolic Forcing Computation**: Automatic derivation using Symbolics.jl - no manual calculus!

julia --project=. test/convergence_reaction_diffusion.jl

```    (x, u),- **Multi-Language Support**: Generate code for C, Fortran, Python, Julia



### Unit Tests    (0.0, 1.0),- **Code Injection**: Automatically inject forcing functions into source files

```bash

# 9 unit tests with different scenarios    expected_order = 2.0,       # e.g., 2nd order method- **1D Problem Support**: Full support for 1D PDEs (elliptic, parabolic, hyperbolic)

julia --project=. test/verification_1d_poisson.jl

    base_resolution = 0.1,- **Error Norms**: L1, L2, L∞ with proper discrete formulations

# Demonstrate error detection

julia --project=. test/verification_reaction_diffusion.jl    num_refinements = 4          # test at h, h/2, h/4, h/8- **Convergence Analysis**: Automatic order of accuracy estimation via log-log regression



# Side-by-side comparison)- **Time-Dependent PDEs**: Support for parabolic/hyperbolic problems

julia --project=. test/compare_solvers.jl

```result = run_convergence_test(test)  # Measures actual convergence rate



All tests pass for the correct solvers and correctly identify bugs in the broken solvers!```#### User Experience



## Solvers Included- **Clean API**: Simple, intuitive function interface



The `test/solvers/` directory contains reference implementations:### Automatic Symbolic Computation- **Comprehensive Results**: Timing, error data, convergence rates, R² goodness-of-fit



**Poisson Equation** (`-u'' = f`):- **Flexible Integration**: Works with your existing solvers

- `poisson_solver_dirichlet`: Homogeneous Dirichlet BCs ✓

- `poisson_solver_periodic`: Periodic BCs ✓The package automatically computes forcing terms using [Symbolics.jl](https://symbolics.juliasymbolics.org/):- **Non-Intrusive**: Original files unchanged (creates new files)

- `poisson_solver_neumann`: Neumann BCs ✓

- **Detailed Logging**: Optional verbose output for debugging

**Reaction-Diffusion** (`-u'' + k*u = f`):

- `reaction_diffusion_solver_correct`: Proper discretization ✓**Spatial PDEs**:

- `reaction_diffusion_solver_broken`: Missing reaction term (for demo)

- `reaction_diffusion_solver_very_broken`: Wrong sign (for demo)```julia### 🔄 Coming Soon



These serve as both reference implementations and test cases.@variables x u(x)



## Key ComponentsDx = Differential(x)- 2D and 3D grid support



### Core Types- Curvilinear/elliptic grid support (for grid generators)



- `Grid1D(a, b, h)` - 1D uniform grid from a to b with spacing h# Any PDE you can write symbolically- More boundary condition types

- `ForcingTerm` - Callable forcing function with symbolic form

- `SimpleMMSTest` - Single resolution test with error tolerancepde = -Dx(Dx(u)) + k*u  # Reaction-diffusion- Power operations in code generation (x^2 → pow(x,2))

- `ConvergenceTest` - Multi-resolution test with order verification

forcing = compute_forcing_term_1d(pde, sin(π*x), x, u)- Automatic compilation and execution

### Main Functions

```- Visualization tools

- `compute_forcing_term_1d(pde, u_mms, x, u)` - Compute forcing for spatial PDEs

- `compute_forcing_term_time_dependent(pde, u_mms, x, t, u)` - For time-dependent PDEs- Symbolic BC extraction

- `run_simple_mms_test(test)` - Run single-resolution verification

- `run_convergence_test(test)` - Run multi-resolution convergence study**Time-dependent PDEs**:



### Error Norms```julia---



- `compute_L1_error(u_num, u_exact, grid)` - Discrete L¹ norm@variables x t u(x,t)

- `compute_L2_error(u_num, u_exact, grid)` - Discrete L² norm

- `compute_Linf_error(u_num, u_exact)` - Maximum norm (L^∞)Dx = Differential(x)## Installation



## How to Use With Your SolverDt = Differential(t)



### 1. Write Your Solver with the Right Signature```julia



Your solver should accept a grid and forcing function:# Heat equation# Clone the repository



```juliapde = Dt(u) - κ*Dx(Dx(u))git clone https://github.com/MarvynBailly/ManufacturedSolutions.jl

function my_solver(grid::Grid1D, forcing::ForcingTerm)

    # Extract grid informationforcing = compute_forcing_term_time_dependent(pde, exp(-t)*sin(π*x), x, t, u)

    x = grid.x        # Grid points

    nx = grid.nx      # Number of points```# Activate and instantiate

    h = grid.dx       # Grid spacing

    using Pkg

    # Build your system of equations

    A = zeros(nx, nx)### Test OrganizationPkg.activate("path/to/ManufacturedSolutions.jl")

    b = zeros(nx)

    Pkg.instantiate()

    for i in 1:nx

        # Evaluate forcing at grid pointTests are organized by domain type for comprehensive coverage:```

        f_i = forcing(x[i])

        

        # Your discretization here...

        # (set up your matrix and RHS)```julia### Dependencies

    end

    # Different domain sizes and locations

    # Solve and return numerical solution

    u_numerical = A \ btest_suite = [- Julia 1.0+

    return u_numerical

end    ("Square domain [0,1]",       (0.0, 1.0)),- Symbolics.jl (for automatic differentiation)

```

    ("Rectangle [0,2]",           (0.0, 2.0)),- LinearAlgebra (stdlib)

### 2. Create a Test

    ("Shifted domain [1,2]",      (1.0, 2.0)),- Printf (stdlib)

```julia

using ManufacturedSolutions    ("Long domain [0,10]",        (0.0, 10.0)),- Statistics (stdlib)



# Define the PDE your solver solves    ("Irregular [-1.5, 3.2]",     (-1.5, 3.2))

@variables x u(x)

Dx = Differential(x)]---

pde = -Dx(Dx(u))  # Your PDE



# Pick a manufactured solution

u_mms = sin(π * x)# Different manufactured solutions## Complete Usage Guide



# Create convergence testsolutions = [

test = ConvergenceTest(

    "My solver verification",    sin(π*x),           # Smooth trigonometric### Method of Manufactured Solutions

    my_solver,

    pde,    x*(1-x),            # Polynomial

    u_mms,

    (x, u),    exp(-x)*sin(π*x),   # Exponential decayMMS is a powerful code verification technique that works by:

    (0.0, 1.0),

    expected_order = 2.0,    sin(3*π*x)          # Higher frequency

    base_resolution = 0.1,

    num_refinements = 4]1. **Choose** a manufactured solution `u(x)` that you know analytically

)

```2. **Substitute** it into your PDE to compute the required forcing term `f(x)`

# Run it!

result = run_convergence_test(test)3. **Solve** your PDE with this forcing term using your numerical method

```

## Example Test Results4. **Compare** the numerical solution to the known manufactured solution

### 3. Interpret Results

5. **Analyze** convergence rates to verify your solver's order of accuracy

- **Test passes**: Your solver is correctly implemented! ✓

- **Test fails (wrong order)**: Discretization error - check your stencilsWhen you run a convergence test, you get detailed output showing the convergence behavior:

- **Test fails (large errors)**: Implementation bug - check your PDE

#### Example: Poisson Equation

## Project Structure

```

```

ManufacturedSolutions.jl/======================================================================**Given PDE**: `-u''(x) = f(x)` on [0,1]

├── src/

│   ├── ManufacturedSolutions.jl    # Main moduleConvergence Test: Square Domain [0,1]

│   ├── core/

│   │   └── types.jl                # Core type definitions======================================================================**Choose**: `u(x) = sin(πx)`

│   ├── symbolic/

│   │   └── forcing_terms.jl        # Automatic forcing computationDomain: [0.0, 1.0]

│   └── verification/

│       ├── error_norms.jl          # L1, L2, Linf normsExpected Order: 2.0**Compute forcing**:

│       ├── simple_unit_test.jl     # Single-resolution tests

│       └── convergence_test.jl     # Multi-resolution testsTesting 4 refinement levels- u'(x) = π cos(πx)

├── test/

│   ├── solvers/                    # Reference solver implementations- u''(x) = -π² sin(πx)

│   ├── convergence_poisson_1d.jl   # Comprehensive Poisson tests

│   └── ...[2/3] Running solver on multiple resolutions...- **Forcing**: f(x) = -(-π² sin(πx)) = π² sin(πx)

└── README.md                       # This file

```



## Requirements  Level 0: h = 0.100000, nx =   11, L2 error = 5.845e-03**Verify**: Your solver should reproduce u(x) = sin(πx) when given f(x) = π² sin(πx).



- Julia 1.0 or later  Level 1: h = 0.050000, nx =   21, L2 error = 1.456e-03

- Symbolics.jl - for symbolic differentiation

- LinearAlgebra - for error norms (standard library)  Level 2: h = 0.025000, nx =   41, L2 error = 3.636e-04**ManufacturedSolutions.jl automates steps 2-5!**

- Test.jl - for running tests (standard library)

  Level 3: h = 0.012500, nx =   81, L2 error = 9.088e-05

## Current Limitations

---

- **1D only**: Currently supports 1D PDEs (2D/3D coming soon)

- **Uniform grids**: Requires uniformly spaced grids[3/3] Computing observed convergence order...

- **Dirichlet BCs**: Best support for homogeneous Dirichlet boundary conditions

- **Julia solvers**: Direct support for Julia solvers### Symbolic Interface



## Contributing  Error Reduction Ratios:



This is an active research project! Contributions welcome:    Level 0→1: 4.015 (expected: 4.000 for order 2.0) ✓The symbolic interface automatically computes forcing terms - no manual calculus needed!

- Bug reports and fixes

- New test cases    Level 1→2: 4.004 (expected: 4.000 for order 2.0) ✓

- Documentation improvements

- Feature requests    Level 2→3: 4.001 (expected: 4.000 for order 2.0) ✓#### 1D Spatial Problems (Elliptic PDEs)



## License



MIT License - see LICENSE file  Observed Order:  2.002```julia



## Citation  Expected Order:  2.000using ManufacturedSolutions



If you use this package in your research, please cite:  Difference:      0.002 < tolerance (0.2)



```bibtex@variables x u(x)

@software{manufactured_solutions_jl,

  author = {Bailly, Marvyn},  ✓ CONVERGENCE TEST PASSEDDx = Differential(x)

  title = {ManufacturedSolutions.jl: Automatic Verification of Numerical PDE Solvers},

  year = {2024},```

  url = {https://github.com/MarvynBailly/ManufacturedSolutions.jl}

}# Define your PDE

```

The key insight: For a 2nd order method, when you halve the grid spacing (h → h/2), the error should decrease by a factor of 4 (since 2² = 4). The test verifies this automatically.pde = -Dx(Dx(u))  # Poisson equation

## References



- **Method of Manufactured Solutions**: Roache, P.J. (2002). "Code Verification by the Method of Manufactured Solutions." *Journal of Fluids Engineering*, 124(1), 4-10.

- **Symbolics.jl**: Gowda, S., et al. (2021). "High-performance symbolic-numerics via multiple dispatch." arXiv:2105.03949## Running the Tests# Choose manufactured solution



## Acknowledgmentsu_mms = sin(π*x)



Developed for verification of numerical PDE solvers and grid generation codes.The package includes comprehensive test suites demonstrating correct and broken implementations:



---# Compute forcing automatically



**Status**: Active Development  ### Convergence Testsforcing = compute_forcing_term_1d(pde, u_mms, x, u)

**Author**: Marvyn Bailly  

**Version**: 0.1.0```bash


# Verify Poisson solver across 8 different domains# Use the forcing:

julia --project=. test/convergence_poisson_1d.jlprintln(forcing.symbolic_form)  # Display: 9.8696sin(πx)

f_value = forcing(0.5)          # Evaluate at x=0.5

# Compare correct vs. broken reaction-diffusion solver```

julia --project=. test/convergence_reaction_diffusion.jl

```#### Time-Dependent Problems (Parabolic/Hyperbolic PDEs)



### Unit Tests```julia

```bash@variables x t u(x,t)

# 9 unit tests with different scenariosDx = Differential(x)

julia --project=. test/verification_1d_poisson.jlDt = Differential(t)



# Demonstrate error detection# Heat equation: u_t - κu_xx = f

julia --project=. test/verification_reaction_diffusion.jlκ = 0.1

pde = Dt(u) - κ*Dx(Dx(u))

# Side-by-side comparison

julia --project=. test/compare_solvers.jl# Manufactured solution

```u_mms = exp(-t) * sin(2π*x)



All tests pass for the correct solvers and correctly identify bugs in the broken solvers!# Compute forcing automatically

forcing = compute_forcing_term_time_dependent(pde, u_mms, x, t, u)

## Solvers Included

# Use in solver:

The `test/solvers/` directory contains reference implementations:f_value = forcing(0.5, 0.1)  # Evaluate at x=0.5, t=0.1

```

**Poisson Equation** (`-u'' = f`):

- `poisson_solver_dirichlet`: Homogeneous Dirichlet BCs (working ✓)#### Benefits

- `poisson_solver_periodic`: Periodic BCs (has known issues)

- `poisson_solver_neumann`: Neumann BCs (under development)✅ **No manual derivation** - Symbolics.jl does the calculus  

✅ **Fewer errors** - No mistakes in derivative computation  

**Reaction-Diffusion** (`-u'' + k*u = f`):✅ **Complex solutions** - Try solutions you wouldn't want to differentiate by hand  

- `reaction_diffusion_solver_correct`: Proper discretization ✓✅ **Rapid iteration** - Change manufactured solution instantly  

- `reaction_diffusion_solver_broken`: Missing reaction term (for demo)✅ **Self-documenting** - Symbolic form shows exactly what forcing is used

- `reaction_diffusion_solver_very_broken`: Wrong sign (for demo)

---

These serve as both reference implementations and test cases.

### Code Injection Workflow

## Key Components

Verify solvers written in **any programming language** by generating and injecting forcing code.

### Core Types

#### Step 1: Prepare Your Solver

- `Grid1D(a, b, h)` - 1D uniform grid from a to b with spacing h

- `ForcingTerm` - Callable forcing function with symbolic formAdd a marker comment where the forcing function should go:

- `SimpleMMSTest` - Single resolution test with error tolerance

- `ConvergenceTest` - Multi-resolution test with order verification**C Example** (`my_poisson.c`):

```c

### Main Functions#include <math.h>



- `compute_forcing_term_1d(pde, u_mms, x, u)` - Compute forcing for spatial PDEs// MMS_FORCING_FUNCTION_HERE

- `compute_forcing_term_time_dependent(pde, u_mms, x, t, u)` - For time-dependent PDEs// Default implementation (will be replaced)

- `run_simple_mms_test(test)` - Run single-resolution verificationdouble forcing(double x) {

- `run_convergence_test(test)` - Run multi-resolution convergence study    return 0.0;

}

### Error Norms

void solve_poisson(int n, double *x, double *u) {

- `compute_L1_error(u_num, u_exact, grid)` - Discrete L¹ norm    // Your solver implementation...

- `compute_L2_error(u_num, u_exact, grid)` - Discrete L² norm    for (int i = 0; i < n; i++) {

- `compute_Linf_error(u_num, u_exact)` - Maximum norm (L^∞)        double f = forcing(x[i]);

        // ... use f in your discretization ...

## How to Use With Your Solver    }

}

### 1. Write Your Solver with the Right Signature```



Your solver should accept a grid and forcing function:**Python Example** (`my_solver.py`):

```python

```juliaimport numpy as np

function my_solver(grid::Grid1D, forcing::ForcingTerm)

    # Extract grid information# MMS_FORCING_FUNCTION_HERE

    x = grid.x        # Grid pointsdef forcing(x):

    nx = grid.nx      # Number of points    return 0.0

    h = grid.dx       # Grid spacing

    def solve(x):

    # Build your system of equations    f = forcing(x)

    A = zeros(nx, nx)    # ... your solver ...

    b = zeros(nx)```

    

    for i in 1:nx**Fortran Example** (`solver.f90`):

        # Evaluate forcing at grid point```fortran

        f_i = forcing(x[i])! MMS_FORCING_FUNCTION_HERE

        function forcing(x) result(f)

        # Your discretization here...    real(8) :: x, f

        # (set up your matrix and RHS)    f = 0.0d0

    endend function

    ```

    # Solve and return numerical solution

    u_numerical = A \ b#### Step 2: Generate and Inject Code

    return u_numerical

end```julia

```using ManufacturedSolutions



### 2. Create a Test# Define problem

@variables x u(x)

```juliaDx = Differential(x)

using ManufacturedSolutionspde = -Dx(Dx(u))

u_mms = sin(π*x)

# Define the PDE your solver solves

@variables x u(x)# Compute forcing

Dx = Differential(x)forcing = compute_forcing_term_1d(pde, u_mms, x, u)

pde = -Dx(Dx(u))  # Your PDE

# Generate C code

# Pick a manufactured solutionc_code = generate_forcing_code(forcing, :c, [:x])

u_mms = sin(π * x)# Produces:

# double forcing(double x) {

# Create convergence test#     return 9.869604401089358*sin(M_PI*x);

test = ConvergenceTest(# }

    "My solver verification",

    my_solver,# Inject into your C file

    pde,inject_forcing_code(

    u_mms,    "my_poisson.c",       # Original file

    (x, u),    c_code,               # Generated code

    (0.0, 1.0),    "my_poisson_mms.c"    # Output file

    expected_order = 2.0,)

    base_resolution = 0.1,# ✓ Code injected successfully!

    num_refinements = 4```

)

#### Step 3: Compile and Run

# Run it!

result = run_convergence_test(test)```bash

```# Compile the modified solver

gcc my_poisson_mms.c -lm -o solver

### 3. Interpret Results

# Run at different resolutions

- **Test passes**: Your solver is correctly implemented! ✓./solver  # Outputs solution to file

- **Test fails (wrong order)**: Discretization error - check your stencils```

- **Test fails (large errors)**: Implementation bug - check your PDE

#### Step 4: Analyze Results

## Project Structure

```julia

```using DelimitedFiles

ManufacturedSolutions.jl/

├── src/# Read numerical solution

│   ├── ManufacturedSolutions.jl    # Main moduledata = readdlm("solution.txt", skipstart=1)

│   ├── core/x_num = data[:, 1]

│   │   └── types.jl                # Core type definitionsu_num = data[:, 2]

│   ├── symbolic/

│   │   └── forcing_terms.jl        # Automatic forcing computation# Compute analytical solution

│   └── verification/u_exact = sin.(π .* x_num)

│       ├── error_norms.jl          # L1, L2, Linf norms

│       ├── simple_unit_test.jl     # Single-resolution tests# Compute error

│       └── convergence_test.jl     # Multi-resolution testsh = x_num[2] - x_num[1]

├── test/L2_error = sqrt(h * sum((u_num .- u_exact).^2))

│   ├── solvers/                    # Reference solver implementationsprintln("L2 error: ", L2_error)

│   ├── convergence_poisson_1d.jl   # Comprehensive Poisson tests

│   └── ...# Repeat for multiple resolutions and check convergence rate

└── README.md                       # This file```

```

#### Supported Languages

## Requirements

| Language | Template | Math Functions | Constants |

- Julia 1.0 or later|----------|----------|----------------|-----------|

- Symbolics.jl - for symbolic differentiation| **C** | `double forcing(double x)` | `sin`, `cos`, `exp` | `M_PI` |

- LinearAlgebra - for error norms (standard library)| **Python** | `def forcing(x):` | `math.sin`, `math.cos` | `math.pi` |

- Test.jl - for running tests (standard library)| **Fortran** | `function forcing(x)` | `dsin`, `dcos`, `dexp` | `4.0d0*datan(1.0d0)` |

| **Julia** | `forcing(x)` | `sin`, `cos`, `exp` | `π` |

## Current Limitations

---

- **1D only**: Currently supports 1D PDEs (2D/3D coming soon)

- **Uniform grids**: Requires uniformly spaced grids## API Reference

- **Dirichlet BCs**: Best support for homogeneous Dirichlet boundary conditions

- **Julia solvers**: Direct support for Julia solvers (multi-language support in development)### Core Types



## Contributing#### `Grid1D`

```julia

This is an active research project! Contributions welcome:Grid1D(left, right, nx)   # Create with number of points

- Bug reports and fixesGrid1D(left, right, dx)   # Create with spacing

- New test cases```

- Documentation improvementsUniform 1D grid with specified domain and resolution.

- Feature requests

#### `ManufacturedSolution`

## License```julia

ManufacturedSolution(expression, domain; name="", description="")

MIT License - see LICENSE file```

Wrapper for analytical solution function.

## Citation

#### `ForcingTerm`

If you use this package in your research, please cite:```julia

forcing.symbolic_form  # Symbolic expression

```bibtexforcing(x)            # Evaluate at point

@software{manufactured_solutions_jl,forcing(x, t)         # Time-dependent evaluation

  author = {Bailly, Marvyn},```

  title = {ManufacturedSolutions.jl: Automatic Verification of Numerical PDE Solvers},Callable forcing function with symbolic representation.

  year = {2024},

  url = {https://github.com/MarvynBailly/ManufacturedSolutions.jl}#### `UserSolver`

}```julia

```UserSolver(solve_function; name="", description="")

```

## ReferencesWrapper for your solver function. Function signature:

```julia

- **Method of Manufactured Solutions**: Roache, P.J. (2002). "Code Verification by the Method of Manufactured Solutions." *Journal of Fluids Engineering*, 124(1), 4-10.function my_solver(grid::Grid1D, mms::ManufacturedSolution; kwargs...)

- **Symbolics.jl**: Gowda, S., et al. (2021). "High-performance symbolic-numerics via multiple dispatch." arXiv:2105.03949    # ... your implementation ...

    return numerical_solution::Vector{Float64}

## Acknowledgmentsend

```

Developed for verification of numerical PDE solvers and grid generation codes.

#### `VerificationResult`

---```julia

results.convergence_data   # Per-resolution error data

**Status**: Active Development  results.convergence_rates  # Dict{Symbol,Float64} with rates

**Author**: Marvyn Bailly  results.problem_info       # Problem metadata

**Version**: 0.4.0results.solver_info        # Solver metadata

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
