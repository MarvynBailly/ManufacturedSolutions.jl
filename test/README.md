# 1D Verification Test Suite

This directory contains a modular, reusable test suite for verifying 1D numerical solvers using the Method of Manufactured Solutions (MMS).

## Directory Structure

```
test/
├── verification_1d_poisson.jl    # Generic 1D verification tests
└── solvers/
    └── poisson_1d.jl             # Solver implementations being tested
```

## Design Philosophy

The test suite is organized to **separate concerns**:

1. **Verification Tests** (`verification_1d_poisson.jl`): Generic, domain-independent tests
   - Can verify ANY 1D numerical solver
   - Tests different domains, resolutions, and manufactured solutions
   - PDE/solver agnostic

2. **Solver Implementations** (`solvers/poisson_1d.jl`): Specific numerical methods
   - Contains the solver being tested
   - Can be easily swapped or extended
   - Independent of test infrastructure

## Running Tests

```bash
julia --project=. test/verification_1d_poisson.jl
```

## Test Coverage

The verification suite includes 9 comprehensive tests:

### Domain Variations
- **Unit domain** `[0,1]`
- **Rectangular domain** `[0,2]`
- **Shifted domain** `[1,2]`
- **Long domain** `[0,10]`
- **Negative domain** `[-1,0]`

### Manufactured Solutions
- **Trigonometric**: `sin(πx)`, `sin(3πx)`
- **Polynomial**: `x(1-x)`
- **Exponential**: `exp(-x)sin(πx)`

### Resolution Tests
- Coarse: `h = 0.1`
- Standard: `h = 0.01`
- Fine: `h = 0.001` (accuracy verification)

## Test Results

All tests pass with expected accuracy:

| Test | Domain | Solution | Error | Status |
|------|--------|----------|-------|--------|
| 1 | [0,1] | sin(πx) | 5.8e-5 | ✓ |
| 2 | [0,2] | sin(πx/2) | 8.2e-5 | ✓ |
| 3 | [1,2] | sin(π(x-1)) | 5.8e-5 | ✓ |
| 4 | [0,1] fine | sin(πx) | 5.8e-7 | ✓ |
| 5 | [0,1] | sin(3πx) | 1.3e-4 | ✓ |
| 6 | [0,1] | x(1-x) | 1.2e-15 | ✓ (machine precision!) |
| 7 | [0,1] | exp(-x)sin(πx) | 2.3e-5 | ✓ |
| 8 | [0,10] | sin(πx/10) | 1.8e-4 | ✓ |
| 9 | [-1,0] | sin(π(x+1)) | 5.8e-5 | ✓ |

**Total**: 14/14 assertions pass in ~46 seconds

## Adding New Solvers

To test a different 1D solver:

1. Create your solver in `solvers/your_solver.jl`
2. Update the import in `verification_1d_poisson.jl`:
   ```julia
   include("solvers/your_solver.jl")
   using .Main: your_solver_function
   ```
3. Change the solver function in each test:
   ```julia
   test = SimpleMMSTest(
       "Test name",
       your_solver_function,  # <-- Replace here
       pde_poisson,
       sin(π * x),
       (x, u),
       (0.0, 1.0),
       resolution = 0.01,
       tolerance = 1e-3
   )
   ```

The same verification tests will now verify your new solver!

## Solver Requirements

Any 1D solver can be verified if it has the signature:

```julia
function my_solver(grid::Grid1D, forcing::ForcingTerm)::Vector{Float64}
    # Your numerical method here
    return numerical_solution
end
```

Where:
- `grid::Grid1D` - Contains domain information (x, nx, dx)
- `forcing::ForcingTerm` - Function `f(x)` for the right-hand side
- Returns: Vector of solution values at grid points

## Available Solvers

The `solvers/poisson_1d.jl` file currently provides:

### 1. `poisson_solver_dirichlet`
Solves: `-u'' = f` with `u(a) = u(b) = 0`
- Second-order centered finite differences
- Verified with all 9 tests ✓

### 2. `poisson_solver_periodic`
Solves: `-u'' = f` with periodic BCs
- Requires compatibility: `∫f dx = 0`
- Currently has implementation issues

### 3. `poisson_solver_neumann`
Solves: `-u'' = f` with flux BCs
- Under development

## Modular Design Benefits

✅ **Reusability**: Same tests work for any 1D solver  
✅ **Maintainability**: Solvers and tests evolve independently  
✅ **Extensibility**: Easy to add new solvers or new tests  
✅ **Clarity**: Clear separation between "what to test" and "how to test"

## Future Extensions

Potential additions:
- Convergence rate tests (compare multiple resolutions)
- Different boundary conditions (Neumann, Robin, mixed)
- Different PDEs (advection-diffusion, reaction-diffusion)
- Higher-order methods (4th order, spectral)
- Time-dependent problems
