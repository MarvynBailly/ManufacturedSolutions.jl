# MMS Verification Demo: Correct vs. Broken Solvers

## 🎯 What We Demonstrated

We created **three different implementations** of a reaction-diffusion ODE solver and used MMS to verify them:

1. ✅ **Correct Implementation** - Proper discretization
2. ❌ **Broken Implementation** - Missing reaction term
3. ❌❌ **Very Broken Implementation** - Wrong sign in discretization

## 📊 Results Summary

### Error Comparison Table

| Solver | L2 Error | Error Ratio | Status |
|--------|----------|-------------|--------|
| ✅ Correct | **5.28e-05** | 1x (baseline) | ✅ PASS |
| ❌ Broken | **7.17e-02** | **1358x larger** | ❌ FAIL |
| ❌❌ Very Broken | **1.60e-01** | **3021x larger** | ❌ FAIL |

### Visual Error Comparison

```
Correct:       █ (5.3e-5)
               ↓
Broken:        ████████████████████ (~1400x larger, 7.2e-2)
               ↓
Very Broken:   ████████████████████████████████████ (~3000x larger, 1.6e-1)
```

## 🔍 What MMS Caught

### Bug #1: Missing Reaction Term
**Problem**: Solver forgot to include `k*u` term in discretization
```julia
# Should be:
A[i, i] = 2.0 + k * h^2

# But was:
A[i, i] = 2.0  # ❌ Missing k*h^2
```

**Detection**: MMS found errors **1358x larger** than correct solver!

### Bug #2: Wrong Sign
**Problem**: Solver used wrong sign for reaction term
```julia
# Should be:
A[i, i] = 2.0 + k * h^2

# But was:
A[i, i] = 2.0 - k * h^2  # ❌ Wrong sign
```

**Detection**: MMS found errors **3021x larger** than correct solver!

## 📈 Detailed Test Results

### Test 1: ✅ Correct Solver
```
L1 error:   4.75e-05  ✓
L2 error:   5.28e-05  ✓  (baseline)
L∞ error:   7.47e-05  ✓
Status:     PASS
```

### Test 2: ❌ Broken Solver (Missing Term)
```
L1 error:   6.46e-02  ✗  (1358x worse)
L2 error:   7.17e-02  ✗  (1358x worse)
L∞ error:   1.01e-01  ✗  (1358x worse)
Status:     FAIL
```

### Test 3: ❌❌ Very Broken Solver (Wrong Sign)
```
L1 error:   1.44e-01  ✗  (3021x worse)
L2 error:   1.60e-01  ✗  (3021x worse)
L∞ error:   2.26e-01  ✗  (3021x worse)
Status:     FAIL
```

## 🎓 Key Lessons

### 1. MMS Works as an Automatic Bug Detector
- Correct implementation: All tests pass ✅
- Buggy implementations: All tests fail ❌
- No manual inspection needed!

### 2. Quantitative Error Analysis
MMS doesn't just say "pass" or "fail" - it tells you **how wrong** you are:
- Small errors (< 1e-3): Likely correct ✅
- Medium errors (~0.01-0.1): Missing terms ❌
- Large errors (> 0.1): Major bugs ❌❌

### 3. Multiple Tests Increase Confidence
We tested with different manufactured solutions:
- `sin(πx)` - Smooth trigonometric
- `x(1-x)` - Polynomial (correct solver achieves machine precision!)
- `exp(-x)sin(πx)` - Mixed exponential/trig

All three caught the same bugs consistently.

### 4. Easy to Add New Solvers
To test your own solver:
```julia
# Just implement this interface:
function my_solver(grid::Grid1D, forcing::ForcingTerm)
    # Your code here
    return solution
end

# Run verification:
test = SimpleMMSTest("My test", my_solver, pde, ...)
result = run_simple_mms_test(test)
```

## 🚀 How to Run These Tests

### Full Verification Suite
```bash
julia --project=. test/verification_reaction_diffusion.jl
```
Runs 6 tests (3 on correct solver, 3 on broken solvers)

### Quick Comparison
```bash
julia --project=. test/compare_solvers.jl
```
Side-by-side comparison with summary table

## 📁 Files Created

```
test/
├── solvers/
│   └── reaction_diffusion_1d.jl        # 3 solver implementations
├── verification_reaction_diffusion.jl   # Full test suite (6 tests)
├── compare_solvers.jl                   # Quick comparison script
├── CORRECT_VS_BROKEN_EXAMPLE.md         # Detailed documentation
└── MMS_DEMO_SUMMARY.md                  # This file
```

## ✨ Conclusion

**MMS verification successfully detected implementation bugs!**

- ✅ Correct solver: Errors ~ 10⁻⁵ (second-order accurate)
- ❌ Broken solvers: Errors 1000-3000x larger
- 🎯 Clear, automatic, quantitative verification

This demonstrates that MMS is a powerful tool for:
1. Verifying numerical implementations
2. Catching subtle discretization errors
3. Providing confidence in solver correctness

**The same framework can now verify ANY 1D solver** - just swap the solver function and manufactured solution!

---

**Next Steps**: Apply this same verification framework to:
- Different PDEs (advection, wave equation, etc.)
- Different boundary conditions (Neumann, Robin, periodic)
- Higher dimensions (2D, 3D)
- Time-dependent problems
