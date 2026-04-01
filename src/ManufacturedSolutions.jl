module ManufacturedSolutions

using LinearAlgebra
using Symbolics
using Printf

# Re-export useful Symbolics functionality
using Symbolics: @variables, Differential, expand_derivatives, substitute, simplify, build_function

export @variables, Differential, expand_derivatives, substitute, simplify, build_function

# Core types
include("core/types.jl")
export AbstractManufacturedSolution, AbstractGrid, AbstractBoundaryCondition
export Grid1D, ManufacturedSolution, ForcingTerm
export VerificationProblem, UserSolver
export ConvergenceData, VerificationResult

# Symbolic processing
include("symbolic/forcing_terms.jl")
export compute_forcing_term, compute_forcing_term_1d, compute_forcing_term_time_dependent

# Code generation and injection
include("codegen/code_injection.jl")
export CodeTemplate, generate_forcing_code, inject_forcing_code
export create_solver_template, symbolic_to_code

# Verification functionality
include("verification/error_norms.jl")
include("verification/convergence.jl")
include("verification/verify.jl")
include("verification/code_injection_workflow.jl")
include("verification/unit_test.jl")
include("verification/simple_unit_test.jl")
include("verification/convergence_test.jl")

# Main interface
export verify_solver
export compute_L1_error, compute_L2_error, compute_Linf_error, compute_all_errors
export estimate_convergence_rate, analyze_convergence

# Code injection workflow
export MMSProblem, SolverSetup, verify_with_code_injection

# Unit testing interface
export MMSUnitTest, run_mms_unit_test, @mms_test

# Simple unit testing (function wrapper approach)
export SimpleMMSTest, run_simple_mms_test, @test_solver

# Convergence-based verification (order of accuracy)
export ConvergenceTest, run_convergence_test

end