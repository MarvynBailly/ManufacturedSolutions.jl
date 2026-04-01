"""
Core type definitions for ManufacturedSolutions.jl

This file defines the fundamental types used throughout the package:
- Grids: Spatial discretizations (currently 1D uniform grids)
- Manufactured Solutions: Analytical solutions used for verification
- Forcing Terms: Right-hand side terms computed from manufactured solutions
- Results: Data structures for verification results

These types provide a clean interface between the verification framework
and user-provided solvers.
"""

# ===================================================================
# Abstract Types
# ===================================================================

"""
    AbstractManufacturedSolution

Abstract base type for manufactured solutions.
"""
abstract type AbstractManufacturedSolution end

"""
    AbstractGrid

Abstract base type for spatial grids/meshes.
"""
abstract type AbstractGrid end

"""
    AbstractBoundaryCondition

Abstract base type for boundary conditions.
"""
abstract type AbstractBoundaryCondition end

# ===================================================================
# Grid Types
# ===================================================================

"""
    Grid1D{T<:Real}

A uniformly-spaced one-dimensional grid for spatial discretization.

This is the fundamental spatial discretization used by solvers. It represents
evenly-spaced points on an interval [a, b].

# Fields
- `x::Vector{T}`: Coordinates of grid points
- `nx::Int`: Number of grid points
- `dx::T`: Spacing between adjacent points (uniform)
- `domain::Tuple{T,T}`: Domain boundaries (left, right)

# Constructors
```julia
# Create grid with specified number of points
grid = Grid1D(0.0, 1.0, 101)  # 101 points from 0 to 1

# Create grid with specified spacing
grid = Grid1D(0.0, 1.0, 0.01)  # spacing h = 0.01
```

# Example
```julia
grid = Grid1D(0.0, 1.0, 0.1)
# grid.x  = [0.0, 0.1, 0.2, ..., 1.0]
# grid.nx = 11
# grid.dx = 0.1
```
"""
struct Grid1D{T<:Real} <: AbstractGrid
    x::Vector{T}
    nx::Int
    dx::T
    domain::Tuple{T,T}
    
    function Grid1D(x_left::T, x_right::T, nx::Int) where {T<:Real}
        # Input validation
        @assert nx > 1 "Grid must have at least 2 points"
        @assert x_right > x_left "Right boundary must be greater than left boundary"
        
        # Create uniformly spaced points
        x = range(x_left, x_right, length=nx) |> collect
        dx = (x_right - x_left) / (nx - 1)
        
        new{T}(x, nx, dx, (x_left, x_right))
    end
end

# Alternate constructor: specify spacing instead of number of points
# Example: Grid1D(0.0, 1.0, 0.01) creates grid with h=0.01
function Grid1D(x_left::Real, x_right::Real, dx::Real)
    nx = round(Int, (x_right - x_left) / dx) + 1
    Grid1D(promote(x_left, x_right)..., nx)
end

# ===================================================================
# Manufactured Solution Types
# ===================================================================

"""
    ManufacturedSolution{F,D}

Represents a manufactured solution for verification.

# Fields
- `expression::F`: Callable function representing the solution
- `domain::D`: Spatial domain
- `metadata::Dict{Symbol,Any}`: Additional information
"""
struct ManufacturedSolution{F,D} <: AbstractManufacturedSolution
    expression::F
    domain::D
    metadata::Dict{Symbol,Any}
    
    function ManufacturedSolution(expr::F, domain::D; 
                                   name::String="",
                                   description::String="") where {F,D}
        metadata = Dict{Symbol,Any}(
            :name => name,
            :description => description
        )
        new{F,D}(expr, domain, metadata)
    end
end

# Evaluate the manufactured solution
(mms::ManufacturedSolution)(args...) = mms.expression(args...)

# ===================================================================
# Forcing Term Types
# ===================================================================

"""
    ForcingTerm{F,S}

Represents a forcing term computed from a manufactured solution.

# Fields
- `function_form::F`: Callable forcing function
- `symbolic_form::S`: Symbolic representation (if available)
"""
struct ForcingTerm{F,S}
    function_form::F
    symbolic_form::S
end

# Evaluate the forcing term
(forcing::ForcingTerm)(args...) = forcing.function_form(args...)

# ===================================================================
# Verification Problem Types
# ===================================================================

"""
    VerificationProblem{P,M,G}

Represents a complete verification problem.

# Fields
- `pde::P`: PDE expression
- `manufactured_solution::M`: Manufactured solution
- `forcing_term::Union{ForcingTerm,Nothing}`: Computed forcing term
- `grid_template::G`: Template for grid generation
"""
struct VerificationProblem{P,M,G}
    pde::P
    manufactured_solution::M
    forcing_term::Union{ForcingTerm,Nothing}
    grid_template::G
    
    function VerificationProblem(pde::P, mms::M, grid_template::G=nothing) where {P,M,G}
        new{P,M,G}(pde, mms, nothing, grid_template)
    end
end

# ===================================================================
# User Solver Types
# ===================================================================

"""
    UserSolver{F}

Wrapper for user's numerical solver function.

# Fields
- `solve_function::F`: User's solver function with signature (grid, forcing, params...) -> solution
- `metadata::Dict{Symbol,Any}`: Solver information
"""
struct UserSolver{F}
    solve_function::F
    metadata::Dict{Symbol,Any}
    
    function UserSolver(solve_fn::F; 
                       name::String="User Solver",
                       description::String="") where {F}
        metadata = Dict{Symbol,Any}(
            :name => name,
            :description => description
        )
        new{F}(solve_fn, metadata)
    end
end

# ===================================================================
# Verification Result Types
# ===================================================================

"""
    ConvergenceData

Stores error data at a single resolution.
"""
struct ConvergenceData{T<:Real}
    resolution::T
    grid_points::Int
    errors::Dict{Symbol,T}
    timing::T
end

"""
    VerificationResult{T<:Real}

Complete verification results across multiple resolutions.

# Fields
- `problem_info::Dict{Symbol,Any}`: Information about the problem
- `solver_info::Dict{Symbol,Any}`: Information about the solver
- `convergence_data::Vector{ConvergenceData{T}}`: Data at each resolution
- `convergence_rates::Dict{Symbol,T}`: Estimated convergence rates
"""
struct VerificationResult{T<:Real}
    problem_info::Dict{Symbol,Any}
    solver_info::Dict{Symbol,Any}
    convergence_data::Vector{ConvergenceData{T}}
    convergence_rates::Dict{Symbol,T}
end
