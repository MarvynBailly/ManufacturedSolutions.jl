using LinearAlgebra
using SparseArrays

# Parameters and grid setup
Nx, Ny = 100, 100  # grid resolution
Lx, Ly = 10.0, 10.0  # physical dimensions
dx, dy = Lx / Nx, Ly / Ny  # grid spacing
dt = 0.01  # time step
T = 1.0  # total simulation time
ν = 1.0  # kinematic viscosity
β = 1.0  # parameter beta

# Initialize variables
p = zeros(Nx, Ny)
u = zeros(Nx, Ny)
v = zeros(Nx, Ny)
Sp = zeros(Nx, Ny)
Su = zeros(Nx, Ny)
Sv = zeros(Nx, Ny)

function compute_rhs!(u, v, p, ν, dx, dy, Sp, Su, Sv)
    # Initialize RHS arrays for u and v
    Nx, Ny = size(u)
    RHSu = zeros(Nx, Ny)
    RHSv = zeros(Nx, Ny)

    for i in 1:Nx, j in 1:Ny
        # Apply periodic indexing
        ip = (i % Nx) + 1     # Right neighbor
        im = (i - 2 + Nx) % Nx + 1  # Left neighbor
        jp = (j % Ny) + 1     # Top neighbor
        jm = (j - 2 + Ny) % Ny + 1  # Bottom neighbor

        # Advection terms for u and v (central difference)
        adv_u = (u[i, j] * (u[ip, j] - u[im, j]) / (2 * dx) +
                 v[i, j] * (u[i, jp] - u[i, jm]) / (2 * dy))
        adv_v = (u[i, j] * (v[ip, j] - v[im, j]) / (2 * dx) +
                 v[i, j] * (v[i, jp] - v[i, jm]) / (2 * dy))

        # Pressure gradient terms for u and v (central difference)
        gradp_u = (p[ip, j] - p[im, j]) / (2 * dx)
        gradp_v = (p[i, jp] - p[i, jm]) / (2 * dy)

        # Diffusion terms for u and v (Laplacian using central difference)
        diff_u = ν * ((u[ip, j] - 2 * u[i, j] + u[im, j]) / dx^2 +
                      (u[i, jp] - 2 * u[i, j] + u[i, jm]) / dy^2)
        diff_v = ν * ((v[ip, j] - 2 * v[i, j] + v[im, j]) / dx^2 +
                      (v[i, jp] - 2 * v[i, j] + v[i, jm]) / dy^2)

        # Construct RHS for u and v by combining all terms
        RHSu[i, j] = -adv_u - gradp_u + diff_u + Su[i, j]
        RHSv[i, j] = -adv_v - gradp_v + diff_v + Sv[i, j]
    end

    return RHSu, RHSv
end


# Function to apply doubly periodic boundary conditions
function apply_periodic!(var)
    var[:, end] = var[:, 1]  # y-periodicity
    var[end, :] = var[1, :]  # x-periodicity
end

# Build Crank-Nicolson matrix for diffusion
function build_diffusion_matrix(Nx, Ny, dx, dy, dt, ν)
    D = spzeros(Nx * Ny, Nx * Ny)  # Sparse matrix for Crank-Nicolson

    # Fill the matrix with finite difference stencil
    for i in 1:Nx, j in 1:Ny
        idx = (i - 1) * Ny + j  # Linear index

        # Center point
        D[idx, idx] = 1 + dt * ν * (2 / dx^2 + 2 / dy^2)

        # Neighbor points with periodic boundary conditions
        D[idx, ((i % Nx) + 1 - 1) * Ny + j] = -0.5 * dt * ν / dx^2  # right neighbor
        D[idx, ((i - 2 + Nx) % Nx + 1 - 1) * Ny + j] = -0.5 * dt * ν / dx^2  # left neighbor
        D[idx, (i - 1) * Ny + (j % Ny) + 1] = -0.5 * dt * ν / dy^2  # top neighbor
        D[idx, (i - 1) * Ny + (j - 2 + Ny) % Ny + 1] = -0.5 * dt * ν / dy^2  # bottom neighbor
    end

    return D
end

# Build the diffusion matrix once
D = build_diffusion_matrix(Nx, Ny, dx, dy, dt, ν)

u .= 0.01 * rand(Nx, Ny)  # random small velocity perturbations
v .= 0.01 * rand(Nx, Ny)  # random small velocity perturbations
p .= 0.0  # zero initial pressure


# Time-stepping loop
for t in 0:dt:T
    # Compute RHS for u and v at the current time step (n)
    RHSu, RHSv = compute_rhs!(u, v, p, ν, dx, dy, Sp, Su, Sv)

    # Flatten RHS for Crank-Nicolson
    rhs_u = reshape(u .+ 0.5 * dt .* RHSu, Nx * Ny)
    rhs_v = reshape(v .+ 0.5 * dt .* RHSv, Nx * Ny)

    # Solve for u and v at the next time step (n+1) using the Crank-Nicolson matrix
    u_next = D \ rhs_u
    v_next = D \ rhs_v

    # Reshape the results back to 2D arrays
    u .= reshape(u_next, Nx, Ny)
    v .= reshape(v_next, Nx, Ny)

    # Apply periodic boundary conditions
    apply_periodic!(u)
    apply_periodic!(v)
end

using Plots
print(u)
heatmap(u, title="Velocity Field u", xlabel="x", ylabel="y", color=:viridis, xticks=0:20:Nx, yticks=0:20:Ny)
savefig("u")
heatmap(v, title="Velocity Field v", xlabel="x", ylabel="y", color=:plasma, xticks=0:20:Nx, yticks=0:20:Ny)
savefig("v")