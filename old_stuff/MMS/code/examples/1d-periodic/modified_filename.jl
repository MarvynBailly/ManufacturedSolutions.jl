using Plots, HDF5
include("solver.jl")


# Lax-Friedrichs
function LxF(u, x, a, k, h)
    up = circshift(u, -1)  # Shift u right for u_{i+1}
    um = circshift(u, 1)   # Shift u left for u_{i-1}
    u_new = 0.5 * (up + um) - 0.5 * (a * k) / (h) .* (up - um)
    return u_new
end

# Lax-Wendroff
function LxW(u, x, a, k, h)
    up = circshift(u, -1)  # Shift u right for u_{i+1}
    um = circshift(u, 1)   # Shift u left for u_{i-1}
    u_new = u - (a*k * 0.5) *(up - um)/(h) + (a^2 * k^2 * 0.5)*(up - 2*u + um)/(h^2)
    return u_new
end

# Initial Condition
f(x) = sin.(x)

forcing(x,t) = -sin(x)*exp(t*x)*sin(t) + 4(cos(t)*cos(x)*exp(t*x) + t*cos(t)*sin(x)*exp(t*x)) + x*cos(t)*sin(x)*exp(t*x)
stopping_time = 1.5
h = 0.005
output_filename = "./code/examples/1d-periodic/data.h5"


# Parameters
a = 4
cfl = 0.25 # CFL number
domain_left = 0;
domain_right = 2*pi;  

# Discretize domain according to h
x = range(domain_left, domain_right, step = h)
k = (h * cfl) / a
println("Using time step k = $k")
total_steps = floor(Int, stopping_time / k)
println("Total time steps = $total_steps")

result = zeros(length(x), total_steps + 1) 

# Get the result
result = solve(f, LxW, h, k, x, stopping_time, a, forcing; deb=true)

if isfile(output_filename)
    rm(output_filename)  # Remove the old file
end
file = h5open(output_filename, "w")

#Save data into the HDF5 file
write(file, "data/results", result[:,end])  
write(file, "data/x_values", collect(x))  
close(file)