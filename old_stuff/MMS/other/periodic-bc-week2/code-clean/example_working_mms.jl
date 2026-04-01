include("solver.jl")
using Plots, LaTeXStrings, HDF5



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

forcing(x,t) = # INSERT_FORCING_TERM_HERE

# Parameters
a = 4
stopping_time = # INSERT_STOP_TIME_HERE
cfl = 0.25 # CFL number
domain_left = 0;
domain_right = 2*pi;  

H_values = # INSERT_H_VALUE_HERE
N_values = H_values^(-1)

result_collection = []

for i in eachindex(N_values)
    # Discretize domain according to h
    x = range(domain_left, domain_right, step = H_values[i])
    k = (H_values[i] * cfl) / a

    # Get the results
    results = solve(f, LxW, H_values[i], k, x, stopping_time, a, forcing)

    push!(result_collection,results)
end


#plot results
i = 1
result = result_collection[i]
x = range(domain_left, domain_right, step = H_values[i])
# plot(x, result[:,end], label="LxF (n=$(N_values[i]), cfl=$cfl)", lw=2, linestyle=:dash)

# savefig("plot_approx.png")  

filename = "data.h5"
if isfile(filename)
    rm(filename)  # Remove the old file
end
file = h5open(filename, "w")

#Save data into the HDF5 file
write(file, "data/results", result[:,end])  
write(file, "data/x_values", collect(x))  
close(file)