include("solver.jl")
using Plots, LaTeXStrings



# Lax-Friedrichs
function LxF(u, x, a, k, h)
    up = circshift(u, -1)  # Shift u right for u_{i+1}
    um = circshift(u, 1)   # Shift u left for u_{i-1}
    u_new = 0.5 * (up + um) - 0.5 * (a * k) / (h) .* (up - um)
    return u_new
end



# Initial Condition
f(x) = sin.(x)

forcing(x,t) = sin(t) .+ 4*cos.(x)

# Parameters
a = 4
stopping_time = 1.5
cfl = 0.85 # CFL number
domain_left = 0;
domain_right = 2*pi;  

N_values = [10]
H_values = N_values.^(-1)

result_collection = []

for i in eachindex(N_values)
    # Discretize domain according to h
    x = range(domain_left, domain_right, step = H_values[i])
    k = (H_values[i] * cfl) / a

    # Get the results
    results = solve(f, LxF, H_values[i], k, x, stopping_time, a, forcing)

    push!(result_collection,results)
end


#plot results
i = 1
result = result_collection[i]
x = range(domain_left, domain_right, step = H_values[i])
plot(x, result[:,end], label="LxF (n=$(N_values[i]), cfl=$cfl)", lw=2, linestyle=:dash)

savefig("plot.png")  