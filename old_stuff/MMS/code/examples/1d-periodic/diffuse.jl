# Import packages
using Oceananigans
using HDF5
using Plots


# Set simulation parameters
domain_z = (-π, π)     
Nz = 128               
h = # INSERT_H_VALUE_HERE
stop_time = # INSERT_STOP_TIME_HERE      
output_filename = # INSERT_FILE_NAME_HERE 

κ = 0.1

Nz = Int((domain_z[2] - domain_z[1]) / h)
Δz = (domain_z[2] - domain_z[1]) / Nz  
Δt = 0.5 * h^2 / κ    

# Create the grid with periodic boundary conditions
grid = RectilinearGrid(size=(Nz,), z=domain_z, topology=(Periodic,))

# Create the closure with specified diffusivity
closure = ScalarDiffusivity(κ=κ)

# Create the model
model = NonhydrostaticModel(; grid, closure, tracers=:T, timestepper=TimeSteppers.RungeKutta3)

# Define the forcing function with placeholder
function forcing_T(c, x, y, z, t)
    return # INSERT_FORCING_TERM_HERE
end

# Apply the forcing to tracer T
model.forcings.T = FieldForcing(forcing_T)

# Set the initial condition to match the manufactured solution at t=0
width = 0.1
initial_temperature(z) = exp(-z^2 / (2width^2))
set!(model, T=initial_temperature)

# Set up and run the simulation
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# Run the simulation
run!(simulation)

# Extract the results at the final time
T_num = model.tracers.T.data[:, 1, 1]  # Numerical solution
z_values = grid.zC                     # Spatial coordinates in z

# Save results to an HDF5 file
if isfile(output_filename)
    rm(output_filename)  # Remove the old file if it exists
end
file = h5open(output_filename, "w")

# Save data to HDF5
write(file, "data/results", T_num)  
write(file, "data/x_values", z_values)  
close(file)

println("Results saved to $output_filename.")
