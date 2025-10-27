include("utils.jl")
using .Utils
using Symbolics, Plots

# Parameters
const DEBUG = true
const PLOT = true
const T_FINAL = 1.5
const H_VALUE = 0.001
const FILE_LOCATION = "example_working_mms.jl"
const MODIFIED_FILE = "modified_script_contents.jl"
const HDF5_FILENAME = "data.h5"

# Define variables
@variables x t u(x, t)

# Define the PDE and manufactured solution
Dx = Differential(x)
Dt = Differential(t)
a = 4
pde = Dt(u) + a * Dx(u)  
modified_solution = exp(x * t) * cos(t) * sin(x)  


######### Main #########

# Create manufactured solution and forcing term
true_sol, forcing_term_str = create_manufactured_solution(pde, modified_solution, u, x, t)

if DEBUG
    println(forcing_term_str)
end

# Modify the script with the forcing term and final time
modify_script(FILE_LOCATION, forcing_term_str, MODIFIED_FILE, T_FINAL, H_VALUE)

# Run the modified script
run_script(MODIFIED_FILE)

# Read the data from HDF5
x_values, results = read_data(HDF5_FILENAME)

# Compute true solution at t_final
true_solution = [true_sol([x, T_FINAL]) for x in x_values]


if PLOT
    @info("plotting results")
    plot(x_values, true_solution, label="True result(x, t=$T_FINAL)", xlabel="x", ylabel="u(x)", linestyle=:dot)
    plot!(x_values, results, label="Approximate result(x, t=$T_FINAL)", xlabel="x", ylabel="u(x)",  legend=:top, linestyle=:dash)
    savefig("plot.png")  
end

# Compute and display errors
l2_error = compute_l2_error(true_solution, results)


println(l2_error)