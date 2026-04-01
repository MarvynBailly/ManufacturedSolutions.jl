include("utils.jl")
using .Utils
using Symbolics, Plots

# Parameters
const DEBUG = true                  #print debug data
const PLOT = true                   #save plots
const T_FINAL = 1.5                 #final time for simulations
const H_VALUES = [0.01,0.005]         #h values to test convergence - large to small!
const FILE_PATH = "./code/examples/1d-periodic/"
const FILE_LOCATION = FILE_PATH*"example1.jl"
# const FILE_LOCATION = FILE_PATH*"example2.jl"
# const FILE_LOCATION = FILE_PATH*"dffuse.jl"
const MODIFIED_FILE = FILE_PATH*"modified_filename.jl"
const OUTPUT_FILENAME = FILE_PATH*"data.h5"

# Define variables
@variables x t u(x, t)

# Define the PDE and manufactured solution
Dx = Differential(x)
Dt = Differential(t)

a = 4
pde = Dt(u) + a * Dx(u)  
modified_solution = exp(x * t) * cos(t) * sin(x)  

# κ = 0.1
# pde = Dt(u) - κ * Dx(Dx(u))
# modified_solution = exp(x * t) * cos(t) * sin(x)



######### Main #########
# Create manufactured solution and forcing term
true_sol, forcing_term_str = create_manufactured_solution(pde, modified_solution, u, x, t)

if DEBUG
    println(forcing_term_str)
end

## Loop through the h values - prob a more efficient way to do this
l2_errors = []

for h in H_VALUES
    @info("Running simulation with h = $h")

    # Modify the script with the forcing term and final time
    modify_script(FILE_LOCATION, OUTPUT_FILENAME, MODIFIED_FILE, T_FINAL,forcing_term_str, h)

    # Run the modified script
    run_script(MODIFIED_FILE)

    # Read the data from HDF5
    x_values, results = read_data(OUTPUT_FILENAME)

    # Compute true solution at t_final
    true_solution = [true_sol([x, T_FINAL]) for x in x_values]


    if PLOT
        @info("plotting results")
        plot(x_values, true_solution, label="True result(x, t=$T_FINAL)", xlabel="x", ylabel="u(x)", linestyle=:dot)
        plot!(x_values, results, label="Approximate result(x, t=$T_FINAL)", xlabel="x", ylabel="u(x)",  legend=:bottomright, linestyle=:dash)
        savefig(FILE_PATH*"manufactured-plot.png")  
    end

    # Compute and display errors
    l2_error = compute_l2_error(true_solution, results, h)
    println(l2_error)
    push!(l2_errors,l2_error)
end

# approximate the order
p, ratio = compute_order(l2_errors, H_VALUES)
println("With a ratio of $ratio, the order of accuracy is $p")