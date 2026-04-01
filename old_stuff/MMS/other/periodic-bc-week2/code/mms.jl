using Symbolics, Plots, HDF5, LinearAlgebra


#### Params

debug = true 
t_final = 1.5
h = 0.001;
file_location = "example_working_mms.jl"

################## Create Manufactured Forcing Term ##################
@info("Beginning symbolic fun")

# Define variables
@variables x,t, u(x,t)

# Define differentials
Dx = Differential(x)
Dt = Differential(t)

# Define the PDE
a = 4
pde = Dt(u) + a*Dx(u)

# Define manufactured solution
manu = exp(x*t)*cos(t)*sin(x)

# The manufactured solution is the true solution to the modified PDE
true_sol = build_function(manu,[x,t], expression=Val{false})

# Substitute manufactured solution into PDE
Q = substitute(pde, u => manu)
# simplified_Q = expand_derivatives(Q)
# println(simplified_Q)

forcing_term = expand_derivatives(Q)

# turn it into a string / Make into whatever we need s.t. we can add it to the lang
forcing_term_str = string(forcing_term)
if debug
    @info("The force terming is given by: Q(x,t)=$forcing_term_str")
end

################## Add Forcing to Script ##################
@info("Adding the forcing term")

# now let's read in the script
script_contents = read(file_location, String)

# and search for '# INSERT_FORCING_TERM_HERE
modified_script_contents = replace(script_contents, "# INSERT_FORCING_TERM_HERE" => "$forcing_term_str")

# AND SERACH FOR '# INSERT_STOP_TIME_HERE
modified_script_contents = replace(modified_script_contents, "# INSERT_STOP_TIME_HERE" => "$t_final")

# AND SERACH FOR '# INSERT_STOP_TIME_HERE
modified_script_contents = replace(modified_script_contents, "# INSERT_H_VALUE_HERE" => "$h")

# write the file with the forcing term
write("modified_script_contents.jl", modified_script_contents)


################## Run the Modified Sciprt ##################
# run the script
@info("Running modified solution")

# Will have to replace with whatever language
run(`julia --project=. modified_script_contents.jl`)
# This is a blocking call so need to worry

################## Read Data ##################
### Question: How does C++ return the data?
@info("reading data")
# Define the file name
filename = "data.h5"

# Read the data from the file 
results = h5read(filename, "data/results")  
x_values = h5read(filename, "data/x_values")  

# Evaluate the manufactured solution for each x value at t 
y_values = [true_sol([x, t_final]) for x in x_values]


################## Analysis ##################
### Question: Best way to do this part? Run the scipt for an increasingly fine mesh ratio and observe the accuracy rate?
# Plot the result
@info("plotting results")
plot(x_values, y_values, label="True result(x, t=$t_final)", xlabel="x", ylabel="u(x)", linestyle=:dot)
plot!(x_values, results, label="Approximate result(x, t=$t_final)", xlabel="x", ylabel="u(x)",  legend=:top, linestyle=:dash)
# plot!(x_values, norm(y_values - results), label="err(x, t=$t_final)", xlabel="x", ylabel="u(x)", lw=2, legend=:top)

savefig("plot.png")  

# Look at the l2 error norm
errs = norm.(y_values - results)
println("The max l2 error is $(maximum(errs))")