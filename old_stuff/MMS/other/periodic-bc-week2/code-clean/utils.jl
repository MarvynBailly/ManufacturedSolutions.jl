module Utils
# make it a module for faster? 

using Symbolics, Plots, HDF5, LinearAlgebra

# add an export cause maybe it makes it faster?
export create_manufactured_solution, modify_script, run_script, read_data, compute_l2_error, compute_order


# Function to set up manufactured solution and forcing term
function create_manufactured_solution(pde, manu, u, x, t)
    @info("Creating manufactured solution and forcing term")
    # Define the true solution
    true_sol = build_function(manu, [x, t], expression = Val{false})
    
    # Substitute manufactured solution into the provided PDE
    forcing_term = expand_derivatives(substitute(pde, u => manu))
    
    return true_sol, string(forcing_term)
end


# Function to modify the script
function modify_script(file_location, forcing_term_str, MODIFIED_FILE, t_final, h_value)
    @info("Modifying the script with forcing term")
    script_contents = read(file_location, String)
    modified_script = replace(script_contents, "# INSERT_FORCING_TERM_HERE" => forcing_term_str)
    modified_script = replace(modified_script, "# INSERT_STOP_TIME_HERE" => "$t_final")
    modified_script = replace(modified_script, "# INSERT_H_VALUE_HERE" => "$h_value")
    write(MODIFIED_FILE, modified_script)
end


# Function to run the modified script
function run_script(MODIFIED_FILE)
    @info("Running modified script")
    run(`julia --project=. $MODIFIED_FILE`)
end


# Function to read data from HDF5
function read_data(filename)
    @info("Reading data from HDF5 file")
    results = h5read(filename, "data/results")
    x_values = h5read(filename, "data/x_values")
    return x_values, results
end


# Function to compute l2 errors
function compute_l2_error(true_solution, results)
    @info("Computing l2 error")
    errs = norm(true_solution - results, 2)
    l2_error = maximum(errs)
    return l2_error
end

# Function to compute the order
function compute_order(errors, h_values)
    ratio = h_values[1]/h_values[2]
    p = log(errors[1]/errors[2]) / log(ratio)
    return p, ratio
end


end