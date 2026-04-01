function solve(inital_condition, method, h, k, x, stopping_time, a, forcing = nothing; deb=false)
    lenX = length(x)

    # Time-stepping loop parameters
    total_steps = floor(Int, stopping_time / k)
    
    results = zeros(lenX, total_steps + 1) 
    u = inital_condition.(x)
    results[:, 1] = u
    
    # Print the initial 2-norm
    if deb
        println("num_steps=$total_steps")
        println("len(x)=$(length(x))")
        println("k=$k")
        println("a=$a")
        println("h=$h")
    end

    # Time-stepping loop
    for n in 1:total_steps
        # Update according to method
        u_new = method(u, x, a, k, h)

        # Add forcing elementwise
        if (forcing !== nothing)
            t = n*k
            for x_spot in 1:lenX
                u_new[x_spot] += k*forcing(x[x_spot],t)
            end
        end 

        results[:, n + 1] = u_new  
        u = u_new  
    end

    return results
end