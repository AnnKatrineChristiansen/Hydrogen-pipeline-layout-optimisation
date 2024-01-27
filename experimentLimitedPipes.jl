include("readData.jl")
include("sweepConstructionHeuristic.jl")
include("writeResults.jl")
import Random
using Combinatorics 


function run(instance, zone, number_of_substations, clockwise, strings, T, n_array_pipes, folder_name)
    
    println("STARTED RUN INSTANCE ", instance, " ZONE ", zone)
    
    # Set seed to reproduce results 
    Random.seed!(678)

    # Get all combinations of pipes
    pipes = collect(1:n_array_pipes)
    combinationer = collect(combinations(pipes, T))

    results = []
    best_sol = nothing
    best_cal = nothing
    best_cost = Inf
    best_combination = nothing

    for combination in combinationer
        
        # Get data
        coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone), combination)

        # Run Sweep Construction heuristic
        substations = collect(1:number_of_substations)
        start_turbines = collect((number_of_substations+1):data.n)
        sweep_sol, sweep_cal = sweepConstructionHeuristic.run(data, coords, substations, start_turbines, clockwise, strings)

        if sweep_sol.cost < best_cost
            best_cost = deepcopy(sweep_sol.cost)
            best_sol = deepcopy(sweep_sol)
            best_cal = deepcopy(sweep_cal)
            best_combination = deepcopy(combination)
        end
        
        # Save solution
        WriteResults.create_folder(folder_name)
        WriteResults.write_json(sweep_sol, sweep_cal, data, folder_name*"/"*instance*zone*"_"*"sweep"* string([element + 1 for element in combination]))

        row = [instance, zone, sweep_sol.cost, string([element + 1 for element in combination]), sweep_cal.feasible, sweep_sol.crossings,  sweep_sol.m_connections, sweep_sol.t_connections]
        push!(results, row)

        println("DONE WITH COMBINATION ", [element + 1 for element in combination])

    end

    WriteResults.write_txt(folder_name*"/"*instance*zone*"_results_"*string(T), results, "exp_limited_pipes")

    best_row = [instance, zone, best_sol.cost, [element + 1 for element in best_combination], best_cal.feasible, best_sol.crossings,  best_sol.m_connections, best_sol.t_connections]
    WriteResults.write_to_exsiting_txt(folder_name*"/all_instances", best_row, false, true)

    println("BEST COMBINATION OF PIPES ", [element + 1 for element in best_combination], " HAS COST ", best_cost)

end 

#= instance = "B"
zone = "1"
number_of_substations = 2
clockwise = true
strings = false
T = 2 # Number of pipes to choose
n_array_pipes = 6
folder_name = "12-12_C_B_exp_limited_pipes_2" =#

#run(instance, zone, number_of_substations, clockwise, strings, T, n_array_pipes, folder_name)
run(ARGS[1], ARGS[2],parse(Int64,ARGS[3]), parse(Bool, ARGS[4]),  parse(Bool, ARGS[5]), parse(Int64,ARGS[6]), 6, ARGS[7])

