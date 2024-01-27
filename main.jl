include("readData.jl")
include("sweepConstructionHeuristic.jl")
include("writeResults.jl")
include("graph.jl")
include("improvementHeuristic.jl")
import Random
using JSON

function elapsed_time(start_time)
    return round((time_ns()-start_time)/1e9,digits = 3)
end 

function run(instance, zone,  number_of_substations, clockwise, strings, timelimit, folder_name, heuristic, balancing = false, digging_cost = 0)
    
    println("STARTED RUN INSTANCE ", instance, " ZONE ", zone)
    # Start time 
    start_time = time_ns()

    # Set seed to reproduce results 
    Random.seed!(678)

    # Get data
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    data.digging_cost = digging_cost


    #p1 = coords[20, :]
    #p2 = coords[21, :]
    #p3 = coords[1, :]
    #p4 = coords[19, :]
    #crossing = improvementHeuristic.do_line_segments_intersect(p1, p2, p3, p4; tol=1e-8)

    # Run Sweep Construction heuristic
    println("START SWEEP")
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):data.n)
    sweep_sol, sweep_cal = sweepConstructionHeuristic.run(data, coords, substations, start_turbines, clockwise, strings, balancing)
    println("DONE SWEEP")

    println("START IMPROVEMENT")
    # Run Improvement Heuristic
    if heuristic == "SA"
        reheat = Inf
        improved_sol, improved_cal, results = improvementHeuristic.simulated_annealing(sweep_sol, sweep_cal, data, substations, strings, start_time, timelimit, reheat, balancing)
        crossings = improvementHeuristic.check_crossings(improved_sol, data, instance, zone)
    
    elseif heuristic == "VNS"
        Kmax = 4
        first_improve = true
        allow_illegals = false
        improved_sol, improved_cal, results = improvementHeuristic.VNS(sweep_sol, sweep_cal, data, substations, start_time, timelimit, Kmax, first_improve, strings, allow_illegals, balancing)
        crossings = improvementHeuristic.check_crossings(improved_sol, data, instance, zone)
    
    end 

    runtime = elapsed_time(start_time)
    println("DONE RUNNING INSTANCE ", instance, " ZONE ", zone, " USING ", heuristic, " BEST COST ", improved_sol.cost, " IN TIME ", runtime)

    # Write results to files 
    WriteResults.create_folder(folder_name)
    WriteResults.write_json(sweep_sol, sweep_cal, data, folder_name*"/"*instance*zone*"_"*"sweep")
    WriteResults.write_json(improved_sol, improved_cal, data, folder_name*"/"*instance*zone*"_"*heuristic) 
    WriteResults.write_txt(folder_name*"/"*instance*zone*"_"*"iterations"*"_"*heuristic, results, heuristic) 
    WriteResults.write_crossings_to_exsiting_txt(folder_name*"/"*"crossings", crossings)

    sweep_LB = improvementHeuristic.cost_function_bounds(sweep_sol, data, "min")
    sweep_UB = improvementHeuristic.cost_function_bounds(sweep_sol, data, "max")

    row = [instance, zone, sweep_sol.cost, sweep_UB, sweep_LB, heuristic, improved_sol.cost, timelimit, runtime, improved_sol.crossings, improved_sol.m_connections, improved_sol.t_connections, (maximum(values(improved_sol.balanced)) - minimum(values(improved_sol.balanced))), improved_sol.digging_cost]
    WriteResults.write_to_exsiting_txt(folder_name*"/all_instances", row)

end 
 
#= instance = "G" 
zone = "1"
number_of_substations = 3
clockwise = true
strings = false
timelimit = 60*30
solution_filename = "12-17_C_B_SA_30min_find_fejl"
heuristic = "SA"   
balacing = false # CHANGE TO TRUE TO ADD BALANCING
digging_price = 0 # CHANGE DIGGING PRICE =#


#run(instance, zone, number_of_substations, clockwise, strings, timelimit, solution_filename, heuristic, balacing, digging_price)
run(ARGS[1], ARGS[2],parse(Int64,ARGS[3]), parse(Bool, ARGS[4]),  parse(Bool, ARGS[5]), parse(Int64,ARGS[6]), ARGS[7], ARGS[8],  parse(Bool, ARGS[9]), parse(Float64,ARGS[10]))