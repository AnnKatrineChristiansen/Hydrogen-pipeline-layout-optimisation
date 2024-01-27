include("readData.jl")
include("sweepConstructionHeuristic.jl")
include("writeResults.jl")
include("graph.jl")
include("improvementHeuristic.jl")
import Random
using JSON
global FOLDER_PATH = "results/" # "solutions/"


#instances = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
#zones = [[1, 2], [1], [1, 2], [1], [1, 2, 3], [1, 2], [1], [1, 2], [1], [1, 2, 3]]
#substations_per_zones = [[1 ,1] ,[2] ,[1 ,1] ,[3] ,[1 ,2 ,1] ,[2 ,1] ,[3] ,[2 ,1] ,[4] ,[2 ,2, 1]]

 instance = "A" #"B" #"A" 
zone = "1"
number_of_substations = 1 #1
clockwise = false
strings = false
balacing = false
digging_cost = 0
folder_name = "16-01_Sweep_plot"

function elapsed_time(start_time)
    return round((time_ns()-start_time)/1e9,digits = 3)
end 


function write_crossings_to_exsiting_txt(filname, rows)

    if !isfile(FOLDER_PATH * filname * "/" *"grupper.txt")
        # The file doesn't exist, so create it
        file = open(FOLDER_PATH * filname * "/" *"grupper.txt", "w")
        
        # Optionally, write initial content to the file
       
        header = ["Substation", "StartingTurbine", "GroupTurbine"]

        header_str = join(string.(header), " ")  # Convert each element to a string
        write(file, header_str * "\n")
        
        # Close the file when you're done
        close(file)
    
    end 
     
    # Open the file in append mode to add a new line
    file = open(FOLDER_PATH * filname * "/" *"grupper.txt", "a")

    # Text you want to write as a new line

    for row in rows
        row_str = join(string.(row), " ")  # Convert each element to a string
        write(file, row_str * "\n")
    end

    # Close the file when you're done
    close(file)

end



function results_sweep(instance, zone, number_of_substations, clockwise, strings, balancing, folder_name, digging_cost)
    Random.seed!(678)
    start_time = time_ns()

    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    data.digging_cost = digging_cost
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):(data.n)) # collect((19):(19))
    sweep_sol, sweep_cal,  worst_sol, worst_cal, grupper  = sweepConstructionHeuristic.run_sweep_implementation(data, coords, substations, start_turbines, clockwise, strings, balancing)

    runtime = elapsed_time(start_time)

    row = [instance, zone, sweep_sol.cost, worst_sol.cost, sweep_sol.crossings, sweep_sol.m_connections, sweep_sol.t_connections, runtime]

    WriteResults.create_folder(folder_name)
    WriteResults.write_json(sweep_sol, sweep_cal, data, folder_name*"/"*instance*zone*"_"*"sweep_best")
    WriteResults.write_json(worst_sol, worst_cal, data, folder_name*"/"*instance*zone*"_"*"sweep_worst")
    WriteResults.write_to_exsiting_txt(folder_name*"/all_instances", row, false, false, true)

    #println(grupper)
    write_crossings_to_exsiting_txt(folder_name, grupper)

end 

results_sweep(instance, zone, number_of_substations, clockwise, strings, balacing, folder_name, digging_cost)
#results_sweep(ARGS[1], ARGS[2],parse(Int64,ARGS[3]), parse(Bool, ARGS[4]),  parse(Bool, ARGS[5]), parse(Bool, ARGS[6]), ARGS[7],  parse(Float64,ARGS[8]))


function evaluate_distance(instance, zone, number_of_substations, clockwise, strings, balancing, folder_name)
    
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):(data.n)) # collect((19):(19))
    sweep_sol, sweep_cal = sweepConstructionHeuristic.run(data, coords, substations, start_turbines, clockwise, strings, balancing)

    #println("Final solution")
    #cost = improvementHeuristic.cost_function(sweep_sol, data, substations)

    WriteResults.create_folder(folder_name)
    WriteResults.write_json(sweep_sol, sweep_cal, data, folder_name*"/"*instance*zone*"_"*"sweep")


    total = 0
    max_dist = 0
    for edge in sweep_sol.edges
        i, j = improvementHeuristic.find_edge_from_index(data, edge)
        dist = data.distance[i,j]
        if  dist > max_dist 
            max_dist = dist
        end 

        total += dist 

    end 

    println("Total distance ", total)
    println("Average distance ", total/data.n)
    println("Max distance ", max_dist)

end 

# Instance A Zone 1 Branching 
#Total distance 35.63932165560401
#Average distance 1.3199748761334817
#Max distance 2.652247535502342


# Instance B Zone 1 Branching 