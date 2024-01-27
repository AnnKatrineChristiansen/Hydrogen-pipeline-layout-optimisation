include("readData.jl")
include("sweepConstructionHeuristic.jl")
include("writeResults.jl")
include("graph.jl")
include("improvementHeuristic.jl")
include("pressureDropCalculations.jl")
include("Kmeans.jl")
import Random
using JSON
using DataFrames
using CSV



function run_test_with_general_pressures(instance, zone,  number_of_substations, clockwise, strings, folder_name)
    
    println("STARTED RUN INSTANCE ", instance, " ZONE ", zone)

    # Get data
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    export_ids, export_rough , export_costs = readData.readExportPipeinfo("Data/Pipeline info and cost data.xlsx")

    # To save results 
    results = DataFrame(start_pressures = Float64[], solution_cost = Float64[],)

    # Run Sweep Construction heuristic for different pressures
    #println("START SWEEP")
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):data.n)

    # Range from 10 to 40
    pressure_range = 10:35

    for start_pressure in pressure_range

        println("Manifold pressures: ", start_pressure)

        sweep_sol, sweep_cal = sweepConstructionHeuristic.run_pressure_experiment(data, coords, substations, start_turbines, ones(number_of_substations)*start_pressure, clockwise, strings)
        
        if sweep_cal.feasible
            push!(results, (start_pressure, sweep_sol.cost))
        else
            push!(results, (start_pressure, NaN))
        end 

    end

    WriteResults.create_folder(folder_name)
    CSV.write(WriteResults.FOLDER_PATH * folder_name*"/"*instance*zone*"_general_pressures.csv", results) 
end

function run_test_with_pipes(instance, zone,  number_of_substations, clockwise, strings, dists_from_shore, pressure_at_shore, folder_name)
    
    println("STARTED RUN INSTANCE ", instance, " ZONE ", zone)

    # Get data
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    export_ids, export_rough , export_costs = readData.readExportPipeinfo("Data/Pipeline info and cost data.xlsx")

    # Save results
    results = DataFrame(shore_pressure = Float64[], dist_from_shore = Float64[], inner_diameter = Float64[], start_pressures = Vector{Float64}[], export_velocities = Vector{Float64}[], average_velocities = Vector{Float64}[], solution_cost = Float64[], export_cost_dist = Float64[])

    # Run Sweep Construction heuristic for different pipes 
    println("START SWEEP")
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):data.n)
    for pressure in pressure_at_shore
        for dist_from_shore in dists_from_shore
            for pipe in 1:length(export_ids)
                # Calculate the start pressure and velocity for each pipe
                inner_dia = export_ids[pipe]
                internal_rough = export_rough[pipe]
                export_cost = export_costs[pipe]
                start_pressures, export_velocities, average_velocities = pressure_at_manifold(coords, data, substations, dist_from_shore, pressure, inner_dia, internal_rough)

                #println("With export pipe diameter: ", inner_dia/25.4, " inches")
                #println("Manifold pressures: ", start_pressures)
                #println("Export velocities: ", export_velocities)

                # Check if pipe can handle pressure and velocity constraint fo each manifold
                feasible_pressure = all(y -> y < data.max_pressure, start_pressures)
                feasible_velocity = all(y -> y < data.max_velocity, export_velocities)  
                # If pipe complies with constraint run experiment      
                if feasible_pressure
                    sweep_sol, sweep_cal = sweepConstructionHeuristic.run_pressure_experiment(data, coords, substations, start_turbines, start_pressures, clockwise, strings)
                    push!(results, (pressure, dist_from_shore, round(inner_dia/25.4), start_pressures, export_velocities, average_velocities, sweep_sol.cost, export_cost*dist_from_shore*1000))
                else 
                    push!(results, (pressure, dist_from_shore, round(inner_dia/25.4), start_pressures, export_velocities, average_velocities, NaN, NaN))
                end

            end
        end
    end

    WriteResults.create_folder(folder_name)
    CSV.write(WriteResults.FOLDER_PATH * folder_name*"/"*instance*zone*"_pipe_pressures.csv", results) 
end


function pressure_at_manifold(coords, data, substations, dist_from_shore, pressure_at_shore, inner_dia, internal_rough)
    
    # Cluster turbines 
    _, clusters = Kmeans.run(data, coords, substations)
    # Define constants 
    c = pressureDropCalculations.initalize_constants()
    # Calculate flow rate in export pipe for each cluster
    flowrates = [length(cluster) for cluster in clusters]*data.turbine_flow_rate_daily
    start_pressures = []
    export_velocities = []
    average_velocities = []
    
    # For each cluster calculate the pressure and velocity of inputted pipe
    for i in 1:length(clusters)
        manifold_pressure, export_velocity = pressureDropCalculations.calculate_inlet_pressure_and_velocity(c, dist_from_shore, pressure_at_shore, inner_dia, flowrates[i], data.temperature, internal_rough)
        average_pressure = (manifold_pressure+pressure_at_shore)/2
        _, average_velocity = pressureDropCalculations.calculate_inlet_pressure_and_velocity(c, dist_from_shore, average_pressure, inner_dia, flowrates[i], data.temperature, internal_rough)
        push!(start_pressures, manifold_pressure)
        push!(export_velocities, export_velocity)
        push!(average_velocities, average_velocity)
    end

    return start_pressures, export_velocities, average_velocities
end

# 
pressure_at_shore = [10,15,20] #bar
dists_from_shore = [25,50]#[25,50,100] #km
instance = "B" 
zone = "1"
number_of_substations = 2
clockwise = true
strings = false
#solution_filename = "Experiment_manifold_pressure" 
solution_filename = "Experiment_ExportPipe"
#


#run_test_with_general_pressures(ARGS[1], ARGS[2],parse(Int64,ARGS[3]), parse(Bool, ARGS[4]),  parse(Bool, ARGS[5]), ARGS[6])

#run_test_with_general_pressures(instance, zone, number_of_substations, clockwise, strings, solution_filename)
run_test_with_pipes(instance, zone, number_of_substations, clockwise, strings, dists_from_shore, pressure_at_shore, solution_filename)