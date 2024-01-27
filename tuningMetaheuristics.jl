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

 

function run(instance, zone,  number_of_substations, clockwise, strings, timelimit, folder_name, heuristic)
    
    println("STARTED RUN INSTANCE ", instance, " ZONE ", zone)

    # Set seed to reproduce results 
    Random.seed!(678)

    # Get data
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))
    balancing = false
    data.digging_cost = 0

    # Run Sweep Construction heuristic
    println("START SWEEP")
    substations = collect(1:number_of_substations)
    start_turbines = collect((number_of_substations+1):data.n)
    sweep_sol, sweep_cal = sweepConstructionHeuristic.run(data, coords, substations, start_turbines, clockwise, strings, balancing)
    println("DONE SWEEP")

    println("START TUNING")
    WriteResults.create_folder(folder_name)
    # Run Improvement Heuristic
    if heuristic == "SA"
        #["Heuristic", "Elapsed Time","Best cost", "Feasible", "Temperature", "Alpha"]
        results_SA = [["Sweep", 0, sweep_sol.cost, sweep_cal.feasible, 0, 0]]

        reheats = [Inf, 5000, 10000, 50000]
        
        for reheat in reheats
            start_time = time_ns()
            improved_sol, improved_cal, results, Tstart, alpha = improvementHeuristic.simulated_annealing(sweep_sol, sweep_cal, data, substations, strings, start_time, timelimit, reheat, balancing)
            runtime = elapsed_time(start_time)
            # SAVE INFORMATION
            SA_row = ["SA", runtime, improved_sol.cost, improved_cal.feasible, Tstart, alpha]
            push!(results_SA, SA_row)
            WriteResults.write_txt(folder_name*"/"*instance*zone*"_"*heuristic*"_"*string(Tstart)*"_"*string(alpha)*"_"*string(reheat), results, heuristic)
            println("COST ", improved_sol.cost, " with Tstart ", Tstart, ", alpha ", alpha, "and reheat ", reheat)
        end 


        WriteResults.write_txt(folder_name*"/"*instance*zone*"_tuning_"*heuristic, results_SA, "TuningSA") 

    elseif heuristic == "VNS"
        Kmaxs = [2, 3, 4, 5, 6]
        first_improves = [false, true]
        allow_illegals = [false, true]
        # ["Heuristic", "Elapsed Time","Best cost", "Feasible","K", "First Improvement"]
        results_VNS = [["Sweep", 0, sweep_sol.cost, sweep_cal.feasible, 0, 0]]
        for Kmax in Kmaxs
            for first_improve in first_improves
                for allow_illegal in allow_illegals
                    start_time = time_ns()
                    improved_sol, improved_cal, results = improvementHeuristic.VNS(sweep_sol, sweep_cal, data, substations, start_time, timelimit, Kmax, first_improve, strings, allow_illegal, balancing)
                    runtime = elapsed_time(start_time)
                    # SAVE INFORMATION
                    VNS_row = ["VNS", runtime, improved_sol.cost, improved_cal.feasible, Kmax, first_improve, allow_illegal]
                    push!(results_VNS, VNS_row)
                    WriteResults.write_txt(folder_name*"/"*instance*zone*"_"*"iterations"*"_"*heuristic*"_"*string(Kmax)*"_"*string(first_improve)*"_"*string(allow_illegal), results, heuristic) 
                    println("COST ", improved_sol.cost, " with first improvement ", first_improve, " and Kmax ", Kmax, " with allow illegal ", allow_illegal)
                end 
            end 
        end 

        WriteResults.write_txt(folder_name*"/"*instance*zone*"_tuning_"*heuristic, results_VNS, "TuningVNS") 
    end 


    println("DONE RUNNING TUNING ON INSTANCE ", instance, " ZONE ", zone, " USING ", heuristic)

end 

#= instance = "B"
zone = "1"
number_of_substations = 2
clockwise = true
strings = false
timelimit = 60*5
solution_filename = "11-22_C_B_tuning_test"
heuristic = "VNS"  =#

#run(instance, zone, number_of_substations, clockwise, strings, timelimit, solution_filename, heuristic)
run(ARGS[1], ARGS[2],parse(Int64,ARGS[3]), parse(Bool, ARGS[4]),  parse(Bool, ARGS[5]), parse(Int64,ARGS[6]), ARGS[7], ARGS[8])