module sweepConstructionHeuristic
include("sweep.jl")
include("minimumSpanningTree.jl")
include("graph.jl")
include("pressureDropCalculations.jl")
include("Kmeans.jl")
include("nearestNeighbor.jl")
include("improvementHeuristic.jl")


function append_vectors(a, b)
    append!(a,b)
    return a
end 

function run(data, coords, substations, starting_turbines, clockwise = true, string = false, balancing = false)
    # Assign turbines to a substation 
    centroids, clusters = Kmeans.run(data, coords, substations) # Currently ssuming fixed substations

    best_sol = []
   
    # for each substation find the best solution from different starting points 
    for(idx, substation) in enumerate(substations)
        cluster_polar_angles = sweep.calculate_polar_angles(coords, clusters[idx], substation)

        best_sol_cluster = nothing
        best_cal_cluster = nothing
        best_cost_cluster = Inf
        worst_cost_cluster = 0
        best_node_cluster  = nothing 

        cluster_starting_turbines = [i for i in starting_turbines if i in clusters[idx]]
        for start_turbine in cluster_starting_turbines

            sol_cluster, cal_cluster =  perform_sweep_on_cluster(data, cluster_polar_angles, substation, start_turbine, clockwise, string, balancing)

            if sol_cluster.cost < best_cost_cluster
                best_cost_cluster = sol_cluster.cost
                best_sol_cluster = sol_cluster
                best_cal_cluster = cal_cluster 
                best_node_cluster = start_turbine
            end 
    
            if sol_cluster.cost > worst_cost_cluster
                worst_cost_cluster = sol_cluster.cost
            end 

        end

        push!(best_sol, best_sol_cluster)

        #println("For substation ", substation, " the best solution starts at node ", best_node_cluster, ". It has a cost of ", best_cost_cluster, " which is ", worst_cost_cluster - best_cost_cluster, " than the worst starting point."  )

    end 

    sol_combined = minimumSpanningTree.solution()
    for sol in best_sol
        sol_combined.adjacency_list = merge!(append_vectors, sol_combined.adjacency_list, sol.adjacency_list)
    end 

    final_sol, final_cal =  pressureDropCalculations.run(sol_combined, data, substations)
    final_sol.cost, final_sol.digging_cost  = improvementHeuristic.cost_function(final_sol, data, substations, final_cal.feasible, balancing)

    return final_sol, final_cal, grupper  
    
end 

function perform_sweep_on_cluster(data, polar_angles, substation, starting_turbine, clockwise = true, strings = false, balancing = false)
    ordered_polar_angles = sweep.order_polar_angles(polar_angles, starting_turbine, clockwise)
    node_list = [substation]
    sol_list = Set()
    previous_sol = nothing
    previous_cal = nothing
    num = 1

    while num <= length(ordered_polar_angles)
        #start_time = time_ns()
        append!(node_list, ordered_polar_angles[num])

        if strings
            sol = nearestNeighbor.run(data, node_list, substation)
        else 
            g = Graph.create(data, node_list)
            sol = minimumSpanningTree.run(g, substation)
            sol.adjacency_list = Graph.convert_to_directed_graph(substation, sol.adjacency_list)
        
        end
        sol, cal  = pressureDropCalculations.run(sol, data, substation)

        # Feasibility check 
        if cal.feasible
            previous_sol = deepcopy(sol)
            previous_cal = deepcopy(cal)
            num += 1

        else
            push!(sol_list, previous_sol)
            node_list = [substation]
        end 

    end 

    push!(sol_list, previous_sol)
    
    # Merge sols from sol_list
    sol_cluster = minimumSpanningTree.solution()
    for sub_sol in sol_list
        sol_cluster.adjacency_list = merge!(append_vectors, sol_cluster.adjacency_list, sub_sol.adjacency_list)
    end 

    sol_final_cluster, cal_final_cluster =  pressureDropCalculations.run(sol_cluster, data, substation)
    sol_final_cluster.cost, sol_final_cluster.digging_cost = improvementHeuristic.cost_function(sol_final_cluster, data, [substation], cal_final_cluster.feasible, balancing)


    return sol_final_cluster, cal_final_cluster

end

function run_pressure_experiment(data, coords, substations, starting_turbines, start_pressures, clockwise = true, string = false)
    # Assign turbines to a substation 
    centroids, clusters = Kmeans.run(data, coords, substations) # Currently ssuming fixed substations

    best_sol = []
   
    # for each substation find the best solution from different starting points 
    for(idx, substation) in enumerate(substations)
        data.start_outlet_pressure = start_pressures[idx]
        cluster_polar_angles = sweep.calculate_polar_angles(coords, clusters[idx], substation)

        best_sol_cluster = nothing
        best_cal_cluster = nothing
        best_cost_cluster = Inf
        worst_cost_cluster = 0
        best_node_cluster  = nothing 

        cluster_starting_turbines = [i for i in starting_turbines if i in clusters[idx]]
        for start_turbine in cluster_starting_turbines

            sol_cluster, cal_cluster =  perform_sweep_on_cluster(data, cluster_polar_angles, substation, start_turbine, clockwise, string)
            sol_cluster.cost = improvementHeuristic.non_penalised_cost_function(sol_cluster, data) # DELETE THIS IF PENALISED

            if sol_cluster.cost < best_cost_cluster
                best_cost_cluster = sol_cluster.cost
                best_sol_cluster = sol_cluster
                best_cal_cluster = cal_cluster 
                best_node_cluster = start_turbine
            end 
    
            if sol_cluster.cost > worst_cost_cluster
                worst_cost_cluster = sol_cluster.cost
            end 

        end

        push!(best_sol, best_sol_cluster)

        #println("For substation ", substation, " the best solution starts at node ", best_node_cluster, ". It has a cost of ", best_cost_cluster, " which is ", worst_cost_cluster - best_cost_cluster, " than the worst starting point."  )

    end 

    sol_combined = minimumSpanningTree.solution()
    for sol in best_sol
        sol_combined.adjacency_list = merge!(append_vectors, sol_combined.adjacency_list, sol.adjacency_list)
    end 

    final_sol, final_cal =  pressureDropCalculations.run(sol_combined, data, substations)
    #final_sol.cost = improvementHeuristic.cost_function(final_sol, data, substations) # Comment out code if PENALISED
    final_sol.cost = improvementHeuristic.non_penalised_cost_function(final_sol, data) # Comment in code if not penalised 

    return final_sol, final_cal 
    
end 


function perform_sweep_on_cluster_implementation(data, polar_angles, substation, starting_turbine, clockwise = true, strings = false, balancing = false, grupper = [])
    ordered_polar_angles = sweep.order_polar_angles(polar_angles, starting_turbine, clockwise)
    node_list = [substation]
    sol_list = Set()
    previous_sol = nothing
    previous_cal = nothing
    num = 1
    grup_num = 1

    while num <= length(ordered_polar_angles)
        #start_time = time_ns()
        append!(node_list, ordered_polar_angles[num])

        if strings
            sol = nearestNeighbor.run(data, node_list, substation)
        else 
            g = Graph.create(data, node_list)
            sol = minimumSpanningTree.run(g, substation)
            sol.adjacency_list = Graph.convert_to_directed_graph(substation, sol.adjacency_list)
        
        end
        sol, cal  = pressureDropCalculations.run(sol, data, substation)

        # Feasibility check 
        if cal.feasible
            previous_sol = deepcopy(sol)
            previous_cal = deepcopy(cal)

            row = [substation, starting_turbine, grup_num, ordered_polar_angles[num]]
            #print(row)
            push!(grupper, row)

            num += 1

        else

            push!(sol_list, previous_sol)
            node_list = [substation]
            grup_num += 1

        end 

    end 

    push!(sol_list, previous_sol)
    
    # Merge sols from sol_list
    sol_cluster = minimumSpanningTree.solution()
    for sub_sol in sol_list
        sol_cluster.adjacency_list = merge!(append_vectors, sol_cluster.adjacency_list, sub_sol.adjacency_list)
    end 

    sol_final_cluster, cal_final_cluster =  pressureDropCalculations.run(sol_cluster, data, substation)
    sol_final_cluster.cost, sol_final_cluster.digging_cost = improvementHeuristic.cost_function(sol_final_cluster, data, [substation], cal_final_cluster.feasible, balancing)


    return sol_final_cluster, cal_final_cluster, grupper

end

function run_sweep_implementation(data, coords, substations, starting_turbines, clockwise = true, string = false, balancing = false)
    # Assign turbines to a substation 
    centroids, clusters = Kmeans.run(data, coords, substations) # Currently ssuming fixed substations

    best_sol = []
    worst_sol = []
    grupper = []
   
    # for each substation find the best solution from different starting points 
    for(idx, substation) in enumerate(substations)
        cluster_polar_angles = sweep.calculate_polar_angles(coords, clusters[idx], substation)

        best_sol_cluster = nothing
        best_cal_cluster = nothing
        best_cost_cluster = Inf
        best_node_cluster  = nothing 

        worst_sol_cluster = nothing
        worst_cal_cluster = nothing
        worst_cost_cluster = 0
        worst_node_cluster  = nothing 

        cluster_starting_turbines = [i for i in starting_turbines if i in clusters[idx]]
        for start_turbine in cluster_starting_turbines

            if start_turbine == 25
                print(start_turbine)
            end 

            

            sol_cluster, cal_cluster, grupper =  perform_sweep_on_cluster_implementation(data, cluster_polar_angles, substation, start_turbine, clockwise, string, balancing, grupper)

            if sol_cluster.cost < best_cost_cluster
                best_cost_cluster = deepcopy(sol_cluster.cost)
                best_sol_cluster = deepcopy(sol_cluster)
                best_cal_cluster = deepcopy(cal_cluster) 
                best_node_cluster = deepcopy(start_turbine)
            end 
    
            if sol_cluster.cost > worst_cost_cluster
                worst_cost_cluster = deepcopy(sol_cluster.cost)
                worst_sol_cluster = deepcopy(sol_cluster)
                worst_cal_cluster = deepcopy(cal_cluster) 
                worst_node_cluster = deepcopy(start_turbine)
            end 

        end

        println("Substation ", substation, " worst starting turbine ",worst_node_cluster)

        push!(best_sol, best_sol_cluster)
        push!(worst_sol, worst_sol_cluster)

        #println("For substation ", substation, " the best solution starts at node ", best_node_cluster, ". It has a cost of ", best_cost_cluster, " which is ", worst_cost_cluster - best_cost_cluster, " than the worst starting point."  )

    end 

    best_sol_combined = minimumSpanningTree.solution()
    for sol in best_sol
        best_sol_combined.adjacency_list = merge!(append_vectors, best_sol_combined.adjacency_list, sol.adjacency_list)
    end

    best_final_sol, best_final_cal =  pressureDropCalculations.run(best_sol_combined, data, substations)
    best_final_sol.cost, best_final_sol.digging_cost  = improvementHeuristic.cost_function(best_final_sol, data, substations, best_final_cal.feasible, balancing)
    
    worst_sol_combined = minimumSpanningTree.solution()
    for sol in worst_sol
        worst_sol_combined.adjacency_list = merge!(append_vectors, worst_sol_combined.adjacency_list, sol.adjacency_list)
    end

    worst_final_sol, worst_final_cal =  pressureDropCalculations.run(worst_sol_combined, data, substations)
    worst_final_sol.cost, worst_final_sol.digging_cost  = improvementHeuristic.cost_function(worst_final_sol, data, substations, worst_final_cal.feasible, balancing)

    return best_final_sol, best_final_cal, worst_final_sol, worst_final_cal, grupper  
    
end 

end 



