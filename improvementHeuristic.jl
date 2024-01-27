module improvementHeuristic
include("pressureDropCalculations.jl")
include("writeResults.jl")

function find_edge_from_index(data,index)
    i = (index-1) % (data.n) + 1
    j = floor(Int64, (index-1)/data.n) + 1
    return i, j
end

function find_index_from_edge(data, i, j)
    index = (j - 1) * data.n + i
    return index
end

function one_opt_move(sol, cal, data, edge, i, j, new_i)
    
    # DELETE FROM SOLUTION
    net = deepcopy(sol)
    
    #delete old edge
    delete!(net.edges, edge)
    delete!(net.pipetype, edge)
    if length(net.adjacency_list[i]) > 1
        deleteat!(net.adjacency_list[i], findall(x->x==j,net.adjacency_list[i]))
    else
        delete!(net.adjacency_list, i)
    end 

    #add new edge
    if haskey(net.adjacency_list, new_i)
        push!(net.adjacency_list[new_i], j)
    else
        net.adjacency_list[new_i] = [j]
    end 

    new_edge = find_index_from_edge(data, new_i, j)
    push!(net.edges, new_edge)

    # DELETE FROM CALCULATION
    calc = deepcopy(cal)
    delete!(calc.velocity, edge)
    delete!(calc.pressuredrop, i)
    delete!(calc.pressuredrop, j)

    return net, calc
end

function find_affected_edges(sol, data, j, i, substations)
    affected_edges = Set{Int64}()

    # Find path from substations to j
    function find_path_to_node(sol, current_node, target_node, visited)
        if current_node == target_node
            return true
        end
        push!(visited, current_node)

        for neighbor in get(sol.adjacency_list, current_node, [])
            if neighbor ∉ visited
                edge_index = find_index_from_edge(data, current_node, neighbor)
                if find_path_to_node(sol, neighbor, target_node, visited)
                    push!(affected_edges, edge_index)
                    return true
                end
            end
        end
        return false
    end

    for substation in substations
        find_path_to_node(sol, substation, j, Set{Int64}())
        find_path_to_node(sol, substation, i, Set{Int64}())
    end

    # Find paths from j to leaves
    function dfs_leaves(sol, node, visited)
        if node in visited
            return
        end
        push!(visited, node)

        if isempty(get(sol.adjacency_list, node, []))  # Node is a leaf
            return
        end

        for neighbor in sol.adjacency_list[node]
            edge_index = find_index_from_edge(data, node, neighbor)
            push!(affected_edges, edge_index)
            dfs_leaves(sol, neighbor, visited)
        end
    end

    nodes_connected_to_substations = [find_edge_from_index(data, i)[2] for i in affected_edges if find_edge_from_index(data, i)[1] in substations ]
    for node in nodes_connected_to_substations 
        dfs_leaves(sol, node, Set{Int64}())
    end 

    return affected_edges
end


#= function evaluate_move(sol, sol_move, data, substations, allow_illegal, balancing)
    improved = false
    best_cal = nothing
    best_sol = deepcopy(sol)
    sol_move, cal_move = pressureDropCalculations.run(sol_move, data, substations)
    if cal_move.feasible || allow_illegal
        sol_move.cost,  sol_move.digging_cost = cost_function(sol_move, data, substations, cal_move.feasible, balancing)
        if sol_move.cost < best_sol.cost
            #println("Improved")
            improved = true
        end 
        return improved, sol_move, cal_move 
    end 

    return improved, best_sol, best_cal 

end  =#

function fast_evaluate_move(sol, sol_move, cal_move, data, substations, affected_edges, allow_illegal, balancing)
    improved = false
    best_cal = nothing
    best_sol = deepcopy(sol)
    sol_move, cal_move = pressureDropCalculations.run_fast(sol_move, cal_move, data, substations, affected_edges, allow_illegal)
    if cal_move.feasible || allow_illegal
        sol_move.cost,  sol_move.digging_cost = cost_function(sol_move, data, substations, cal_move.feasible, balancing)
        if sol_move.cost < best_sol.cost
            #println("Improved")
            improved = true
        end 
        return improved, sol_move, cal_move 
    end 

    return improved, best_sol, best_cal 

end 

function find_descendants_of_j(sol, j)
    descendants = Set(j)

    if haskey(sol.adjacency_list, j) == false
        return descendants
    end 

    neighbor_list =  deepcopy(sol.adjacency_list[j])

    it = 1
    while (it <= length(neighbor_list))
        neighbor = neighbor_list[it]
        push!(descendants, neighbor)

        if haskey(sol.adjacency_list, neighbor)
            next_neighbor = deepcopy(sol.adjacency_list[neighbor]) 
            append!(neighbor_list, next_neighbor)
        end 

        it += 1

    end 
    return descendants
end

function find_non_descendants_of_j(data, sol, j, substations)    
    descendants = find_descendants_of_j(sol, j)
    #non_descendants_old = [i for i in 1:data.n if i ∉ descendants] # includes the substations
    non_descendants = []
    for i in 1:data.n
        if ((i ∉ descendants) && data.distance[i,j] <= 5) || (i in substations)
            append!(non_descendants, i)
        end

    end 

    #println("Antal non descendants sparet ", length(non_descendants_old) - length(non_descendants))
    return non_descendants
end

function find_potential_leafs(data, sol, j, substations)
    
    leafs = Set()
    decendants = find_descendants_of_j(sol, j)
    
    for i in 1:data.n
        if (haskey(sol.adjacency_list, i) == false) && (i ∉ decendants) 
            push!(leafs, i)

        elseif i in substations # includes the substations
            push!(leafs, i)
        end 
    end 

    return leafs 
end

function first_improvement(sol, cal, data, substations, string, allow_illegal, balancing, start_time = time_ns(), timelimit = 60)
    its = 0
    for edge in collect(sol.edges)

        improved = false
        i, j = find_edge_from_index(data, edge)

        if string
            possible_i = find_potential_leafs(data, sol, j, substations)
        else
            possible_i = find_non_descendants_of_j(data, sol, j, substations)
        end 
        
        for new_i in possible_i
            # Make one opt move and evaluate solution
            sol_move, cal_move = one_opt_move(sol, cal, data, edge, i, j , new_i)
            affected_edges = find_affected_edges(sol_move, data, j, i, substations)
            improved, sol_move, cal_move = fast_evaluate_move(sol, sol_move, cal_move, data, substations, affected_edges, allow_illegal, balancing)
            
            #OLD METHOD
            #improved, sol_move, cal_move = evaluate_move(sol, sol_move, data, substations, allow_illegal)
    
            if improved || (elapsed_time(start_time) > timelimit)
                return sol_move, cal_move
            end
        end

        its += 1
    end
    return sol, cal
end

function best_improvement(sol, cal, data, substations, string, allow_illegal, balancing, start_time = time_ns(), timelimit = 60 )
    best_sol = deepcopy(sol)
    best_cost_improvement = 0.0 
    best_cal = deepcopy(cal)
    its = 0
    for edge in collect(sol.edges)
        i, j = find_edge_from_index(data, edge)

        if string
            possible_i = find_potential_leafs(data, sol, j, substations)
        else
            possible_i = find_non_descendants_of_j(data, sol, j, substations)
        end

        for new_i in possible_i
                
            sol_move, cal_move = one_opt_move(sol, cal, data, edge, i, j, new_i)
            affected_edges = find_affected_edges(sol_move, data, j, i, substations)
            improved, sol_move, cal_move = fast_evaluate_move(sol, sol_move, cal_move, data, substations, affected_edges, allow_illegal, balancing)

            #println("BEST IMPROVEMENT: from ", (i,j), "to ", (new_i, j))
            
            #OLD METHOD
            #improved, sol_move, cal_move = evaluate_move(sol, sol_move, data, substations, allow_illegal)
            
            cost_improvement = sol.cost - sol_move.cost
            
            if improved && cost_improvement > best_cost_improvement
                best_sol = deepcopy(sol_move)
                best_cost_improvement = cost_improvement
                best_cal = deepcopy(cal_move)
            end


            if (elapsed_time(start_time) > timelimit)     
                return best_sol, best_cost_improvement, best_cal
            end 

        end

        its += 1

    end

    return best_sol, best_cost_improvement, best_cal
end

function cost_function(sol, data, substations, feasible = true, balancing = false )
    ###### Cost of current chosen pipes   
    total_cost = 0.0
    sol.t_connections = 0
    sol.m_connections = 0
    sol.crossings = 0 
    sol.digging_cost = 0
    total_digging_cost =0
    for (num, arc) in enumerate(keys(sol.pipetype))
        i, j = find_edge_from_index(data, arc)
        # Find the index of the assigned diameter
        diameter_index = sol.pipetype[arc]
        # Calculate cost for this arc
        arc_digging_cost = data.digging_cost * data.distance[i, j]*1000 # from [km] to [m]
        arc_cost = data.pipe_costs[diameter_index] * data.distance[i, j]*1000 # from [km] to [m]
        total_cost += (arc_cost + arc_digging_cost)
        total_digging_cost += arc_digging_cost
        
        ###### Cost of crossing edges
        edge_length = length(keys(sol.pipetype))
        if (edge_length > 1) && (num+1<= edge_length)
            for oarc in collect(keys(sol.pipetype))[(num+1):edge_length]

                i_oarc, j_oarc = find_edge_from_index(data, oarc)

                if (j ≠ i_oarc) && (i ≠ j_oarc) && (i ≠ i_oarc)
                    # Edge (p1, p2)
                    p1 = data.coords[i,:]
                    p2 = data.coords[j,:]
                    # Edge (p3, p4)
                    p3 = data.coords[i_oarc,:]
                    p4 = data.coords[j_oarc,:]
                    
                    result = do_line_segments_intersect(p1, p2, p3, p4)
                    
                    sol.crossings += result * 1
                    total_cost += data.bigM * result 

                    #if result == true
                    #    println("Edge between ", i, " and ", j, " intersects with edge ", i_oarc, " and ", j_oarc)
                    #end 
                end 
            end
        end 
    end

    ##### Cost due to max connection at manifold and turbine 
    for node in keys(sol.adjacency_list)
        
        # Calculatue connections to manifold
        if node in substations
            m_connections = length(sol.adjacency_list[node])

            if m_connections > data.max_m_connections
                total_cost += data.bigM * (m_connections - data.max_m_connections)
                sol.m_connections +=  (m_connections - data.max_m_connections)
            end 
        
        # Calculate connections to turbines
        else
            t_connections = length(sol.adjacency_list[node])

            if t_connections > data.max_t_connections
                total_cost += data.bigM * (t_connections - data.max_t_connections)
                sol.t_connections += t_connections - data.max_t_connections
            end
        end 
    end 

    ##### Calculate cost of non-balanced strings
    if balancing
        for substation in substations
            m_connections = sol.adjacency_list[substation]
            for connection in m_connections
                sol.balanced[connection] = length(find_descendants_of_j(sol, connection)) + 1  # comment out code in experiment with balanced strings
                #sol.balanced[connection] = 0 # Comment code in experiment with balanced strings 
            end 
            
            total_cost += (maximum(values(sol.balanced)) - minimum(values(sol.balanced))) * 10e5 # comment out code in experiemtn with balanced strings 
            
        end
    else
        #=for substation in substations
            m_connections = sol.adjacency_list[substation]
            for connection in m_connections
                sol.balanced[connection] = 0 # Comment code in experiment with balanced strings 
            end 
        end=#
        for substation in substations
            #if haskey(sol.adjacency_list, substation)
            m_connections = sol.adjacency_list[substation]
            for connection in m_connections
                sol.balanced[connection] = 0
                end
            #end
        end

    end  

    # Add extra cost for illegal solution
    if !feasible
        total_cost += data.bigM * 100
    end 

    return total_cost, total_digging_cost 
end


function check_crossings(sol, data, instance, zone)
    results = []
    for (num, arc) in enumerate(keys(sol.pipetype))
        i, j = find_edge_from_index(data, arc)
        edge_length = length(keys(sol.pipetype))
        if (edge_length > 1) && (num+1<= edge_length)
            for oarc in collect(keys(sol.pipetype))[(num+1):edge_length]

                i_oarc, j_oarc = find_edge_from_index(data, oarc)

                if (j ≠ i_oarc) && (i ≠ j_oarc) && (i ≠ i_oarc)
                    # Edge (p1, p2)
                    p1 = data.coords[i,:]
                    p2 = data.coords[j,:]
                    # Edge (p3, p4)
                    p3 = data.coords[i_oarc,:]
                    p4 = data.coords[j_oarc,:]
                    
                    result = do_line_segments_intersect(p1, p2, p3, p4)

                    if result == true
                        its_row = [instance, zone, i, j, i_oarc, j_oarc]
                        push!(results, its_row)
                        println("Edge between ", i, " and ", j, " intersects with edge ", i_oarc, " and ", j_oarc)
                    end 
                end 
            end
        end 
    end
    return results
end 

function non_penalised_cost_function(sol, data)
    ###### Cost of current chosen pipes   
    total_cost = 0.0
    for (num, arc) in enumerate(keys(sol.pipetype))
        i, j = find_edge_from_index(data, arc)
        # Find the index of the assigned diameter
        diameter_index = sol.pipetype[arc]
        # Calculate cost for this arc
        arc_cost = data.pipe_costs[diameter_index] * data.distance[i, j]*1000 # from [km] to [m]
        total_cost += arc_cost
        
    end
    return total_cost 
end

function cost_function_bounds(sol, data, bound)
    # Calculate cost of current chosen pipes   
    total_cost = 0.0
    for (num, arc) in enumerate(keys(sol.pipetype))
        i, j = find_edge_from_index(data, arc)
        
        # Find the index of min or max pipe size
        if bound == "min"
            diameter_index = 1
        elseif bound == "max"
            diameter_index = length(data.pipetype_inner_diameters)
        end 
        
        # Calculate cost for this arc
        arc_cost = data.pipe_costs[diameter_index] * data.distance[i, j]*1000 # from [km] to [m]
        total_cost += arc_cost
        
    end

    return total_cost
end

function cross_product(p1, p2)
    return p1[1] * p2[2] - p1[2] * p2[1]
end


function do_line_segments_intersect(p1, p2, p3, p4; tol=1e-8)
    A = (p1[1], p1[2])
    B = (p2[1], p2[2])
    C = (p3[1], p3[2])
    D = (p4[1], p4[2])

    # Check if the line segments share an endpoint
    if A == C || A == D || B == C || B == D
        return false
    end

    line1 = (A, B)
    line2 = (C, D)

    function det(a, b)
        return a[1] * b[2] - a[2] * b[1]
    end

    xdiff = (line1[1][1] - line1[2][1], line2[1][1] - line2[2][1])
    ydiff = (line1[1][2] - line1[2][2], line2[1][2] - line2[2][2])

    div = det(xdiff, ydiff)
    if div == 0
        return false  # lines do not intersect
    end

    d = (det(line1[1:2]...), det(line2[1:2]...))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div

    # Check if the intersection point is within the line segments with tolerance
    if min(line1[1][1], line1[2][1]) - tol ≤ x ≤ max(line1[1][1], line1[2][1]) + tol &&
        min(line1[1][2], line1[2][2]) - tol ≤ y ≤ max(line1[1][2], line1[2][2]) + tol &&
        min(line2[1][1], line2[2][1]) - tol ≤ x ≤ max(line2[1][1], line2[2][1]) + tol &&
        min(line2[1][2], line2[2][2]) - tol ≤ y ≤ max(line2[1][2], line2[2][2]) + tol
        return true
    else
        return false
    end
end


function do_line_segments_intersect_old(p1, p2, p3, p4)
    v1 = p2 - p1
    v2 = p4 - p3
    cross1 = cross_product(p3 - p1, v1)
    cross2 = cross_product(p4 - p1, v1)
    cross3 = cross_product(p1 - p3, v2)
    cross4 = cross_product(p2 - p3, v2)

    if cross1 * cross2 ≤ 0 && cross3 * cross4 ≤ 0
        return true
    end
    return false
end

function initalize_result_table(start_time, sol, metaheuristic = "SA")

    if metaheuristic == "SA"
        #[["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections",  "BalanceDiff", "Digging Cost", 'Temperature', 'Alpha','Reheats]]
        results = [[0, "Sweep", elapsed_time(start_time), sol.cost, sol.cost, sol.crossings, sol.m_connections, sol.t_connections, (maximum(values(sol.balanced)) - minimum(values(sol.balanced))), sol.digging_cost ,0, 0, 0]]

    else 
        #[["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections", "BalanceDiff", "Digging Cost", "K", "First Improvement", "Feasible"]]
        results = [[0, "Sweep", elapsed_time(start_time), sol.cost, sol.cost, sol.crossings, sol.m_connections, sol.t_connections, (maximum(values(sol.balanced)) - minimum(values(sol.balanced))), sol.digging_cost, 0, 0, true]]

    end 
    return results 
end 


function initialise_temperature_and_alpha(data, sol, cal, K, string, substations, start_time, timelimit, reheat, balancing)

    tuning_starttime = time_ns()
    gamma = 0.9
    beta_m = 0
    beta_p = 0
    its = 0
    best_sol = deepcopy(sol)
    best_cal = deepcopy(cal)
    delta_total = 0
    allow_illegal = false

    while (its < K) && (elapsed_time(start_time) < timelimit)

        # Pick edge to remove
        edge_remove = rand(collect(sol.edges))
        i, j = find_edge_from_index(data, edge_remove) 
        
        if string
            possible_i = find_potential_leafs(data, sol, j, substations)
        else
            possible_i = find_non_descendants_of_j(data, sol, j, substations)
        end

        # Find new edge to add 
        new_i_idx = rand(1:length(possible_i))
        new_i = collect(possible_i)[new_i_idx]

        # Evaluate solution
        #net = one_opt_move(sol, cal, data, edge_remove, i, j , new_i)
        #_, sol_move, cal_move = evaluate_move(sol, net, data, substations)

        sol_move, cal_move = one_opt_move(sol, cal, data, edge_remove, i, j, new_i)
        affected_edges = find_affected_edges(sol_move, data, j, i, substations)
        _, sol_move, cal_move = fast_evaluate_move(sol, sol_move, cal_move, data, substations, affected_edges, allow_illegal, balancing)

        delta_non_penalised = non_penalised_cost_function(sol_move, data) - sol.cost 
        delta = sol_move.cost - sol.cost
        
        if (delta < 0 || rand() < gamma) && (cal_move ≠ nothing) 

            sol = deepcopy(sol_move)
            cal = deepcopy(cal_move)

            if delta < 0 
                beta_p += 1

                if sol.cost < best_sol.cost 
                    best_sol = deepcopy(sol)
                    best_cal = deepcopy(cal)
                end 

            else
                beta_m += 1
                #delta_total += delta
                delta_total += max(0, delta_non_penalised)
            end 

        else
            beta_m += 1
            #delta_total += delta
            delta_total += max(0, delta_non_penalised)
        end 
        its += 1
    end 

    eta = 0.0000001

    if reheat == Inf
        expected_its = timelimit/elapsed_time(tuning_starttime) * K 
        Tstart = delta_total/(log(beta_m/(beta_m*gamma-beta_p*(1-gamma))))
        alpha = exp(log(-1/(Tstart*log(eta)))/expected_its)
    
    else
        expected_its = reheat
        Tstart = delta_total/(log(beta_m/(beta_m*gamma-beta_p*(1-gamma))))
        alpha = exp(log(-1/(Tstart*log(eta)))/expected_its)

    end 

    #println("SA --  Its ", its, " Temperature ", Tstart, " alpha ", alpha)

    return Tstart, alpha, best_sol, best_cal 

end

function simulated_annealing(sol, cal, data, substations, string, start_time, timelimit, reheat, balancing)
    
    results = initalize_result_table(start_time, sol, "SA")
    best_sol = deepcopy(sol)
    best_cal = deepcopy(cal)


    Tstart, alpha, best_sol, best_cal = initialise_temperature_and_alpha(data, sol, cal, 10000, string, substations, start_time, timelimit, reheat, balancing)    
    its_row = [0, "SA", elapsed_time(start_time), best_sol.cost, sol.cost, sol.crossings, sol.m_connections, sol.t_connections, (maximum(values(sol.balanced)) - minimum(values(sol.balanced))), sol.digging_cost, Tstart, alpha, 0 ]
    push!(results, its_row)

    t = Tstart 
    its = 1
    it_not_improved = 0
    antal_reheats = 0
    allow_illegal = false
    while elapsed_time(start_time) < timelimit
        # Pick edge to remove
        edge_remove = rand(collect(sol.edges))
        i, j = find_edge_from_index(data,edge_remove) 
        
        if string
            possible_i = find_potential_leafs(data, sol, j, substations)
        else
            possible_i = find_non_descendants_of_j(data, sol, j, substations)
        end

        # Find new edge to add 
        new_i_idx = rand(1:length(possible_i))
        new_i = collect(possible_i)[new_i_idx]

        # Evaluate solution
        #net = one_opt_move(sol, cal, data, edge_remove, i, j , new_i)
        #_, sol_move, cal_move = evaluate_move(sol, net, data, substations)
        sol_move, cal_move = one_opt_move(sol, cal, data, edge_remove, i, j, new_i)
        affected_edges = find_affected_edges(sol_move, data, j, i, substations)
        _, sol_move, cal_move = fast_evaluate_move(sol, sol_move, cal_move, data, substations, affected_edges, allow_illegal, balancing)
        
        delta = sol_move.cost - sol.cost
        # Pick solution if better or with a probability
        if (delta < 0 || rand() < exp(-delta/t)) && (cal_move ≠ nothing) 
            sol = deepcopy(sol_move)
            cal = deepcopy(cal_move)

            if sol.cost < best_sol.cost
                #println("SA: Improved from ", best_sol.cost, " to ", sol.cost)
                best_sol = deepcopy(sol)
                best_cal = deepcopy(cal)
            else
                #println("SA: Picked a worse solution ", sol.cost)
                it_not_improved += 1
            end

        else
            it_not_improved += 1
        end 

        if (its % 10 == 0) 
            #["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections", "BalanceDiff", "Temperature", "Alpha","Reheats"]
            its_row = [its, "SA", elapsed_time(start_time), best_sol.cost, sol.cost, sol.crossings, sol.m_connections, sol.t_connections,  (maximum(values(sol.balanced)) - minimum(values(sol.balanced))), sol.digging_cost, t, alpha, antal_reheats ]
            push!(results, its_row)
        end 

        t = t * alpha 
        its += 1 

        if it_not_improved > reheat # husk at lave det om til iterationer uden forbedringer
            t = Tstart
            it_not_improved = 0
            antal_reheats += 1
        end
    end 
    #println("DONE ", elapsed_time(start_time) )
    println("SA its ", its)

    #_, test_cal = pressureDropCalculations.run(best_sol, data, substations)
    #println("Velocity ens: ", test_cal.velocity == best_cal.velocity)
    #println("Pressuredrop ens: ", test_cal.velocity == best_cal.velocity)

    return best_sol, best_cal, results, Tstart, alpha
end 

function shaking_phase(sol, cal, k, data, substations, string, allow_illegal, balancing, start_time = time_ns(), timelimit = 60)
    sol_prime = deepcopy(sol)
    cal_prime = deepcopy(cal)
    count_k = 0
    its = 0
    #for _ in 1:k
    while (count_k < k) && (elapsed_time(start_time) < timelimit) # Makes sure we make k alternations
        edge_to_remove = rand(collect(sol_prime.edges))
        i, j = find_edge_from_index(data, edge_to_remove)

        if string
            possible_i = find_potential_leafs(data, sol_prime, j, substations)
        else
            possible_i = find_non_descendants_of_j(data, sol_prime, j, substations)
        end
        
        # Find new edge to add 
        new_i_idx = rand(1:length(possible_i))
        new_i = collect(possible_i)[new_i_idx]
        
        sol_net, cal_net = one_opt_move(sol_prime, cal_prime, data, edge_to_remove, i, j, new_i)
        affected_edges = find_affected_edges(sol_net, data, j, i, substations)
        sol_net, cal_net = pressureDropCalculations.run_fast(sol_net, cal_net, data, substations, affected_edges, allow_illegal)

        if allow_illegal || cal_net.feasible
            sol_net.cost, sol_net.digging_cost = cost_function(sol_net, data, substations, cal_net.feasible, balancing )
            sol_prime = deepcopy(sol_net)
            cal_prime = deepcopy(cal_net)
            count_k += 1

            if cal_net.feasible == false
                #print("Shaking not feasible")
            end 
            #println("SHAKING: from ", (i,j), "to ", (new_i, j))
        end
        its += 1

    end

    return sol_prime, cal_prime
end

function local_search(sol_prime, cal_prime, data, substations, string, allow_illegal, start_time, first_improve, balancing, timelimit, results, k)
    best_sol = deepcopy(sol_prime)
    best_cal = deepcopy(cal_prime)
    improved = true
    its = 0
    while improved && (elapsed_time(start_time) < timelimit)

        if first_improve     
            sol_double_prime, cal_double_prime = first_improvement(best_sol, best_cal, data, substations, string, allow_illegal, balancing, start_time, timelimit)   
        else   
            sol_double_prime, _, cal_double_prime = best_improvement(best_sol, best_cal, data, substations, string, allow_illegal, balancing, start_time, timelimit)    
        end 

        if sol_double_prime.cost < best_sol.cost
            best_sol = deepcopy(sol_double_prime)
            best_cal = deepcopy(cal_double_prime)
            improved = true
            #println("Improved to ", best_sol.cost)
            its_row = [its, "LocalSearch", elapsed_time(start_time), best_sol.cost, sol_double_prime.cost, sol_double_prime.crossings, sol_double_prime.m_connections, sol_double_prime.t_connections,(maximum(values(sol_double_prime.balanced)) - minimum(values(sol_double_prime.balanced))), k, first_improve, cal_double_prime.feasible]
            push!(results, its_row)
        else
            improved = false
        end 

        its += 1
        #println(its)

        
    end
    #println("Iterations local search ", its)
    return best_sol, best_cal, results
end 

function VNS(sol, cal, data, substations, start_time, timelimit, Kmax=4, first_improvement = true, string = true, allow_illegal = true, balancing = false)
    results = initalize_result_table(start_time, sol, "VNS")
    best_sol = deepcopy(sol)
    best_cal = deepcopy(cal)
    k = 1    
    its = 0
    #allow_illegal = true
    while elapsed_time(start_time) < timelimit
        # Shaking
        sol_prime, cal_prime = shaking_phase(best_sol, best_cal, k, data, substations, string, allow_illegal, balancing, start_time, timelimit)
        # Local Search 
        sol_double_prime, cal_double_prime, results = local_search(sol_prime, cal_prime, data, substations, string, allow_illegal, start_time, first_improvement, balancing, timelimit, results, k)

        # Move or Not
        if sol_double_prime.cost < best_sol.cost
            best_sol = deepcopy(sol_double_prime)
            best_cal = deepcopy(cal_double_prime)
            k = 1
        else
            k += 1
            if k > Kmax
                k = 1
            end
        end

        its += 1
        
        #["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections","BalanceDiff", 'K', "First Improvement", "Feasible"]
        its_row = [its, "VNS",elapsed_time(start_time), best_sol.cost, sol_double_prime.cost, sol_double_prime.crossings, sol_double_prime.m_connections, sol_double_prime.t_connections,(maximum(values(sol_double_prime.balanced)) - minimum(values(sol_double_prime.balanced))), sol_double_prime.digging_cost, k, first_improvement, cal_double_prime.feasible]
        push!(results, its_row)

        t_it = elapsed_time(start_time)
        #println("Time new it ", t_it-t_start)
    
    end
    println("Antal iterationer ", its)

    #_, test_cal = pressureDropCalculations.run(best_sol, data, substations)
    #println("Velocity ens: ", test_cal.velocity == best_cal.velocity)
    #println("Pressuredrop ens: ", test_cal.velocity == best_cal.velocity)

    return best_sol, best_cal, results
end

function elapsed_time(start_time)
    return round((time_ns()-start_time)/1e9,digits = 3)
end 

end