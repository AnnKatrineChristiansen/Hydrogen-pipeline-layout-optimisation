using JuMP
using CoolProp
using Gurobi

include("readData.jl")
include("graph.jl")
include("pressureDropCalculations.jl")
include("writeResults.jl")
include("improvementHeuristic.jl")

mutable struct Solution
    nodes::Set{Int64}
    adjacency_list::Dict{Int64, Vector{Int64}}
    edges::Set{Int64}
    pipetype::Dict{Int64, Int64}
    flowrate::Dict{Int64, Int64}
    cost::Float64
    crossings::Int64
    m_connections::Int64 #Dict{Int64, Int64}
    t_connections::Int64
    balanced::Dict{Int64, Int64}
    digging_cost::Float64
    Solution() = new(Set{Int64}(), Dict{Int64, Vector{Int64}}(), Set{Int64}(), Dict{Int64, Int64}(), Dict{Int64, Int64}(), 0.0, 0, 0, 0, Dict(),0.0)
end 

function run_mip(instance, zone, number_of_manifolds, pressure, strings, timelimit, folder_name, crossings, it = 0)

    if isempty(crossings)
        println("RUN MIP ON INSTANCE ", instance, " ZONE ", zone)

    else
        println("OPTMISING MIP WITHOUT CROSSINGS INSTANCE ", instance, " ZONE ", zone)
    end 

    if isfile("results/"*folder_name*"/crossings.txt")
        all_crossings = WriteResults.read_txt_file(folder_name*"/crossings.txt")
        crossings = [row for row in all_crossings if row[1] == instance && row[2] == parse(Int64,zone)]
        println("Found crossings", crossings)
    end 

    # Load data
    filename = "Data/instances/"*instance*"/layout.txt"
    coords, data = readData.run(filename, parse(Int, zone))
    ########### Data
    N = data.n
    S = number_of_manifolds # number of substations
    T = length(data.pipe_costs) # number of pipetypes
    P = 363 * 24#/1000 # mass flow production for each turbine
    dist = data.distance # dist matrix
    u = data.pipe_costs*1000 # pipe costs / km
    p_start = pressure
    density = CoolProp.PropsSI("D", "T", 3 + 273.15 , "P", p_start * 1e5, "H2")

    C = ceil(((N-S)/S)/9)
    N_m = ceil((N)/S)
    balanced_u = P * ceil((N_m-1)/C)
    balanced_l = P * floor((N_m-1)/C)

    manifold_array = collect(1:number_of_manifolds)

    #max mass flow rates to keep velocity below 25 m/S
    # V = Q/A  => V * A 
    k = zeros(Float64, length(data.pipetype_inner_diameters))
    for (idx_d, d) in enumerate(data.pipetype_inner_diameters)
        k[idx_d] = 25 * (π * (d / 2 / 1000)^2) * (24 * 3600) * density
    end 

    ########### Model
    model = Model(Gurobi.Optimizer);

    # set timelimit
    set_time_limit_sec(model, timelimit)

    # 1 if arc i, j is of pipetype t, otherwise 0
    @variable(model, x[1:N,1:N,1:T], Bin); 
    # flow in arc i,j
    @variable(model, f[1:N, 1:N] >= 0, Int); 
    # 1 if an arc is built
    @variable(model, y[1:N, 1:N], Bin);
    # 
    @variable(model, c[1:S], Bin);
    # aux
    @variable(model, aux, Int);


    # Sum costs of chosen pipes
    @objective(model, Min, sum(u[t] * dist[i,j] * x[i,j,t] for i=1:N, j=1:N, t=1:T));

    # Only one pipe type can be selected for a built arc
    @constraint(model,[i = 1:N, j = 1:N], sum(x[i,j,t] for t=1:T) == y[i,j]);

    # The flow exiting each node h is equal to the energy entering h plus the production of that node
    @constraint(model,[h = S+1:N],  sum(i == h ? 0 : f[h,i] - f[i,h] for i=1:N) == P);

    # Instead ensure that the flow does not exceed the max velocity of the installed pipe
    @constraint(model,[i = 1:N, j = 1:N], sum(k[t] * x[i,j,t] for t=1:T) >= f[i,j]);

    # Ensure that exactly one arc leaves a turbine
    @constraint(model,[h = S+1:N], sum(j == h ? 0 : y[h,j] for j=1:N) == 1) 

    # Nothing exits a manifold
    @constraint(model,[h = 1:S], sum(j == h ? 0 : y[h,j] for j=1:N) == 0)

    # At most 3 enters a turbine
    if strings
        @constraint(model,[h = S+1:N], sum(i == h ? 0 : y[i,h] for i=1:N) <= 1)
    else
        @constraint(model,[h = S+1:N], sum(i == h ? 0 : y[i,h] for i=1:N) <= 3)
    end
    
    # At exactly C pipes enters a manifold
    @constraint(model,[h = 1:S], sum(i == h ? 0 : y[i,h] for i=1:N) == C)

    # Balanced constraint
    @constraint(model, [h=1:S, i=1:N], f[i, h] <= balanced_u)
    @constraint(model, [h=1:S, i=1:N], balanced_l * y[i,h] <=  f[i,h])
    
    #@constraint(model, [h=1:S, j=1:N], f[h,j] <= balanced_u)
    #@constraint(model, [h=1:S, j=1:N], balanced_l * y[h,j] <=  f[h,j])

    for crossing in crossings
        i_1 = crossing[3]
        j_1 = crossing[4]
        i_2 = crossing[5]
        j_2 = crossing[6]
 
        println("Fix crossing between (", i_1, ",", j_1, ") and (", i_2, ",", j_2, ")")
        @constraint(model, y[i_1, j_1] + y[j_1, i_1] + y[i_2, j_2] + y[j_2, i_2] <= 1)
    
    end 

    optimize!(model)

    println("DONE OPTMISING")

    new_crossings = evaluate_termination(timelimit, model, x, f, y, N, T, instance, zone, data, manifold_array, folder_name, balanced_u, balanced_l)
    println("Crossings", new_crossings)
    
    if !isempty(new_crossings) && it < 3
        combined_crossings = deepcopy(crossings)

        for crossing in new_crossings
            push!(combined_crossings, crossing)

        end 

        WriteResults.write_crossings_to_exsiting_txt(folder_name*"/crossings", new_crossings)

        run_mip(instance, zone, number_of_manifolds, pressure, strings, timelimit, folder_name, combined_crossings, it+1)

    end

end


function evaluate_termination(timelimit, model, x, f, y, N, T, instance, zone, data, manifold_array, folder_name, balanced_u, balanced_l)

    # Write results
    termination = ""
    best_objective = NaN
    best_bound = NaN
    sol_mip = nothing  # Default to nothing in case no solution is found
    sol_sizing = nothing
    cal_sizing = nothing
    cal_mip = nothing
    crossings = []

    WriteResults.create_folder(folder_name)

    if termination_status(model) == MOI.TIME_LIMIT && primal_status(model) == MOI.FEASIBLE_POINT
        termination = "Time"
        # Retrieve the best objective value (upper bound)
        best_objective = objective_value(model)
        # Retrieve the best bound (lower bound)
        best_bound = objective_bound(model)

        # save solution function
        sol_mip, termination, best_objective, best_bound = write_solution(model, data, x, f, termination, best_objective, best_bound, N, T)
        sol_mip, cal_mip = pressureDropCalculations.calculate_mip(sol_mip, data, manifold_array)
        WriteResults.write_json(sol_mip, cal_mip, data, folder_name*"/"*instance*zone*"_"*"mip") 
       
        sol_sizing, cal_sizing = pressureDropCalculations.run(sol_mip, data, manifold_array, true)
        println("FEASIBLE: ", cal_sizing.feasible)
        sol_sizing.cost, _ = improvementHeuristic.cost_function(sol_sizing, data, manifold_array)
        WriteResults.write_json(sol_sizing, cal_sizing, data, folder_name*"/"*instance*zone*"_"*"sizing")

        crossings = improvementHeuristic.check_crossings(sol_sizing, data, instance, zone)

        if isempty(crossings)
            sol_balancing = deepcopy(sol_sizing)  
            sol_balancing.cost, _ = improvementHeuristic.cost_function(sol_balancing, data, manifold_array, true, true) 
            row = [instance, zone, termination, best_bound, timelimit, best_objective, sol_sizing.cost, cal_mip.feasible, cal_sizing.feasible, balanced_u, balanced_l, (maximum(values(sol_balancing.balanced)) - minimum(values(sol_balancing.balanced)))]
        else
            return crossings 
        end 


    elseif termination_status(model) == MOI.OPTIMAL
        termination = "Optimal"
        # Retrieve the best objective value (upper bound)
        best_objective = objective_value(model)

        # save solution function
        sol_mip, termination, best_objective, _ = write_solution(model, data, x, f, termination, best_objective, best_bound, N, T)
        sol_mip, cal_mip = pressureDropCalculations.calculate_mip(sol_mip, data, manifold_array)
        WriteResults.write_json(sol_mip, cal_mip, data, folder_name*"/"*instance*zone*"_"*"mip") 

        sol_sizing, cal_sizing = pressureDropCalculations.run(sol_mip, data, manifold_array, true)
        println("FEASIBLE: ", cal_sizing.feasible)
        sol_sizing.cost, _ = improvementHeuristic.cost_function(sol_sizing, data, manifold_array)
        WriteResults.write_json(sol_sizing, cal_sizing, data, folder_name*"/"*instance*zone*"_"*"sizing")

        crossings = improvementHeuristic.check_crossings(sol_sizing, data, instance, zone)

        if isempty(crossings) 
            sol_balancing = deepcopy(sol_sizing)  
            sol_balancing.cost, _ = improvementHeuristic.cost_function(sol_balancing, data, manifold_array, true, true) 
            row = [instance, zone, termination, best_bound, timelimit, best_objective, sol_sizing.cost, cal_mip.feasible, cal_sizing.feasible, balanced_u, balanced_l, (maximum(values(sol_balancing.balanced)) - minimum(values(sol_balancing.balanced)))]       
        else
            return crossings 
        end 

    else # Handle infeasibility
        termination = string(termination_status(model))
        sol_mip, termination, _, _ = write_solution(model, data, x, f, termination, best_objective, best_bound, N, T)
        row = [instance, zone, termination, best_bound, timelimit, best_objective, NaN, NaN, NaN, balanced_u, balanced_l, NaN]
    
    end

    WriteResults.write_to_exsiting_txt(folder_name*"/all_instances", row, true)

    return crossings

end 

function find_index_from_edge(data, i, j)
    index = (j - 1) * data.n + i
    return index
end

function write_solution(model, data, x, f, termination, best_objective, best_bound, N, T)
    sol = Solution()  # Initialize a new Solution object

    # Populate the Solution object based on x, f, y
    for i in 1:N
        for j in 1:N
            for t in 1:T
                if value(x[i, j, t]) > 0
                    idx = find_index_from_edge(data, j, i)
                    push!(sol.edges, idx)
                    #haskey(sol.adjacency_list, i) ? push!(sol.adjacency_list[i], j) : (sol.adjacency_list[i] = [j])
                    haskey(sol.adjacency_list, j) ? push!(sol.adjacency_list[j], i) : (sol.adjacency_list[j] = [i])
                    sol.pipetype[idx] = t
                    sol.flowrate[idx] = round(value(f[i, j]))
                end
            end
        end 
    end


    sol.cost = objective_value(model)
    sol.nodes = Set(1:N)

    return sol, termination, best_objective, best_bound
end


#= folder_name = "11-12_S_MIP_5min_test"
instance = "A"
zone = "2"
number_of_manifolds = 1
pressure = 26
timelimit = 60*30
strings = true 
crossings = [] =#

#run_mip(instance, zone, number_of_manifolds, pressure, strings, timelimit, folder_name, [])
run_mip(ARGS[1], ARGS[2], parse(Int64,ARGS[3]), parse(Int64,ARGS[4]), parse(Bool, ARGS[5]), parse(Int64,ARGS[6]), ARGS[7], [])
