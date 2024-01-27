module nearestNeighbor

mutable struct solution
    nodes::Set{Int64}
    adjacency_list::Dict{Int64, Vector{Int64}} 
    edges::Set{Int64}
    pipetype::Dict{Int64, Int64} #Array{Int32,1} # pipetype used between node i and j
    flowrate::Dict{Int64, Int64} #Array{Int32,1}  # [kg H2/day] % The flowrate which the pipe needs to transport
    cost:: Float64
    crossings::Int64
    m_connections::Int64 #Dict{Int64, Int64}
    t_connections::Int64
    balanced::Dict{Int64, Int64}
    digging_cost::Float64
    solution() = new(Set{Int32}(),Dict{Int64, Vector{Int64}}(), Set{Int64}(),  Dict(), Dict(), 0, 0, 0, 0, Dict(), 0)
end 

function get_nearest_neighbor(data, visited, node, node_list)

    distance = Inf
    next = 0
    for i in node_list 
        if !(i in visited) && distance > data.distance[node, i]
            distance = data.distance[node, i]
            next = i
        end
    end

    return next
end 

function run(data, node_list, substation)
    
    sol = solution()
    visited = Set{Int64}()
    route = zeros(Int64,length(node_list))

    push!(sol.nodes, substation)
    push!(visited, substation)
    route[1] = substation

    for i in 1:length(node_list)-1
       j = get_nearest_neighbor(data, visited, route[i], node_list)

       route[i+1] = j
       push!(visited, j)
       sol.adjacency_list[route[i]] = Vector{Int64}()
       append!(sol.adjacency_list[route[i]], j)

    end

    return sol
end 
end 