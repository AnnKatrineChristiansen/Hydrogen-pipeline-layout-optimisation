
module minimumSpanningTree

mutable struct solution
    nodes::Set{Int64}
    adjacency_list::Dict{Int64, Vector{Int64}} 
    edges::Set{Int64}
    pipetype::Dict{Int64, Int64} #Array{Int32,1} # pipetype used between node i and j
    flowrate::Dict{Int64, Int64} #Array{Int32,1}  # [kg H2/day] % The flowrate which the pipe needs to transport
    cost::Float64
    crossings::Int64
    m_connections::Int64 #Dict{Int64, Int64}
    t_connections::Int64
    balanced::Dict{Int64, Int64}
    digging_cost::Float64
    solution() = new(Set{Int32}(),Dict{Int64, Vector{Int64}}(), Set{Int64}(),  Dict(), Dict(), 0, 0, 0, 0, Dict(), 0)
end 

function find_min_edge(g, sol::solution)
    fromNode = nothing
    toNode = nothing 
    minWeight = Inf
    
    for node in sol.nodes 
        # Find key with smallest distance 
        k = reduce((x, y) -> g.weighted_adjacency_list[node][x] ≤ g.weighted_adjacency_list[node][y] ? x : y, keys(g.weighted_adjacency_list[node])) # smallest key
        weight_k = g.weighted_adjacency_list[node][k]
        if (weight_k < minWeight) 
            fromNode = node
            toNode = k
            minWeight = weight_k
        end

    end 
    return (fromNode, toNode, minWeight)
end

function run(g, substation)
    # Define solution
    sol = solution()
    # Select the first node to start the tree
    push!(sol.nodes,substation)
    sol.adjacency_list[substation] = Vector{Int64}()

    while (length(sol.nodes) < g.nodes)

        #Find the minimum edge in the set of edges
        (fromNode, toNode, minWeight) = find_min_edge(g, sol)

        # Remove edges in both directions from nodes connected
        for node in sol.nodes
            delete!(g.weighted_adjacency_list[toNode], node)
            delete!(g.weighted_adjacency_list[node], toNode)
        end 

        # Add the toNode to the MST - fromNode is already added
        push!(sol.nodes, toNode)

        # Add edge in both directions to adjacency list
        if haskey(sol.adjacency_list, toNode) == false
            sol.adjacency_list[toNode] = Vector{Int64}()
        end 

        append!(sol.adjacency_list[toNode], fromNode)
        append!(sol.adjacency_list[fromNode], toNode)

    end 

    return sol 

end 

end 



 