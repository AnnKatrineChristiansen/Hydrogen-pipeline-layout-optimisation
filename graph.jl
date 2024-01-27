module Graph
### GRAPH
#using Graphs, GraphPlot #, Colors, ColorSchemes, DelimitedFiles, SimpleWeightedGraphs

mutable struct graph
    weighted_adjacency_list::Dict{Int64, Dict{Int64, Float64}} 
    nodes::Int64
    graph(nodes) = new(Dict(key => Dict{Int64, Float64}() for key in nodes) , length(nodes))
end 

function create(data, node_list)
    g = graph(node_list)

    for i in node_list
        for j in node_list
            if i!= j 
                dist = data.distance[i,j]
                g.weighted_adjacency_list[i][j] = dist
                g.weighted_adjacency_list[j][i] = dist 
            end 
        end 
    end 
    return g
end


using Graphs, GraphPlot, SimpleWeightedGraphs
function plot(coords, sol)
        # Create an empty undirected graph with the number of nodes
    g = SimpleGraph(size(coords,1))

    # Add edges from the adjacency list
    for (node, neighbors) in sol.adjacency_list
        for neighbor in neighbors
            add_edge!(g, node, neighbor)
        end
    end

    x_coords = coords[:, 1]
    y_coords = coords[:, 2]

    gplot(g, x_coords, y_coords, nodelabel=1:nv(g), nodefillc="skyblue")
end 



function convert_to_directed_graph(root_node::Int64,adjacency_list::Dict{Int64, Vector{Int64}})
    directed_graph = Dict{Int64, Vector{Int64}}()
    visited = Set{Int64}()
    
    function dfs(node, parent = nothing)
        push!(visited, node)
        for neighbor in get(adjacency_list, node, [])
            if !(neighbor in visited)
                push!(get!(directed_graph, node, []), neighbor)
                dfs(neighbor, node)
            end
        end
    end
    
    dfs(root_node)  
    return directed_graph
end

end 