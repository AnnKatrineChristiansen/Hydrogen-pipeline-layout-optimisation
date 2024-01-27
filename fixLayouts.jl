include("minimumSpanningTree.jl")
include("readData.jl")
include("improvementHeuristic.jl")
include("pressureDropCalculations.jl")
include("writeResults.jl")


using JSON
global FOLDER_PATH = "results/" # "solutions/"

function find_edge_from_index(data,index)
    i = (index-1) % (data.n) + 1
    j = floor(Int64, (index-1)/data.n) + 1
    return i, j
end

function find_index_from_edge(data, i, j)
    index = (j - 1) * data.n + i
    return index
end

function from_dict_to_struct(dict, data)
   
    sol = minimumSpanningTree.solution()
    cal = pressureDropCalculations.calculations()
    
    # Adjacency list
    sol.adjacency_list = Dict(parse(Int, key) => value for (key, value) in dict["adjacency_list"])

    # Cost
    sol.cost = dict["cost"]

    # Pipetypes + velocity + flowrate
    pipetype_dict = Dict{Int64, Int64}()
    flowrate_dict = Dict{Int64, Int64}()
    velocity_dict = Dict{Int64, Float64}()
    edges = Set{Int64}()
    for (edge, pipetype) in dict["pipetype"]
        e = split(replace(edge, r"[() ]" => ""), ",")
        i = parse(Int, e[1])
        j = parse(Int, e[2])
        idx = find_index_from_edge(data, i, j)
        pipetype_dict[idx] = pipetype
        flowrate_dict[idx] = dict["flowrate"][edge]
        velocity_dict[idx] = dict["velocity"][edge]
        push!(edges, idx)
    end 
    sol.pipetype = pipetype_dict
    sol.flowrate = flowrate_dict
    sol.edges = edges
    cal.velocity = velocity_dict

    
    pressuredrop_dict = Dict{Int64, Float64}()
    for (edge, pressure) in dict["pressuredrop"]
        pressuredrop_dict[parse(Int, edge)] = pressure
    end 
    cal.pressuredrop = pressuredrop_dict
    cal.feasible = true

    return sol, cal 
end

function load_json_file(filename, instance, zone)
    # Read the contents of the JSON file
    json_content = read(FOLDER_PATH * filename * "/" *instance * string(zone) * "_sizing.json", String)

    # Parse the JSON content
    dict = JSON.parse(json_content)

    return dict
end


function find_path_to_node(data, sol, current_node, target_node, visited)
    if current_node == target_node
        return true
    end
    push!(visited, current_node)

    for neighbor in get(sol.adjacency_list, current_node, [])
        if neighbor ∉ visited
            edge_index = find_index_from_edge(data, current_node, neighbor)
            if find_path_to_node(data, sol, neighbor, target_node, visited)
                push!(path, edge_index)
                return true
            end
        end
    end
    return false
end

function find_leafs(data, sol)
    
    leafs = Set()
    
    for i in 1:data.n
        if (haskey(sol.adjacency_list, i) == false) 
            push!(leafs, i)
        end
    end 

    return leafs 
end


function pressure_limited_path(path, sol)

    previous_edge = 6
    l = length(path)
    for (num, edge) in enumerate(reverse(path))
        current_edge = sol.pipetype[edge]

        if (num == l) & (current_edge > 1)
            return false#true

        elseif previous_edge < current_edge
            return true
        else
            previous_edge = deepcopy(current_edge)
        end 

    end 

    return false

end

function find_paths_limited_by_pressure(sol, data, substation)

    affected_paths = []
    leafs = find_leafs(data, sol)

    for leaf in leafs
        global path = Vector{Int64}()
        path_exsists = find_path_to_node(data, sol, substation, leaf, Set{Int64}())

        if path_exsists
            result = pressure_limited_path(path, sol)

            if result
                push!(affected_paths, reverse(path))
            end
        end 

    end 
    return affected_paths
end



function process_paths(paths, data, sol, cal, substation)
    temp_sol = deepcopy(sol)
    temp_cal = deepcopy(cal)
    c = pressureDropCalculations.initalize_constants()


    for path in paths
        for edge in path
            if temp_sol.pipetype[edge] < 6
                temp_sol.pipetype[edge] += 1

                # Perform calculations for current edge
                i, j = find_edge_from_index(data, edge)
                temp_sol, temp_cal = perform_calculations(path, edge, i, j, data, temp_sol, temp_cal, c)

                #check if the path is still affected
                new_affected_paths = find_paths_limited_by_pressure(temp_sol, data, substation)
                
                if path in new_affected_paths
                    # If still affected, continue to the next edge
                    continue
                else
                    #if not affected, break and move to the next path
                    break
                end

            end
        end

        return temp_sol, temp_cal
    end
end

function perform_calculations(path, edge, i, j, data, temp_sol, temp_cal, c)
    
    segment_length = data.distance[i, j]
    outlet_pressure = temp_cal.pressuredrop[i]
    flowrate = temp_sol.flowrate[edge]
    temperature = data.temperature
    pipetype = temp_sol.pipetype[edge]
    pipe_inner_diameter = data.pipetype_inner_diameters[pipetype]
    internal_roughness = data.internal_roughness[pipetype]
    
    inlet_pressure, velocity = pressureDropCalculations.calculate_inlet_pressure_and_velocity(c, segment_length, outlet_pressure, pipe_inner_diameter, flowrate, temperature, internal_roughness)

    temp_cal.pressuredrop[j] = inlet_pressure
    temp_cal.velocity[edge] = velocity

    current_index = findfirst(==(edge), path)
    edges_affected = path[current_index+1:end]

    for edge in edges_affected

        i, j = find_edge_from_index(data,edge)

        segment_length = data.distance[i, j] 
        outlet_pressure = temp_cal.pressuredrop[i]
        flowrate = temp_sol.flowrate[edge]
        temperature = data.temperature

        inlet_pressure, velocity, pipetype = pressureDropCalculations.minimize_pipe_calculate_velocity_and_pressure(c, data, segment_length, outlet_pressure, flowrate, temperature)

        temp_cal.pressuredrop[j] = inlet_pressure
        temp_cal.velocity[edge] = velocity
        temp_sol.pipetype[edge] = pipetype
    end

    return temp_sol, temp_cal
end


function write_to_exsiting_txt(filname, row)

    if !isfile(FOLDER_PATH * filname * "/" *"layoutFix.txt")
        # The file doesn't exist, so create it
        file = open(FOLDER_PATH * filname * "/" *"layoutFix.txt", "w")
        
        # Optionally, write initial content to the file
       
        header = ["Instance", "Zone", "NumberAffectedPaths", "CurrentCost", "NewCost"]

        header_str = join(string.(header), " ")  # Convert each element to a string
        write(file, header_str * "\n")
        
        # Close the file when you're done
        close(file)
    
    end 
     
    # Open the file in append mode to add a new line
    file = open(FOLDER_PATH * filname * "/" *"layoutFix.txt", "a")

    # Text you want to write as a new line
    
    row_str = join(string.(row), " ")  # Convert each element to a string
    write(file, row_str * "\n")


    # Close the file when you're done
    close(file)

end


# read data
filename = "12-12_S_MIP_20bar_30min"
instance = "B"
zone = "1"
number_of_substations = 2


function run(filename, instance, zone, number_of_substations)

    substations = collect(1:number_of_substations)
    coords, data = readData.run("Data/instances/"*instance*"/layout.txt",  parse(Int, zone))

    # Load sol and cal 
    layout_dict = load_json_file(filename, instance, zone)
    sol, cal  = from_dict_to_struct(layout_dict, data)

    temp_sol = deepcopy(sol)
    temp_cal = deepcopy(cal)
    number_of_paths_affected = 0

    for substation in substations

        affected_paths = find_paths_limited_by_pressure(temp_sol, data, substation)
        
        if !isempty(affected_paths)
            number_of_paths_affected += length(affected_paths)
            temp_sol, temp_cal = process_paths(affected_paths, data, temp_sol, temp_cal, substation)
        end 

    end 
    # check cost
    current_cost, _ = improvementHeuristic.cost_function(sol, data, substations)
    new_cost, _ = improvementHeuristic.cost_function(temp_sol, data, substations)

    row = [instance, zone, number_of_paths_affected, current_cost, new_cost]
    write_to_exsiting_txt(filename, row)

    if current_cost != new_cost
        WriteResults.write_json(temp_sol, temp_cal, data, filename*"/"*instance*zone*"_"*"layoutFix")
    end 

    println("Instance ", instance, " zone ", zone)
    println("Current cost ", current_cost)
    println("New cost ", new_cost)
    println("")

end 


#= # read data
filename = "12-12_S_MIP_20bar_30min"
instance = "B"
zone = "1"
number_of_substations = 2 

run(filename, instance, zone, number_of_substations) =#

run(ARGS[1], ARGS[2], ARGS[3], parse(Int64,ARGS[4]))