module WriteResults
using JSON

global FOLDER_PATH = "results/" # "solutions/"

function find_edge_from_index(data,index)
    i = (index-1) % (data.n) + 1
    j = floor(Int64, (index-1)/data.n) + 1
    return i, j
end

function to_dict(sol, cal, data)
    # Use find_edge_from_index to convert integer indices back to edges (i, j)
    pipetype_dict = Dict(find_edge_from_index(data, k) => v for (k, v) in sol.pipetype)
    flowrate_dict = Dict(find_edge_from_index(data, k) => v for (k, v) in sol.flowrate)
    velocity_dict = Dict(string(find_edge_from_index(data, k)) => v for (k, v) in cal.velocity)
 
    
    return Dict(
        "adjacency_list" => sol.adjacency_list,
        "pipetype" => pipetype_dict,
        "flowrate" => flowrate_dict,
        "cost" => sol.cost,
        "pressuredrop" => cal.pressuredrop,
        "velocity" => velocity_dict,
        "feasible" => cal.feasible,
        "edges" => sol.edges
    )
end

function write_json(sol, cal, data, instance_name)
    json_str = JSON.json(to_dict(sol, cal, data))
    open(FOLDER_PATH * instance_name * ".json", "w") do f
        write(f, json_str)
    end
end


function write_txt(filename, results, metaheuristic)
    outFile = open(FOLDER_PATH * filename * ".txt", "w")

    # Write the header line           
    if metaheuristic == "SA"
        #header = ["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections", "BalanceDiff", "Temperature", "Alpha","Reheats"]
        header = join(["It", "Heuristic", "ElapsedTime","Bestcost", "Currentcost", "Crossing", "MConnections", "TConnections", "BalanceDiff","CostDigging", "Temperature", "Alpha","Reheats", "Balancing"], " ")
    
    elseif metaheuristic == "VNS"
    #   header = ["It", "Heuristic", "Elapsed Time","Best cost", "Current cost", "Crossing", "MConnections", "TConnections","BalanceDiff", 'K', "First Improvement", "Feasible"]
        header = join(["It", "Heuristic", "ElapsedTime","Bestcost", "Currentcost", "Crossing", "MConnections", "TConnections", "BalanceDiff", "CostDigging", "K", "FirstImprovement", "Feasible", "Balancing"], " ")
    
    elseif metaheuristic == "TuningVNS"
        header = join(["Heuristic", "ElapsedTime","Bestcost", "Feasible","K", "FirstImprovement", "AllowIllegal"], " ")

    elseif metaheuristic == "TuningSA"
        header = join(["Heuristic", "ElapsedTime","Bestcost", "Feasible", "Temperature", "Alpha"], " ")
   
    elseif metaheuristic == "exp_limited_pipes"
        header = join(["Instance", "Zone", "Cost", "Pipes", "Feasible", "Crossing", "MConnections", "TConnections",], " ")
        
    end 
    
    
    write(outFile, header * "\n")

    # Write each row
    for row in results
        row_str = join(string.(row), " ")  # Convert each element to a string
        write(outFile, row_str * "\n")
    end

    close(outFile)
end

function write_to_exsiting_txt(filname, row, mip=false, limited_pipes=false, sweep_test = false)

    if !isfile(FOLDER_PATH * filname * ".txt")
        # The file doesn't exist, so create it
        file = open(FOLDER_PATH * filname * ".txt", "w")
        
        # Optionally, write initial content to the file
        if (mip == false) && (limited_pipes == false) && (sweep_test == false)
            header = ["Instance", "Zone" , "SweepCost", "SweepCostUB", "SweepCostLB", "Heuristic", "BestCost", "ElapsedTime", "Timelimit", "Crossings", "MConnections", "TConnections", "BalanceDiff", "CostDigging"] 
        
        elseif (mip == false) && (limited_pipes == true)
            header = ["Instance", "Zone", "Cost", "Pipes", "Feasible", "Crossing", "MConnections", "TConnections"]
        
        elseif (sweep_test = true)
            header = ["Instance", "Zone", "Cost", "WorstCost", "Crossing", "MConnections", "TConnections", "RunTime"]
        
        else
            header = ["Instance", "Zone", "Termination", "LB", "Timelimit", "MIPObjective(UB)", "SizingCost", "MIPFeasible", "SizingFeasible"]
        end
        header_str = join(string.(header), " ")  # Convert each element to a string
        write(file, header_str * "\n")
        
        # Close the file when you're done
        close(file)
    
    end 
     
    # Open the file in append mode to add a new line
    file = open(FOLDER_PATH * filname * ".txt", "a")

    # Text you want to write as a new line
    row_str = join(string.(row), " ")  # Convert each element to a string
    write(file, row_str * "\n")

    # Close the file when you're done
    close(file)

end

function write_crossings_to_exsiting_txt(filname, rows)

    if !isfile(FOLDER_PATH * filname * ".txt")
        # The file doesn't exist, so create it
        file = open(FOLDER_PATH * filname * ".txt", "w")
        
        # Optionally, write initial content to the file
       
        header = ["Instance", "Zone", "i_1", "j_1", "i_2", "j_2"]

        header_str = join(string.(header), " ")  # Convert each element to a string
        write(file, header_str * "\n")
        
        # Close the file when you're done
        close(file)
    
    end 
     
    # Open the file in append mode to add a new line
    file = open(FOLDER_PATH * filname * ".txt", "a")

    # Text you want to write as a new line
    for row in rows
        row_str = join(string.(row), " ")  # Convert each element to a string
        write(file, row_str * "\n")
    end

    # Close the file when you're done
    close(file)

end


function create_folder(folder_name)
    # Create the folder if it doesn't exist
    if !isdir(FOLDER_PATH * folder_name)
        mkdir(FOLDER_PATH * folder_name)
    end
end 


function read_txt_file(file_path)
    try
        # Open the file
        file = open(FOLDER_PATH * file_path, "r")

        # Read lines from the file
        lines = readlines(file)

        # Close the file
        close(file)

        # Initialize a list to store the data
        data = []

        # Skip the first line (header)
        for line in lines[2:end]
            # Split the line into words
            words = split(line)

            # Convert numeric columns to integers
            row = [i == 1 ? words[i] : parse(Int, words[i]) for i in 1:length(words)]

            # Append the row to the data list
            push!(data, row)
        end

        return data
    catch e
        println("Error: $e")
        return []
    end
end

end 