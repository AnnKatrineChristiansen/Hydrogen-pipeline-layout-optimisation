module pressureDropCalculations
using CoolProp # You need to have the CoolProp.jl package installed
include("writeResults.jl")

mutable struct calculations
    pressuredrop::Dict{Int64, Float64}
    velocity::Dict{Int64, Float64}
    feasible::Bool
    calculations() = new(Dict(), Dict(), true)
end 

struct constants
    s::Int64  # Dimensionless elevation adjustment parameter. Set to 0.
    P_B::Int64 # Base pressure [Pa]
    T_B::Float64 # Base temperature [K]
    M_H2::Float64 # Molar mass H2 [kg/mol]
    M_air::Float64 # Molar mass Air [kg/mol]
    G::Float64 #  Gas specific gravity
    constants(s, P_B, T_B, M_H2, M_air, G) = new(s, P_B, T_B, M_H2, M_air, G)
end 

struct properties
    z::Float64
    density::Float64
    viscosity::Float64
    Q::Float64
    Q_normal::Float64
    properties(z, density, viscosity, Q, Q_normal) = new(z, density, viscosity, Q, Q_normal)
end 

function initalize_constants()
    s = 0 # % Dimensionless elevation adjustment parameter. Set to 0.
    P_B = 101353 # Base pressure [Pa]
    T_B = 15.56 + 273.15 # Base temperature [K]
    M_H2 = CoolProp.PropsSI("M", "T", T_B, "P", P_B, "H2") # Molar mass H2 [kg/mol]
    M_air = CoolProp.PropsSI("M", "T", T_B, "P", P_B, "Air") # Molar mass Air [kg/mol]
    G = M_H2 / M_air # Gas specific gravity
    return constants(s, P_B, T_B, M_H2, M_air, G)

end 

function define_properties(c::constants, temperature, outlet_pressure, flowrate)
    # properties
    z = CoolProp.PropsSI("Z", "T", temperature, "P", outlet_pressure * 1e5, "H2") # Compressibility factor of Hydrogen at this temp and pressure
    density = CoolProp.PropsSI("D", "T", temperature, "P", outlet_pressure * 1e5, "H2") # Compute the density at the segment. 
    viscosity = CoolProp.PropsSI("V", "T", temperature, "P", outlet_pressure * 1e5, "H2") # Dynamic viscosity of hydrogen
    density_normal = CoolProp.PropsSI("D", "T", c.T_B, "P", c.P_B, "H2") # [kg/m3] compute the density at the segment. Could also be calculated as " pressure [Pa] / (temperature * 4124.2) "
    # Actual volumetric flow rate
    Q = flowrate / density # [kg H2/day] / [kg/m3] = [m3/day] 
    # Normal flow rate 
    Q_normal = flowrate / density_normal # [kg H2/day] / [kg/m3] = [m3/day] 
    return properties(z, density, viscosity, Q, Q_normal)

end

function cumulative_flow(data, adjacency_list::Dict{Int64, Vector{Int64}}, node::Int64, visited::Set{Int64}=Set{Int64}())
    push!(visited, node)
    total_flow = data.turbine_flow_rate_daily    # Each node produces a flow of 1 unit
    for predecessor in get(adjacency_list, node, Int64[])
        if !(predecessor in visited)
            total_flow += cumulative_flow(data, adjacency_list, predecessor, visited)
        end
    end

    return total_flow
end

function calculate_flowrate(data, sol)
    for (node, neighbors) in sol.adjacency_list
        for neighbor in neighbors
            idx = find_index_from_edge(data, node, neighbor)
            sol.flowrate[idx] = cumulative_flow(data, sol.adjacency_list, neighbor)
            push!(sol.edges, idx)
        end
    end 
    return sol
end 

function calculate_friction_factor(internal_roughness, reynolds_number, pipe_inner_diameter)
    # Calculate friction factor by trial and error
    diff_percentage = 1000.0
    f_guess = 0.02
    tol = 0.5
    while abs(diff_percentage) > tol
        left_side = -2 * log10(internal_roughness / (3.7 * (pipe_inner_diameter / 1000)) + 2.51 / (reynolds_number * sqrt(f_guess)))
        global f = (1 / left_side)^2
        diff_percentage = ((f_guess - f) / f) * 100
        f_guess = f
    end 
    return f
end

function calculate_inlet_pressure_and_velocity(c::constants, segment_length, outlet_pressure, pipe_inner_diameter, flowrate, temperature, internal_roughness)
    # Define properties for given outlet pressure
    p = define_properties(c, temperature, outlet_pressure, flowrate)
    # Calculate the velocity at the outlet
    velocity = p.Q / (24 * 3600) / (π * (pipe_inner_diameter / 2 / 1000)^2) # Velocity at the outlet
    # Calculate Reynolds number
    reynolds_number = p.density * velocity * (pipe_inner_diameter / 1000) / p.viscosity # [-]
    # Calculate the friction factor
    f = calculate_friction_factor(internal_roughness, reynolds_number, pipe_inner_diameter)
    # Calculate the inlet pressure
    inlet_pressure = sqrt((p.Q_normal / (1.1494e-3 * (c.T_B / c.P_B) * (pipe_inner_diameter^2.5)))^2 * c.G * temperature * segment_length * p.z * f + exp(c.s) * (outlet_pressure * 1e5)^2)
    # Converting the units to bar_a
    inlet_pressure /= 1e5
    return inlet_pressure, velocity 

end

function time_pressure(c::constants, p, segment_length, outlet_pressure, pipe_inner_diameter, temperature, internal_roughness)
    # Define properties for given outlet pressure
    # Calculate the velocity at the outlet
    velocity = p.Q / (24 * 3600) / (π * (pipe_inner_diameter / 2 / 1000)^2) # Velocity at the outlet
    # Calculate Reynolds number
    reynolds_number = p.density * velocity * (pipe_inner_diameter / 1000) / p.viscosity # [-]
    # Calculate the friction factor
    f = calculate_friction_factor(internal_roughness, reynolds_number, pipe_inner_diameter)
    # Calculate the inlet pressure
    inlet_pressure = sqrt((p.Q_normal / (1.1494e-3 * (c.T_B / c.P_B) * (pipe_inner_diameter^2.5)))^2 * c.G * temperature * segment_length * p.z * f + exp(c.s) * (outlet_pressure * 1e5)^2)
    # Converting the units to bar_a
    inlet_pressure /= 1e5
    return inlet_pressure, velocity 

end

function time_velocity(c::constants, p, pipe_inner_diameter)
    # Calculate the velocity at the outlet
    velocity = p.Q / (24 * 3600) / (π * (pipe_inner_diameter / 2 / 1000)^2) # Velocity at the outlet
    
    return velocity 
end

function time_velocity2(c::constants, pipe_inner_diameter, temperature, outlet_pressure, flowrate)
    p = define_properties(c, temperature, outlet_pressure, flowrate)
    # Calculate the velocity at the outlet
    velocity = p.Q / (24 * 3600) / (π * (pipe_inner_diameter / 2 / 1000)^2) # Velocity at the outlet
    
    return velocity 
end

function minimize_pipe_calculate_velocity_and_pressure(c::constants, data, segment_length, outlet_pressure, flowrate, temperature)

    # Define properties for given outlet pressure for all diameters
    p = define_properties(c, temperature, outlet_pressure, flowrate)

    # Calculate the velocities at the outlet for all diameters using vectorized calculation
    velocities = p.Q ./ (24 * 3600) ./ (π * ((data.pipetype_inner_diameters ./ 2) ./ 1000).^2)

    # Indices under 25m/s
    valid_indices = findall(x -> x < data.max_velocity, velocities)

    inlet_pressure = Inf 
    while !isempty(valid_indices) && (inlet_pressure > data.max_pressure)
        
        # Choose optimal index as smallest possible pipe
        smallest_possible_index = argmax(velocities[valid_indices])
        pipe_index = valid_indices[smallest_possible_index]
        velocity = velocities[pipe_index]

        # Choose the diameter that matches optimal index
        optimal_diameter_mm = data.pipetype_inner_diameters[pipe_index]
        internal_roughness = data.internal_roughness[pipe_index]

        # Calculate Reynolds number, friction factor, and inlet pressure for the optimal diameter
        reynolds_number = p.density * velocity * (optimal_diameter_mm / 1000) / p.viscosity
        f = calculate_friction_factor(internal_roughness, reynolds_number, optimal_diameter_mm)
        inlet_pressure = sqrt((p.Q_normal / (1.1494e-3 * (c.T_B / c.P_B) * (optimal_diameter_mm^2.5)))^2 * c.G * temperature * segment_length * p.z * f + exp(c.s) * (outlet_pressure * 1e5)^2)
        inlet_pressure /= 1e5

        deleteat!(valid_indices, smallest_possible_index)
    
    end 

    # If not valid with any pipes - ensure that it returns values that are not feasible
    if isempty(valid_indices) && (inlet_pressure > data.max_pressure)
        velocity = minimum(velocities) 
        pipe_index = length(data.pipetype_inner_diameters)

        # Choose the diameter that matches optimal index
        optimal_diameter_mm = data.pipetype_inner_diameters[pipe_index]
        internal_roughness = data.internal_roughness[pipe_index]

        # Calculate Reynolds number, friction factor, and inlet pressure for the optimal diameter
        reynolds_number = p.density * velocity * (optimal_diameter_mm / 1000) / p.viscosity
        f = calculate_friction_factor(internal_roughness, reynolds_number, optimal_diameter_mm)
        inlet_pressure = sqrt((p.Q_normal / (1.1494e-3 * (c.T_B / c.P_B) * (optimal_diameter_mm^2.5)))^2 * c.G * temperature * segment_length * p.z * f + exp(c.s) * (outlet_pressure * 1e5)^2)
        inlet_pressure /= 1e5

    end 

    return inlet_pressure, velocity, pipe_index
end

function find_edge_from_index(data,index)
    i = (index-1) % (data.n) + 1
    j = floor(Int64, (index-1)/data.n) + 1
    return i, j
end

function find_index_from_edge(data, i, j)
    index = (j - 1) * data.n + i
    return index
end

function calculate_mip(sol, data, substations)
    # Initliaze constants and solution
    c = initalize_constants()
    sol = calculate_flowrate(data, sol)
    cal = calculations()

    for substation in substations
        # Define the substation as the first index, assign known pressure and find connected turbines
        cal.pressuredrop[substation] = data.start_outlet_pressure
        neighbor_list = deepcopy(sol.adjacency_list[substation]) # O
        node_list = ones(Int32,length(neighbor_list))*substation # I 

        it = 1
        while  (it <= length(neighbor_list))
            # find inlet and outlet index of segment based on iteration
            neighbor = neighbor_list[it]
            node = node_list[it]
            idx = find_index_from_edge(data, node, neighbor)

            # Find segment length, outlet_pressure, and pipe diameter
            pipe = sol.pipetype[idx]
            segment_length = data.distance[node, neighbor] 
            outlet_pressure = cal.pressuredrop[node]
            pipe_inner_diameter = data.pipetype_inner_diameters[pipe]
            flowrate = sol.flowrate[idx]
            temperature = data.temperature
            internal_roughness = data.internal_roughness[pipe]
            # Calculate the inlet pressure and outlet velocity
            inlet_pressure, velocity = calculate_inlet_pressure_and_velocity(c, segment_length, outlet_pressure, pipe_inner_diameter, flowrate, temperature, internal_roughness)
            
            cal.pressuredrop[neighbor] = inlet_pressure
            cal.velocity[idx] = velocity

            if (inlet_pressure > data.max_pressure) || (velocity > data.max_velocity)
                cal.feasible = false
                #return sol, cal
            end  

            # Find turbines connected to current inlet index
            next_node  = deepcopy(neighbor)
            # If there are turbines connected, append to i and o list
            if haskey(sol.adjacency_list, next_node)
                next_neighbor = deepcopy(sol.adjacency_list[next_node]) # her var den
                append!(neighbor_list, next_neighbor)
                append!(node_list, ones(length(next_neighbor))*next_node)
            end 
            it += 1
            
             
        end 
    end
    
    return sol, cal
end


function run(sol, data, substations, mip = false)
    # Initliaze constants and solution
    c = initalize_constants()
    sol = calculate_flowrate(data, sol)
    cal = calculations()

    for substation in substations

        # Define the substation as the first index, assign known pressure and find connected turbines
        cal.pressuredrop[substation] = data.start_outlet_pressure
        neighbor_list = deepcopy(sol.adjacency_list[substation]) 
        node_list = ones(Int32,length(neighbor_list))*substation  

        it = 1
        while  (it <= length(neighbor_list))
            # find inlet and outlet index of segment based on iteration
            neighbor = neighbor_list[it]
            node = node_list[it]
            idx = find_index_from_edge(data, node, neighbor)

            # Find segment length, outlet_pressure, and pipe diameter
            segment_length = data.distance[node, neighbor] 
            outlet_pressure = cal.pressuredrop[node]
            flowrate = sol.flowrate[idx]
            temperature = data.temperature
            inlet_pressure, velocity, pipetype = minimize_pipe_calculate_velocity_and_pressure(c, data, segment_length, outlet_pressure, flowrate, temperature)
            
            cal.pressuredrop[neighbor] = inlet_pressure
            cal.velocity[idx] = velocity
            sol.pipetype[idx] = pipetype

            if (inlet_pressure > data.max_pressure) || (velocity > data.max_velocity)
                cal.feasible = false
                
                if !mip
                    return sol, cal 
                end 

            end  

            # Find turbines connected to current inlet index
            next_node  = deepcopy(neighbor)
            # If there are turbines connected, append to neighbor and node list
            if haskey(sol.adjacency_list, next_node)
                next_neighbor = deepcopy(sol.adjacency_list[next_node]) # her var den
                append!(neighbor_list, next_neighbor)
                append!(node_list, ones(length(next_neighbor))*next_node)
            end 
            it += 1   
 
        end 
    end 
    
    return sol, cal

end



function run_fast(sol, cal, data, substations, edges_affected, allow_illegal)
    # Initliaze constants and solution
    c = initalize_constants()
    sol = calculate_flowrate(data, sol)
    #cal = calculations()
    neighbor_list = []
    node_list = []

    for substation in substations

        # Define the substation as the first index, assign known pressure and find connected turbines
        if haskey(sol.adjacency_list, substation)
            cal.pressuredrop[substation] = data.start_outlet_pressure
            neighbor_list = deepcopy(sol.adjacency_list[substation]) 
            node_list = ones(Int32,length(neighbor_list))*substation 
        
        else
            cal.feasible = false
            folder_name = "SA_error"
            println(folder_name)
            WriteResults.create_folder(folder_name)
            WriteResults.write_json(sol, cal, data, folder_name*"/substation_error")
            return sol, cal

        end 

        it = 1
        while  (it <= length(neighbor_list))
            # find inlet and outlet index of segment based on iteration
            neighbor = neighbor_list[it]
            node = node_list[it]
            idx = find_index_from_edge(data, node, neighbor)
            
            if idx in edges_affected
                # Find segment length, outlet_pressure, and pipe diameter
                segment_length = data.distance[node, neighbor] 
                outlet_pressure = cal.pressuredrop[node]
                flowrate = sol.flowrate[idx]
                temperature = data.temperature
                inlet_pressure, velocity, pipetype = minimize_pipe_calculate_velocity_and_pressure(c, data, segment_length, outlet_pressure, flowrate, temperature)
                
                cal.pressuredrop[neighbor] = inlet_pressure
                cal.velocity[idx] = velocity
                sol.pipetype[idx] = pipetype

                if (inlet_pressure > data.max_pressure) || (velocity > data.max_velocity)
                    cal.feasible = false

                    if allow_illegal == false
                        return sol, cal
                    end 
                end  

                # Find turbines connected to current inlet index
                next_node  = deepcopy(neighbor)
                # If there are turbines connected, append to neighbor and node list
                if haskey(sol.adjacency_list, next_node)
                    next_neighbor = deepcopy(sol.adjacency_list[next_node]) # her var den
                    append!(neighbor_list, next_neighbor)
                    append!(node_list, ones(length(next_neighbor))*next_node)
                end 
            end
            it += 1
            
 
        end 
    end 
    
    return sol, cal

end


end 

 

