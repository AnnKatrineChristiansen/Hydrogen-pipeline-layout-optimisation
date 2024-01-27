module readData
using DelimitedFiles, DataFrames, XLSX, Distances 

mutable struct data
    n::Int32 # number of turbines
    start_outlet_pressure::Float64 # [bar a]
    pipetype_inner_diameters::Array{Float64,1}
    internal_roughness::Array{Float64,1} # [m]  Internal roughness of a steel pipe
    pipe_costs::Array{Float64,1} # cost of different diameters
    distance::Array{Float64,2}
    temperature::Float64 # [K] 
    coords::Array{Float64,2}
    bigM::Int64
    max_pressure::Float64
    max_velocity::Float64
    turbine_flow_rate_daily::Int64
    max_m_connections::Int64
    max_t_connections::Int64
    digging_cost::Int64
    data(n, start_outlet_pressure, pipetype_inner_diameters, internal_roughness, pipe_costs, distance,temperature, coords, bigM) = new(n, start_outlet_pressure, pipetype_inner_diameters, internal_roughness, pipe_costs, distance,temperature, coords, bigM,  40, 25,  363*24, 10, 3, 0)
end 

function readInstance(filename, zone = 1)
    data = readdlm(filename, ',', Float64, skipstart=1)
    zone_data = data[data[:,5].==zone,:]
    coords = zone_data[:,1:2]
    n = size(coords, 1)
    dist_matrix = pairwise(Euclidean(), coords', coords')./1000 # From [m] to [km]
    elevation = zone_data[:,3]
    
    return n, coords, dist_matrix, elevation
end

function readPipeinfo(filename, pipe_idx = [1,2,3,4,5,6])
    pipe_info = DataFrame(XLSX.readtable(filename,"Sheet1"))
    rename!(pipe_info, :"inner diameter (inch)" => :"diameter_inch", :"absolute internal roughness (mm)" => :"abs_roughness_mm", :"cost pr. meter" => :"cost_per_m")
    pipe_info = pipe_info[pipe_idx,:] #without export pipes
    pipe_diameters = pipe_info.diameter_inch * 25.4 # From [inch] to [mm]
    pipe_roughness = pipe_info.abs_roughness_mm * 10^(-3) # From [mm] to [m]
    pipe_costs = pipe_info.cost_per_m

    return pipe_diameters, pipe_roughness, pipe_costs
end

function readExportPipeinfo(filename)
    pipe_info = DataFrame(XLSX.readtable(filename,"Sheet1"))
    rename!(pipe_info, :"inner diameter (inch)" => :"diameter_inch", :"absolute internal roughness (mm)" => :"abs_roughness_mm", :"cost pr. meter" => :"cost_per_m")
    pipe_info = pipe_info[7:12,:] #with export pipes
    pipe_diameters = pipe_info.diameter_inch * 25.4 # From [inch] to [mm]
    pipe_roughness = pipe_info.abs_roughness_mm * 10^(-3) # From [mm] to [m]
    pipe_costs = pipe_info.cost_per_m

    return pipe_diameters, pipe_roughness, pipe_costs
end

function run(filename, zone = 1, pipe_idx = [1,2,3,4,5,6])
    n, coords, dist_matrix, elevation = readInstance(filename, zone)
    pipetype_inner_diameters, internal_roughness, pipe_costs = readPipeinfo("Data/Pipeline info and cost data.xlsx", pipe_idx)
    pipe_costs = pipe_costs 
    start_outlet_pressure = 26
    temperature = 3 + 273.15
    
    bigM = 5*10e5
    
    return coords, data(n, start_outlet_pressure, pipetype_inner_diameters, internal_roughness, pipe_costs, dist_matrix, temperature, coords, bigM)
end


#pipe_info = DataFrame(XLSX.readtable("Data/Pipeline info and cost data.xlsx","Sheet1"))
#rename!(pipe_info, :"inner diameter (inch)" => :"diameter_inch", :"absolute internal roughness (mm)" => :"abs_roughness_mm", :"cost pr. meter" => :"cost_per_m")

#filename = "Data/instances/A/layout.txt"
#n, coords, dist_matrix, elevation = readInstance(filename)

end