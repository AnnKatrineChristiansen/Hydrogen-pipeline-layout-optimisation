module sweep

function calculate_polar_angles(coords, cluster, substation)
    x_sub = coords[substation, 1] 
    y_sub = coords[substation, 2]
    x = coords[:, 1] # Longitude
    y = coords[:, 2] # Latitude
    
    polar_angles = Dict{Int64, Float64}()
    for node in cluster
        if node != substation
            polar_angles[node] = atan(y[node] - y_sub, x[node] - x_sub)+pi/2
        end
    end  

    return sort(polar_angles; byvalue=true) #sort(polar_angles) #
end 


function order_polar_angles(polar_angles, starting_turbine, clockwise)
    
    nodes = collect(keys(polar_angles))
    idx_starting_node =  findall(x->x==starting_turbine, nodes)[1]
    ordered_polar_angles = []
    if clockwise 
        last = nodes[1:idx_starting_node-1]
        first = nodes[idx_starting_node:length(nodes)]
        append!(ordered_polar_angles, first)
        append!(ordered_polar_angles, last)
    else 
        first = reverse(nodes[1:(idx_starting_node-1)])
        last = reverse(nodes[(idx_starting_node+1):length(nodes)])
        append!(ordered_polar_angles, starting_turbine)
        append!(ordered_polar_angles, first)
        append!(ordered_polar_angles, last)
    end 

    return ordered_polar_angles
end 

end 