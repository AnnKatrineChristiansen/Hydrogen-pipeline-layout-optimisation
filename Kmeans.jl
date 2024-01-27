module Kmeans

# Function to calculate the Euclidean distance between two points
function euclidean_distance(point1, point2)
    return sqrt(sum((x - y) ^ 2 for (x, y) in zip(point1, point2)))
end

function initialize_centroids(coords, data, k)
    substations = rand(1:data.n,k)
    centroids = [coords[s,:] for s in substations]
    return centroids
end

# Function to assign each data point to the nearest centroid
function assign_to_clusters(coords, data, substations, centroids, K, fixed_centroids)

    clusters = [Vector{Int64}() for k=1:K]
    
    for node in 1:data.n

        if fixed_centroids
            closest_centroid_index = argmin([data.distance[node, substation] for substation in substations])
        else
            closest_centroid_index = argmin([euclidean_distance(coords[node,:], centroid) for centroid in centroids])
        end 

        push!(clusters[closest_centroid_index], node)

    end
    
    return clusters
end

# Function to update centroid locations
function update_centroids(coords, clusters, centroids)
    new_centroids = Vector{Vector{Float64}}(undef, length(clusters))
    
    for (i, cluster) in enumerate(clusters)
        cluster_size = length(cluster)
        if cluster_size == 0
            new_centroids[i] = cluster[1]
        else
            new_centroid = sum([coords[i,:] for i in cluster])/cluster_size
            new_centroids[i] = new_centroid
        end
    end
    
    return new_centroids
end

# Function to check if the centroids have converged
function centroids_converged(old_centroids, new_centroids, tolerance=1e-4)
    return all(euclidean_distance(old, new) < tolerance for (old, new) in zip(old_centroids, new_centroids))
end


function run(data, coords, substations = [], fixed_centroids = true, k = 2)

    # Initialise centroids
    if  isempty(substations)
        centroids = initialize_centroids(coords, data, k)
    else
        centroids = [coords[s,:] for s in substations]
        k = length(centroids)
    end 
    
    clusters = nothing
    converged = false 
    while converged == false 
        
        # Assign clusters 
        clusters = assign_to_clusters(coords, data, substations, centroids, k, fixed_centroids)

        # Update centroids
        if fixed_centroids 
            new_centroids = deepcopy(centroids)
        else
            new_centroids = update_centroids(coords, clusters, centroids)
        end

        # Check if centroids have converged
        converged = centroids_converged(centroids, new_centroids)

        # Assign new centroids
        centroids = deepcopy(new_centroids)
    end  

    return centroids, clusters 
end 

end 
