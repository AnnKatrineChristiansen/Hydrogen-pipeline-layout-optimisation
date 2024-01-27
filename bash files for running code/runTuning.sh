#!/bin/bash

#Gå ind i mappen med cd “”
#source ~/.bashrc
#Initialize: chmod +x runTuning.sh
#Kør: ./runTuning.sh      


# Define the instances (name, zone, no. substations)
#instances=("A 1 1" "A 2 1") 
instances=("B 1 2")

# Common parameters
clockwise=true
strings=false
heuristic="SA"
timelimit=1800
solution_filename="12-18_C_B_SA_tuning"

# Loop through the instances and run the Julia script
for instance in "${instances[@]}"; do
    # Split the instance string into filename, zone, and substations
    read -r -a params <<< "$instance"

    filename="${params[0]}"
    zone="${params[1]}"
    number_of_substations="${params[2]}"

    # Run the Julia script with the specified parameters
    start_time=$SECONDS
    julia tuningMetaheuristics.jl "$filename" "$zone" "$number_of_substations" "$clockwise" "$strings" "$timelimit" "$solution_filename" "$heuristic"
    elapsed=$(( SECONDS - start_time ))
    echo "Elapsed time for $filename $zone: $elapsed seconds"
done
