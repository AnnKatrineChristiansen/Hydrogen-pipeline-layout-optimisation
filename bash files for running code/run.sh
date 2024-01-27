#!/bin/bash

#Gå ind i mappen med cd “”
#source ~/.bashrc
#Initialize: chmod +x run.sh
#Kør: ./run.sh      


# Define the instances (name, zone, no. substations)
#instances=("A 1 1" "A 2 1") 
instances=("A 1 1" "A 2 1" "B 1 2" "C 1 1" "C 2 1" "D 1 3" "E 1 1" "E 2 2" "E 3 1" "F 1 2" "F 2 1" "G 1 3" "H 1 2" "H 2 1" "I 1 4" "J 1 2" "J 2 2" "J 3 1")
#instances=("B 1 2" "D 1 3" "E 2 2" "G 1 3" "H 1 2" "J 2 2")
# Common parameters
clockwise=true
strings=false
heuristic="VNS"
timelimit=1800
solution_filename="01-02_C_B_VNS_Digging_30min_510"
balacing=false
digging_price=510

# Loop through the instances and run the Julia script
for instance in "${instances[@]}"; do
    # Split the instance string into filename, zone, and substations
    read -r -a params <<< "$instance"

    filename="${params[0]}"
    zone="${params[1]}"
    number_of_substations="${params[2]}"

    # Run the Julia script with the specified parameters
    start_time=$SECONDS
    julia main.jl "$filename" "$zone" "$number_of_substations" "$clockwise" "$strings" "$timelimit" "$solution_filename" "$heuristic" "$balacing" "$digging_price"
    elapsed=$(( SECONDS - start_time ))
    echo "Elapsed time for $filename $zone: $elapsed seconds"
done
