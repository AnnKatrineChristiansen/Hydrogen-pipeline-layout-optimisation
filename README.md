# Optimisation of pipeline layout for hydrogen-producing wind farms

This repository contains the code for a master thesis, optimising the pipeline layouts of hydrogen-producing wind farms using metaheuristics and exact methods. 

The repository contains the following files: 
* A gas pipeline hydraulic model with calculations of gas velocity and pressure drop.
    * `pressureDropCalculations.jl`: Calculations of gas velocity and pressuredrop combined to a pipe sizing heuristic.
    * `fixLayouts.jl`: Solves limitations in current pipe sizing heuristic.
* A sweep multistart heuristic to find initial soltions, it uses either the Minimum Spanning Tree or Nearest Neighbour algorithm to connect turbines based on the topology considered.
    * `mainSweep.jl`: Executes the sweep multistart heuristic.
    * `sweepConstructionHeuristic.jl`: Combines code from below files into a sweep multistart construction heuristic.
        * `Kmeans.jl`: Groups turbines to *K* manifolds.
        * `sweep.jl`: Orders turbines by angular position to starting turbine and sweep line direction.
        * `minimumSpanningTree.jl`: Connects a group of turbines using minimum spanning tree.
        * `nearestNeighbour.jl`: Connects a group of turbines using nearest neighbour.
        * `graph.jl`: Converts undirected graph to directed graph based on a root node.
* Two metaheuristics used to improve upon the solutions optained from the sweep multistart heuristic.
    * `main.jl`: Executes the metaheuristics Simulated Annealing or Variable Neighbourhood Search.
    * `improvementHeuristics.jl`: Comprises Simulated Annealing and Variable Neighbourhood Search.
    * `tuningMetaheuristics.jl`: Tuning of the metaheuristics.
* An exact Mixed Integer Linear Programming (MILP) model to optimise pipeline layouts.
    * `MIP.jl`: MILP modelformulation using JuMP.
* Experiments to align problem with real-world wind farms, extending the baseline scenario.
    * `experimentLimitedPipes.jl`: Limiting available pipe choices to 2-3 pipes.
    * `experimentManifoldPressure.jl`: Includes export pipes in addition to inter-array pipes.
    * `MIP_balanced.jl`: Enforces balanced layouts in the MILP formulation.
* Miscellaneous.
    * `readData.jl`: Reads wind farm instances and available pipe types.
    * `writeResults.jl`: Writes results to log files.



