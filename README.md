# Code for paper: Machine-Learning Component for Multi-Start Metaheuristics to Solve the Capacitated Vehicle Routing Problem



## Overview
This repository contains the code and the data of the paper:
**Machine-Learning Component for Multi-Start Metaheuristics to Solve the Capacitated Vehicle Routing Problem**. *Pre-print available at: [https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4120102](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4120102)*. 

# Folder structure
Files and shell scripts are structured as follows: 
``` bash
${repo_root_directory}         # e.g. amazon-sagemaker-amazon-routing-challenge-sol
└── GRASP_MLP/                       # GRASP with VND local search and ISC component
    └── src/                 # C++ files
└── ONLY_GRASP/                       # GRASP only with VND local search
    └── src/                 # C++ files
└── RR_GRASP_MLP/                       # GRASP with Ruin-and-Recreate local search and ISC component
    └── src/                 # C++ files
└── RR_ONLY_GRASP/                       # GRASP only with Ruin-and-Recreate local search
    └── src/                 # C++ files   
└── instances/                      # Instances of CVRP with CVRPLIB format
└── results/                       # Results folder

├── generate_timed_experiment.sh   # Bash script to reproduce our experiments of Ruin-and-Recreate local search results on a Linux Computer with Bash and SLURM
├── generate_vnds_timed_experiment.sh       # Bash script to reproduce our experiments of VND local search results on a Linux Computer with Bash and SLURM
```
## Requirements

Intel oneAPI
DPC++/C++ Compiler, Boost 1.78 library

## Running the Experiment

To run the experiment first compile the GRASPs that you want to use (In this case, we show the code examples for the Ruin-and-Recreate variants.
```
cd RR_GRASP_MLP
make
cd ..
cd RR_ONLU_GRASP
make
```

Now run the GRASP models to solve a CVRP instance by running the commands

The format is:
```
SEED = Number of seed to use with the random numbers generator
METHOD = 1 Always
INSTANCE = name of the cvrp instance to solve, such as X-n979-k58.vrp
./grasp_mlp $SEED $METHOD $INSTANCE
```


First, GRASP with only Ruin-and-Recreate
```
cd RR_ONLY_GRASP/
./only_grasp 10 1 X-n101-k25.vrp
cd ..
```
Second, GRASP with Ruin-and-Recreate and ISC
```
cd RR_GRASP_MLP/
./grasp_mlp 10 1 X-n101-k25.vrp
cd ..
```

The results will be stored in the 
``` bash
└── results/ 
```

Folder


