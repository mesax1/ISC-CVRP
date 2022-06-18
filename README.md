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
    └── complete_grasp/            #Store results for each iteration of the GRASPS with only local search
    └── consolidated_results/      #Store the consolidated results of the experiments executed using the bash scripts "generate_timed_experiment.sh" or "generate_vnds_timed_experiment.sh
    └── evolution/                 #Store the evolution/convergence of the improvement of a solution during the GRASP execution
    └── features_datasets/         #Store the features that describe initial solutions generated on each iteration of the GRASPS
    └── solutions_datasets/        #Store results for each iteration of the GRASPS with local search and ISC during training phase
    └── validation_solutions/      #Store results for each iteration of the GRASPS with local search and ISC during classification phase
    

├── generate_timed_experiment.sh            # Bash script to reproduce our experiments of Ruin-and-Recreate local search results on a Linux Computer with Bash and SLURM
├── generate_vnds_timed_experiment.sh       # Bash script to reproduce our experiments of VND local search results on a Linux Computer with Bash and SLURM
├── consolidate_rr_results.py       #Python script to consolidate results of the generate_timed_experiment.sh 
├── consolidate_vnd_results.py      #Python script to consolidate results of the generate_vnds_timed_experiment.sh
```
## Requirements

Intel oneAPI DPC++/C++ Compiler

Boost 1.78 library

## Running the Experiment

To run the experiment first compile the GRASPs that you want to use (In this case, we show the code examples for the Ruin-and-Recreate variants.
```
cd RR_GRASP_MLP
make
cd ..
cd RR_ONLY_GRASP
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

The results will be stored in the folder
``` bash
└── results/ 
```

## Reproducing the results of our work

On a Linux computer, with Bash and Slurm, such as the Apolo supercomputer from Universidad EAFIT [https://www.eafit.edu.co/apolo](https://www.eafit.edu.co/apolo):

1. Clone this repository
```
mkdir GRASP_ISC
cd GRASP_ISC
git clone https://github.com/mesax1/ISC-CVRP.git
```

If necessary, load C++ modules, with the Intel compiler, boost library, Python
```
module load intel/2022_oneAPI-update1 boost-1.78.0-gcc-11.2.0-ga5km6r
```

2. Compile the different GRASP codes
```
cd RR_GRASP_MLP
make
cd ..
cd RR_ONLY_GRASP
make
cd GRASP_MLP
make
cd ..
cd ONLY_GRASP
make
```

3. Execute the experiments bash script
3.1 Ruin-and-recreate local search
```
generate_timed_experiment.sh
```
3.2 VND local search
```
generate_vnds_timed_experiment.sh
```

4. When the experiments end, execute the results consolidation script
4.1 Ruin-and-recreate local search
```
consolidate_rr_results.py
```
4.2 VND local search
```
consolidate_vnd_results.py
```
5. The results are stored in the /results/consolidated_results/ folder
