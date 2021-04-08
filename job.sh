#!/bin/bash
# Run snakemake with a singularity container
## Sylvain SCHMITT
## 25/02/2020

#SBATCH -J snakemake
#SBATCH --time=00:05:00
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 128000MB
#SBATCH -o snakemake.%N.%j.out
#SBATCH -e snakemake.%N.%j.err

# Module
module load snakemake

# Configuration
CLUSTER_CONFIG=cluster.json # Cluster config file with parameters for each rule of the Snakefile
CLUSTER="sbatch --mem={cluster.mem} --ntasks-per-node {cluster.npernode} -t {cluster.time} -n {cluster.ntasks} -c {cluster.c} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err"
# Slurm command to launch each new jobs (aka rule) create by snakemake. 
# Arguments values are parameters of the cluster config file.
# Use at most N cores in parallel (default: 1). If N is omitted or 'all', the limit is set to the number of available cores.
MAX_CORES=100

# Preparation
mkdir -p snake_subjob_log # Create a log directory for all the slurm output files
snakemake -s Snakefile --dag | dot -Tpng > dag.png # Create the directive acyclic graph of the workflow

# Run
snakemake -s Snakefile --use-singularity -j $MAX_CORES --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER"
## Launch the workflow : -s Snakefile --use-singularity launch the container 
## --cluster-config cluster configuration file for each rule in the cluster --cluster sbacth shell command

# Report 
snakemake -s Snakefile --report report.html 

## Informations
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job ID:' $SLURM_JOB_ID
echo 'Array task ID:' ${SLURM_ARRAY_TASK_ID}
echo 'Number of nodes assigned to job:' $SLURM_JOB_NUM_NODES
echo 'Total number of cores for job:' $SLURM_NTASKS
echo 'Number of requested cores per node:' $SLURM_NTASKS_PER_NODE
echo 'Nodes assigned to job:' $SLURM_JOB_NODELIST
echo 'Directory:' $(pwd)
## Detail Information:
echo 'scontrol show job:'
scontrol show job $SLURM_JOB_ID
echo '########################################'
