#!/bin/bash

#SBATCH --job-name=SIMULATION_H0-1TTE-nm  # Job name
#SBATCH --array=1-200     # Number of jobs
#SBATCH --partition=standard    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/SIMULATION_H0-1TTE-nm/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/SIMULATION_H0-1TTE-nm/slurm-%x-%a-%J.out # Same file for output and error

R CMD BATCH --vanilla BATCH_SIMULATION_H0-1TTE-nm.R output/SIMULATION_H0-1TTE-nm/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd p:/Cluster/GPC/Article-inference-Ustatistic-Rao/
## clean outputs      ## rm -r ./output/SIMULATION_H0-1TTE-nm/*
## clean results      ## rm -r ./Results/SIMULATION_H0-1TTE-nm/*
## submission command ## sbatch SUBM_SIMULATION_H0-1TTE-nm.slurm

## submission output  ## Submitted batch job 43689 (time using C-c t: 03/30/21 9:34   ) 

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server

## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
