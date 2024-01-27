#!/bin/bash
#SBATCH --job-name=Measuringtime_002
#SBATCH --output=/rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/output/output_%5a.txt
#SBATCH --error=/rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/output/output_%5a.txt
#SBATCH --time=23:00:00
#SBATCH --partition=c18m
#SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --mem-per-cpu=3882M
#SBATCH --account=jara0052


bash /rwthfs/rz/cluster/home/bo791269/Software/AntiprotonAnalysis/Recipes/Measuringtime/step002/jobs/slurmjob_$(printf "%05d" ${SLURM_ARRAY_TASK_ID}).sh
