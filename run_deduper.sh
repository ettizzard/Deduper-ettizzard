#!/bin/bash

#SBATCH --account=bgmp   ### change this to your actual account for charging
#SBATCH --partition=bgmp     ### queue to submit to
#SBATCH --job-name=dedupe   ### job name
#SBATCH --output=dedupe%j.out   ### file in which to store job stdout
#SBATCH --error=dedupe%j.err    ### file in which to store job stderr
#SBATCH --time=1-0                ### wall-clock time limit, in minutes
#SBATCH --nodes=1               ### number of nodes to use
#SBATCH --cpus-per-task=1       ### number of cores for each task

set -e

/usr/bin/time -v ./Tizzard_deduper.py -f sorted_survey_input.sam -o survey_output.sam -u STL96.txt 