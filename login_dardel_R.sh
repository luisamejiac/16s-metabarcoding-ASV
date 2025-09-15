#!/bin/bash

#Go to directory 
cd scripts/16s-metabarcoding-ASV/output

# Load modules
module load PDC/24.11


# Set safe locale to avoid R crashing
export LANG=C
export LC_ALL=C


# Start R interactively
R

ENDSSH