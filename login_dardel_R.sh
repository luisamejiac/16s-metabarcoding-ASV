#!/bin/bash

# Load modules
module load PDC/24.11

# Set safe locale to avoid R crashing
export LANG=C
export LC_ALL=C


# Start R interactively
R

ENDSSH
