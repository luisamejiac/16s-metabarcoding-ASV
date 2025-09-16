#!/bin/bash
ssh almcama@dardel.pdc.kth.se <<'ENDSSH'
cd scripts/16s-metabarcoding-ASV/output #Go to directory 
module load PDC/24.11 # Load modules
export LANG=C # Set safe locale to avoid R crashing
export LC_ALL=C
R # Start R interactively
ENDSSH
