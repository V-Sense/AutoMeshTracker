#!/bin/bash

matlab_script_path="/MeshTracker/resource/scripts"

# Run MATLAB script that takes feature vectors as input and outputs keyframe indices
echo "-------------------------------------------------"
echo "KEYFRAME EXTRACTION VIA SHAPE SIMILARITY"
echo "-------------------------------------------------"

descriptors=$1
feas=$2
min_region=10
matlab "$@" -nodesktop -nosplash -r "addpath(\"${matlab_script_path}\"); dot_similarity(\"${descriptors}\",\"${feas}\",$min_region); quit"

exit;
