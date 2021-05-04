#!/bin/bash

# Perform an affine transformation to a sequence of meshes

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

MSERVER="/meshlab/distrib/meshlabserver"
MLABSCRIPT="/MeshTracker/resource/scripts/abs_tform.mlx"

first_frame=00001
last_frame=00075
ext="ply"
length=$(expr $last_frame - $first_frame)
file_pref="Frame_"
file_postf="_hd_t."${ext}

for frame in $(seq -w $first_frame $last_frame)
do

    fname=${file_pref}${frame}${file_postf}
    xvfb-run -a -s "-screen 0 800x600x24"  $MSERVER -i $fname -o $fname -s $MLABSCRIPT

done
