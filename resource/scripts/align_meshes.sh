#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

MSERVER="/meshlab/src/build/distrib/meshlabserver"
MLABSCRIPT="/MeshTracker/resource/scripts/abs_tform.mlx"

# To actually use meshlabserver you need to use a virtual framebuffer and export the display 
# Xvfb :100 &
# export DISPLAY=:100.0

first_frame=00001
last_frame=00075
ext="ply"
length=$(expr $last_frame - $first_frame)
file_pref="Frame_"
file_postf="_hd_t."${ext}

for frame in $(seq -w $first_frame $last_frame)
do

    fname=${file_pref}${frame}${file_postf}
    $MSERVER -i $fname -o $fname -s $MLABSCRIPT

done