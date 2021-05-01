#!/bin/bash

MSERVER="/meshlab/src/build/distrib/meshlabserver"
MLABSCRIPT="/MeshTracker/resource/scripts/meshlab_decimate.mlx"

# To actually use meshlabserver you need to use a virtual framebuffer and export the display 
# Xvfb :100 &
# export DISPLAY=:100.0

first_frame=00101
last_frame=00886
ext="obj"
length=$(expr $last_frame - $first_frame)
file_pref="Mesh-F"
file_postf="."${ext}

for frame in $(seq -w $first_frame $last_frame)
do

    fname=${file_pref}${frame}${file_postf}
    $MSERVER -i $fname -o $fname -s $MLABSCRIPT

done