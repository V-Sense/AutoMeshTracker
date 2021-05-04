#!/bin/bash

MSERVER="/meshlab/distrib/meshlabserver"
MLABSCRIPT="/MeshTracker/resource/scripts/meshlab_decimate.mlx"

first_frame=00101
last_frame=00886
ext="obj"
length=$(expr $last_frame - $first_frame)
file_pref="Mesh-F"
file_postf="."${ext}

for frame in $(seq -w $first_frame $last_frame)
do

    fname=${file_pref}${frame}${file_postf}
    xvfb-run -a -s "-screen 0 800x600x24" $MSERVER -i $fname -o $fname -s $MLABSCRIPT

done
