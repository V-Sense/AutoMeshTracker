#!/bin/bash

# set meshlabserver dir and script location
MSERVER="/meshlab/distrib/meshlabserver"
MLABSCRIPT_POISSON="/MeshTracker/resource/scripts/meshlab_poisson_6.mlx"

# perform poisson remesh
xvfb-run -a -s "-screen 0 800x600x24" $MSERVER -i $1 -o $2 -s $MLABSCRIPT_POISSON

