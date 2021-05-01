#!/bin/bash

# set meshlabserver dir and script location
MSERVER="/meshlab/distrib/meshlabserver"
MLABSCRIPT_POISSON="/MeshTracker/resource/scripts/meshlab_poisson.mlx"
MLABSCRIPT_ISOTROPIC="/MeshTracker/resource/scripts/meshlab_isotropic_remesh.mlx"

# set instant meshes dir
IMESH="/MeshTracker/resource/thirdparty/Instant_Meshes"

# perform poisson mesh
xvfb-run -a -s "-screen 0 800x600x24" $MSERVER -i $1 -o /tmp/remeshed_graph.ply -s $MLABSCRIPT_POISSON

# perform uniform remesh
# $MSERVER -i /tmp/remeshed_graph.ply -o $2 -s $MLABSCRIPT_ISOTROPIC
$IMESH -o $2 -i -s $3 -r 6 -p 6 /tmp/remeshed_graph.ply

# rm remeshed_graph.ply
