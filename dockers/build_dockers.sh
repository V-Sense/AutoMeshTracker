#!/bin/bash

# enter the number of jobs to run in parallel for make commands(default=3)
NUM_JOBS=$1

## Build dependencies in stages
docker build --build-arg NUM_JOBS=${NUM_JOBS} -t mesh_tracker:base -f ./dockerfile_base .
docker build --build-arg NUM_JOBS=${NUM_JOBS} -t mesh_tracker:meshlab -f ./dockerfile_meshlab .
docker build --build-arg NUM_JOBS=${NUM_JOBS} -t mesh_tracker:pcl -f ./dockerfile_pcl .
docker build --build-arg NUM_JOBS=${NUM_JOBS} -t mesh_tracker:cgal -f ./dockerfile_cgal .

## Main code and immediatey thirdparty deps contained here
cd ..
docker build --build-arg NUM_JOBS=${NUM_JOBS} -t mesh_tracker:latest -f ./dockers/dockerfile_build .