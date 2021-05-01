# MeshTracker

A segmentation-based tracking algorithm for registering volumetric video meshes (ply/obj) in C++.
This is the official implementation of the paper: 
>Moynihan, M., Ruano, S., Pagés, R. and Smolic, A., 2021. Autonomous Tracking For Volumetric Video Sequences. In Proceedings of the IEEE/CVF Winter Conference on Applications of Computer Vision (pp. 1660-1669).

![](/resource/main.jpg)
[![Watch the video](https://img.youtube.com/vi/JwO2obk0tJM/maxresdefault.jpg)](https://youtu.be/JwO2obk0tJM)

## Getting Started

* Git clone the repository into your workspace and see instructions below for install:

## Installing

1. Install Docker

1. Download this repo and run the following from the same location as Dockerfile. 
NUM_JOBS is just an integer which will tell Make how many threads to use (Default=3). 

```
cd dockers && sh ./build_dockers.sh NUM_JOBS
```

** Use as many jobs as your cpu/ram will allow as the docker building is quite slow! **

## Usage 

To run the code on your own input simply run the following command 
 ```
 docker run --rm -v <folder containing .plys or .objs>:/input -v <output folder path>:/result --name meshtracker_instance mesh_tracker:latest MeshTracker -r
 ```
N.B! The code currently does not support keyframe detection yet. For custom data you must provide the 
regions.txt files in the same format as we have provided in the sample data

### Configuration

In lieu of a configuration file (coming soon), the following list of parameters can be tweaked to
accommodate various topology characteristics:
In general most default parameters will be in the header files. Some important ones...
* constants.h: kDG_KNN, amount of graph nodes which influence a given point on the mesh. Higher is more stable but less able to model small scale deformations
* cpd_icp.h: kSamplesUpperLimit, kSamplesLowerLimit: Bounds for segment size and high-res alignment trigger. 
* gauss_newton_solver.h: Most of the parameters for the minimization criteria of the GN solver e.g. rigidity/smoothness
* mesh_segmentation.h
* mesh_tracking.h
* mesh_tracking.cpp:344, set this to true if your input meshes are truly awful. It'll remesh them and hopefully add some stability to processing
mesh_tracking.cpp:461, by default the detail synthesis is disabled, see [this issue](). Set true to enable and tweak parames in mesh_trackpering.cpp:GetDetailLayer() to suit the input. 

## Testing 

We provide two sequences to test the system, Floss and Matt_Hello.
In particular Matt_Hello demonstrates a tracked edit on frame 30 where the
missing hand details in the original have been added back in. 

## Limitations

* Automatic keyframe selection needs to be ported from it's original MATLAB implementation and 
is not yet provided as part of the public release. Keyframes must be provided manually. 
* Will likely fail on any self-intersecting meshes
* Can track fast motion but should be avoided. Ideally there should be a good 
amount of overlap between meshes
* Per the usual apology, this is research code and as such it features some coffee-fuelled, late night
hacks and quick fixes. We'll be supporting it over time but PR's and patience are welcome :) 

## Issues
* cloud_processing.cpp: UniformsamplePointCloud(), broken with update to PCL. Thin structures are removed from output, probably due to voxel tree size. Temporarily replaced by random sampling
* matching.cpp:311, temporarily removed alignment criteria to prevent cases where no influencing node could be found 
* Some components of the detail synthesis step will be store in memory and scale with distance between region boundary and keyframe. It's not a high priority as kfs should be frequent enough but it should be fixed
* mesh_trackpering.cpp:GetDetailLayer, is still very experimental and needs to be adjusted on a per-sequence basis

## Preparing keyframes for tracking user-edits

* In blender, use boolean modifier with union to merge the edited objects into one mesh
* Then apply the remesh modifier adjusting the octree depth and smoothness parameters to best suit your output
* Export the result, then load it in meshlab 
* Perform Quadratic edge collapse. Set desired faces to be 2x the final output
* Clean up non-manifold edges and vertices
* Perfom isometric remeshing, 1 iteration, adaptive, world unit 0.5% 
* Load up the original unedited frame and use Meshlab's align tool to realing it
* Save with same name as original keyframe

## Authors

* **Matt Moynihan** - [GitHub](https://github.com/mjkmoynihan), [LinkedIn](https://www.linkedin.com/in/mjkmoynihan/)

* **Peter Gadomski** - [Coherent Point Drift Implementation (CPD-ICP)](https://github.com/gadomski/cpd)

* **Yizhi Tang** - [Optimal Step Non-Rigid Registration Implementation](https://github.com/Tonsty/Non-Rigid-Registar)

* **Kaiwen Guo** - [Robust non-rigid motion tracking and surface reconstruction using l0 regularization](https://www.guokaiwen.com/svr.html)

Please see LICENSE.txt for usage and further attributions. 

## Paper Errata 

Figure 9 (a)-(d) is missing the following reference to the source of the mesh used.
```
@inproceedings{casas20144d,
  title={4d video textures for interactive character appearance},
  author={Casas, Dan and Volino, Marco and Collomosse, John and Hilton, Adrian},
  booktitle={Computer Graphics Forum},
  volume={33},
  number={2},
  pages={371--380},
  year={2014},
  organization={Wiley Online Library}
}
```
