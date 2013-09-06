Mesh and Principal Component Silhouette Projection
==========================

This is the demo code for the following paper, 

Benjamin Irving, Tania Douglas, Paul Taylor. 2D X-ray airway tree segmentation by 
3D deformable model projection and registration

To be presented at the Fifth International Workshop on Pulmonary Image Analysis (www.lungworkshop.org)


***Creation of examples, instructions and cleaning up of the code is still in progress***

***Dependencies***

Python libraries:
- time, copy, os
- pickle
- matplotlib
- numpy
- scipy
- mayavi


***Installation***

The c++ component needs to first be compiled for your os and wrapped for python using swig. 
(Only tested on linux so far)

- Requirements: swig

Using linux:

``` bash
cd mesh3D_mod
swig -c++ -python -o MeshProject_wrap.cpp MeshProject.i
gcc -fPIC $(python-config --includes) -c MeshProject_wrap.cpp MeshProject.cpp
g++ -shared MeshProject_wrap.o MeshProject.o -o _MeshProject.so
```
This c++ object is now callable from the python 3D processing class. 

***Running Examples***
``` bash
python Example1_mesh_silhouette_and_projection.py
```

``` bash
python Example2_pca_mesh.py
```

``` bash
python Example3_centreline_landmark_projection.py
```

...
***Example Outputs***
