"""
run_pcaproj
Main file for running the projection algorithm
1) import statistical model
2) project statistical model to 2D
3) align statistical model with 2D airway
4) optimise the fit of the statistical components

Benjamin Irving
2013 / 08 / 25

"""
print __doc__

#Python modules
import numpy as np
import time, copy, os, pickle
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.ndimage as nd
import mayavi.mlab as mlab

#My modules
import point_man_class as pmc
import mesh3D_mod.project_module as pm
import mesh3D_mod.mesh_class as mc
from mesh3D_mod import meshplot_module
from importdata.matinteract_module import H5fileData
import image2D.imagefilt2D_module as ft

        
#Input case
fname='B08'
file_to_load=fname + ' -APLo'
reg_weight=0.003
circ_weight=0.6
edge_weight=0.4

    
# Viewing direction
#Lodox source to detector distance 1299mm
# Projecting onto a plane
ViewP=[1000.0, 0.0, 0.0]
p0=[-100.0, 0.0, 0.0]
proj_type="Lodox"

# reading in mesh and statistical model
math5=H5fileData("data_model/meshdata.h5")
scale1=1.0/math5.scale1.mean()

#PCA operations on the mesh
MH=mc.PcaMesh(faces=math5.pytempmesh_faces, vertices=math5.pytempmesh_vertices,
           style="matlab")

# P - Eigenvectors
# V - Contribution of eigenvalues
MH.input_pca(P=math5.P, V=math5.V, points=math5.pyairpoints, skel=math5.pyairskel)

###################################################################################
#3D mesh manipulation and silhouette detection

#rotation and update of vertices
n11=0.05*np.pi
n22=0.05*np.pi

#Moving along 2 mode of variation (One extreme)
MH.reset_mesh()
#
eigenvalues=np.array([0.0, 1.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.])
#
eigenvalues=eigenvalues.reshape((eigenvalues.shape[0], 1))
eigenvalues=eigenvalues*(3*np.sqrt(MH.V[1:12]))
MH.prin_comp_trans(eigenvalues)
MH.scale_mesh(scale1)
MH.rotate_mesh(n11, "xy")
MH.rotate_mesh(n22, "zx")
#Initialisation of edges
MH.unique_edges()
#Identification of edges
MH.silhouette_vertices("Lodox", ViewP)
#projecting alignment points
points_proj=pm.project_points_ss(MH.silvert1, l0=ViewP, p0=p0) 
#PLOT 1: Plotting the projection (deformed statistical model)
meshplot_module.plot_big(MH, [], points_proj, [], [], ViewP, proj_type, 
	plot_type="simple", fnum=1)

#Moving along 2nd mode of variation (Other extreme)
MH.reset_mesh()

#eigenvalues
eigenvalues=np.array([0.0, -1.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.])
#
eigenvalues=eigenvalues.reshape((eigenvalues.shape[0], 1))
eigenvalues=eigenvalues*(3*np.sqrt(MH.V[1:12]))
MH.prin_comp_trans(eigenvalues)
MH.scale_mesh(scale1)
MH.rotate_mesh(n11, "xy")
MH.rotate_mesh(n22, "zx")
#Initialisation of edges
MH.unique_edges()
#Identification of edges
MH.silhouette_vertices("Lodox", ViewP)
#projecting alignment points
points_proj=pm.project_points_ss(MH.silvert1, l0=ViewP, p0=p0) 
#PLOT 1: Plotting the projection (deformed statistical model)
meshplot_module.plot_big(MH, [], points_proj, [], [], ViewP, proj_type, 
	plot_type="simple", fnum=2)


mlab.show()

