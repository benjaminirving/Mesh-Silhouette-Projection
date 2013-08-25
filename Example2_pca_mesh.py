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

#My modules
import Optimisation_Fitting as ofit
import point_man_class as pmc
import mesh3D_mod.project_module as pm
import mesh3D_mod.mesh_class as mc
from mesh3D_mod import meshplot_module
from importdata.matinteract_module import H5fileData
import image2D.imagefilt2D_module as ft
import mayavi.mlab as mlab
        
#Input case
fname='B08'
file_to_load=fname + ' -APLo'
reg_weight=0.003
circ_weight=0.6
edge_weight=0.4

#Folders to load
filedir1='/home/benjamin/Dropbox/PHD/CAD/CADdata/LodoxTBanon/dicoms/'
if not os.path.isdir(filedir1):
    print "Desktop folders"
    filedir1='/home/b1e5n98jamin/Dropbox/PHD/CAD/CADdata/LodoxTBanon/dicoms/'
    
# Viewing direction
#Lodox source to detector distance 1299mm
# Projecting onto a plane
ViewP=[1000.0, 0.0, 0.0]
p0=[-100.0, 0.0, 0.0]
proj_type="Lodox"#scipy.nd.filters.convolve

# reading in mesh and statistical model
math5=H5fileData("data_model/meshdata.h5")
scale1=1.0/math5.scale1.mean()

#Full example mesh
MHexam=mc.Mesh(faces=math5.pyexammesh_faces, vertices=math5.pyexammesh_vertices,
               style="matlab")
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

#3D mesh
MH.scale_mesh(scale1)
MHexam.scale_mesh(scale1)
MH.rotate_mesh(n11, "xy")
MHexam.rotate_mesh(n11, "xy")    
MH.rotate_mesh(n22, "zx")
MHexam.rotate_mesh(n22,"zx")
#Initialisation of edges
MH.unique_edges()
#Identification of edges
MH.silhouette_vertices("Lodox", ViewP)
#projecting alignment points
points_bif3Dproj=pm.project_points_ss(MH.points, l0=ViewP, p0=p0) 
skel_3Dproj=pm.project_points_ss(MH.skel, l0=ViewP, p0=p0)

####################################################################################
#2D image processing

#X-ray processing
Im=ft.XrayImage(filedir1, file_to_load)
#Cropping black border
Im.crop_border()
Im.select_points(1)
Im2=copy.deepcopy(Im)

#Localised normalisation
Im.doub_norm()
Im.I=ft.local_norm(Im.I, 20)
Im3=copy.deepcopy(Im)
Im.I=nd.gaussian_filter(Im.I,5)
Im.I=-Im.Igrad_mag


Im2.plot()
plt.show()
Im3.plot()
plt.show()

#######################################################################
#Points of interest for manipulation

#manually chosen bifurcation points
points_bif2D=(np.array([Im.xx, Im.yy])).transpose()
# (z, y) representation in 3D space becomes (x,y) representation in 2D
   
Pnt=pmc.PointMan(points_bif3Dproj[:,(2,1)], points_bif2D)
Pnt.get_skel(skel_3Dproj[:,(2,1)]) #inputting skeleton

#####################################################################
#creation of test image (not currently used)
mg=open(filedir1 + file_to_load + "_edge.pkl", 'rb')
man_pts=pickle.load(mg)
mg.close()
from image2D.check_fit_class import PhantomCreate
pc=PhantomCreate(man_pts, Im.I)
pc.create_phantom(-2)
pc.smooth(10)
#pc.plot()
#plt.show()

#--------------------------------------------------------------------
#####################################################################
# Projecting and mapping points
Er1=Pnt.similarity_transform()
Pnt.trans_skel() #transforming the skeleton to the aligned points
points_proj=pm.project_points_ss(MH.silvert1, l0=ViewP, p0=p0) 
Pnt.get_p3D(points_proj[:,(2,1)])
#map 3d points onto local coordinates
Pnt.trans_p3D()
Pnt.plot_landmarks(image=Im3.I, figurenum=1)

######################################################################
# Outward direction of each landmark point
Pnt.get_directions()
Pnt.plot_get_dir()
        
######################################################################
#Initialising the energy function
print "Error from Similariy transform: ", Er1
En=ofit.EnergyFun(Im.I, Im3.I, MH, Pnt, scale1, ViewP, p0, n11, n22, 
                  w_reg=reg_weight, w_circ=circ_weight, w_edge=edge_weight)
# number of eigenvectors/eigenvalues used is defined by the length of e
eigenvalues=np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])  
# number of eigenvalues and regularisation needs to be optimised


#######################################################################
# OPTIMISATION 
print "Energy before optimisation :", En(eigenvalues)
tic=time.clock()#scipy.nd.filters.convolve
eigenvalues_opt=opt.fmin(En, eigenvalues, maxiter=1000)
#eigenvalues_opt=opt.fmin_bfgs(En, eigenvalues, maxiter=20000)

print "Optimisation time: %.5f" % (time.clock()-tic)

#Im=Im2

########################################################################
# Figures and Evaluation (initial input)
print "Standard plot eigenvalues"
EE=En(eigenvalues)
MH, Pnt, points_proj = En.out()
Pnt_init=copy.deepcopy(Pnt) # making a copy for later
# gets the direction of each silhouette point relative to the bif points

#PLOT 1: Plotting
meshplot_module.plot_big(MH, MHexam, points_proj, points_bif3Dproj, skel_3Dproj, ViewP, proj_type, fnum=1)
#mlab.savefig("test1.x3d") #exporting as x3d
#mlab.savefig("figures_thesis/mean_model.png", magnification=2)
#mlab.show()

#PLOT 2: Plotting in 2D
#meshplot_module.plot_im(Im, Pnt)
#plt.show() #plt.ion() allows interactive mode. For updating plot.

#Comparing the method to the manual annotations 
mg=open(filedir1 + file_to_load + "_edge.pkl", 'rb')
man_pts=pickle.load(mg)
mg.close()

#Alternative: Comparing the method to the manual annotations from imagej
#Comment out one of these two inputs
#from image2D.check_fit_class import load_ij_text
#man_pts=load_ij_text(filedir1 + file_to_load + "_edge_reader2.txt")
###

from image2D.check_fit_class import CheckFit
chFit=CheckFit(Pnt.point2d, man_pts)
chFit.find_match()
#PLOT 3: comparing the plot to manual annotations
chFit.plot(Im3, 2)
#plt.savefig(fname +"a.pdf")
#chFit.plot(Im, 5,type1='notall')
#plt.savefig(fname +"c.pdf")

##########################################################################
#Figures and Evaluation (after optimisation)
#rerunning with best parameters
print "Optimal eigenvalues:"
print eigenvalues_opt
EE=En(eigenvalues_opt)
MH, Pnt, points_proj = En.out()

#plotting
#plt.ion()

#PLOT 1: Plotting the projection (deformed statistical model)
meshplot_module.plot_big(MH, [], points_proj, points_bif3Dproj, [], ViewP, proj_type, plot_type="simple", fnum=2)
#mlab.savefig("test1.x3d")
#mlab.savefig("figures_thesis/mean_model2.png", magnification=2)
mlab.show()

#PLOT 2: Plotting in 2D (deformed statistical model)
#meshplot_module.plot_im(Im, Pnt)
#plt.show() #plt.ion() allows interactive mode. For updating plot.

#Comparing to manual segmentation
chFit=CheckFit(Pnt.point2d, man_pts)
chFit.find_match()
chFit.dmeas
#PLOT 3: compare to manual annotations
chFit.plot(Im, 3)
chFit.plot(Im, 4, type1='comp', extra1=Pnt_init.point2d)
plt.savefig("figures_thesis/"+fname +"comp.pdf")

#chFit.plot(Im, 4,type1='notall')
#plt.savefig("fname +"d.pdf")

plt.show()
#Saving Pnt for reuse
#output=open('dat_pts.pkl', 'wb')
#pickle.dump(Pnt, output)
#output.close()


