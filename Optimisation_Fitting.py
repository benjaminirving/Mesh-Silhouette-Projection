# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 12:25:28 2012

@author: benjamin
"""
import mesh3D_mod.project_module as pm
import scipy.ndimage as nd
import numpy as np
import matplotlib.pyplot as plt


# TODO
#    Slighty more global search?
#    Edge finding perp to points
#    Energy using intensity
#    Regularisation
#    Bias towards main modes?
#Energy optimisation (with regularisation, intensity mean)
#(once final soln is reached then project entire mesh onto 2d)  

        #Translation component

class EnergyFun:

    def __init__(self, gradI, imI, MH, Pnt, scale1, ViewP, p0, n11, n22,
                  w_reg=0, w_circ=1, w_edge=1, sradius=7, sdist=9):
        
        self.gradI=gradI
        self.imI=imI
        self.MH=MH
        self.Pnt=Pnt
        self.scale1=scale1
        self.ViewP=ViewP
        self.p0=p0
        self.n11=n11
        self.n22=n22
        self.w_circ=w_circ
        self.w_reg=w_reg
        self.w_edge=w_edge
        self.sradius=sradius
        self.sdist=sdist
        
        #generate a standard sphere using the radius
        self.circlelist=self.sphere_coords(self.sradius)
        # getting directional components

        
    
    def __call__(self, eigenvalues): # class is made callable directly for optimisation algorithm
        
        self.eigenvalues=eigenvalues.reshape((eigenvalues.shape[0], 1))
        self.MH.reset_mesh()
        self.MH.prin_comp_trans(self.eigenvalues)
        self.MH.scale_mesh(self.scale1)
        self.MH.rotate_mesh(self.n11, "xy")   
        self.MH.rotate_mesh(self.n22, "zx")
        
        #Projection of points from 3D to 2D
        self.points_proj=pm.project_points_ss(self.MH.silvert1, l0=self.ViewP, p0=self.p0) 
        self.Pnt.get_p3D(self.points_proj[:,(2,1)])
        #map 3d points onto local coordinates
        self.Pnt.trans_p3D()
        #points on the image
        self.xpoints=self.Pnt.point2d        
        
        En=self.E_inten()
        return En
        
    def E_inten(self):
        
        """
        Energy of the system
        Energy based on greyscale value (active contours)
        Using bilinear interpolation to calculate energy at each point     
        """
        
        self.__E_edge()
        self.__E_reg()
        self.__E_circle()

        #cp=np.int16(self.x.round())
        #En=sum(self.gradI[cp[:,0], cp[:,1]])
        
        #print En_edge, w*En_reg
        #En=self.En_edge + self.w*self.En_reg
        En=self.w_edge*self.En_edge + self.w_circ*self.En_circ + self.w_reg*self.En_reg
        #print self.w_edge*self.En_edge, " ", self.w_circ*self.En_circ, " ", self.w_reg*self.En_reg
        return En
    
    def __E_edge(self):
        '''
        Energy from edge term
        '''
        # saving time by not running the alg if weight is 0
        if self.w_edge == 0:
            self.En_edge=0
            return
            
        
        cp=self.xpoints.transpose()
        #finds the grayscale value at each point 
        #on the line using bilinear interpolation of the greyscale value
        #achieving a subvoxel greyscale accuracy for each landmark point
        ec=nd.map_coordinates(self.gradI, cp, order=1)
        self.En_edge=sum(ec)

    
    def __E_reg(self):
        '''
        Regularistation energy
        '''
        if self.w_reg==0:
            self.En_reg=0
            return
        
        
        ll=self.eigenvalues.size
        lambda1 = self.MH.V[0:ll,0]
        # should maybe be sqrt(lambda)
        reg1=self.eigenvalues.transpose()/lambda1
        self.En_reg=np.absolute(reg1).sum()
        
    
    def __E_circle(self):
        '''
        Energy based inward and outward circular sampling
        '''
        if self.w_circ == 0:
            self.En_circ=0
            return

        #finding circle of interest for each point
        ecirc=np.zeros(len(self.xpoints))
        for ii in range(len(self.xpoints)):
            #circle outside airway
            circ1=self.circlelist + self.xpoints[ii,:] + self.sdist*self.Pnt.dir1[ii,:]
            #circle inside airway
            circ2=self.circlelist + self.xpoints[ii,:] - self.sdist*self.Pnt.dir1[ii,:]
            

            
            cp1=circ1.transpose()
            ec1=nd.map_coordinates(self.imI, cp1, order=1)
            cp2=circ2.transpose()
            ec2=nd.map_coordinates(self.imI, cp2, order=1)
            
            mean1=ec1.mean()
            mean2=ec2.mean()
            #mean1=self.imI[circ1b[:,0], circ1b[:,1]].mean()
            #mean2=self.imI[circ2b[:,0], circ2b[:,1]].mean()
            ecirc[ii]=(mean2-mean1)/mean1
            #print ecirc[ii]
            
            #-----------------FOR PLOTTING CIRCLES ------------------
            #circ1b=circ1.astype("int")
            #circ2b=circ2.astype("int")
            #self.imI[circ1b[:,0], circ1b[:,1]]=1
            #self.imI[circ2b[:,0], circ2b[:,1]]=1
            #--------------------------------------------------------
        
        self.En_circ=ecirc.sum()
        #print self.xpoints
        
        #-----------------FOR PLOTTING CIRCLES-----------------------
        #plt.figure(13)
        #plt.imshow(self.imI, cmap=plt.cm.gray)
        #plt.show()
        #------------------------------------------------------------
            
            
        
    def out(self):
        
        return self.MH, self.Pnt, self.points_proj
    
    @staticmethod
    def sphere_coords(radius):
        '''
        Finds the coordinates for all pixels in a circle centred at zero
        Only has to be run once and can then just be translated
        :param radius:
        '''
        circlelist=[]
        for ii in range(-radius, radius):
            for jj in range(-radius, radius):
                if (ii**2 +jj**2) <= radius**2:
                    circlelist.append([ii,jj])
        return circlelist
