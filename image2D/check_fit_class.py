# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 16:49:18 2012
@author: b1e5n98jamin
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as nd

class PhantomCreate:
    """
    Creates an image phantom (of the same size as the original image)
    from the manual annotated points.
    Used for testing the optimisation procedure
    """

    
    def __init__(self, dpnt, image):
        self.dpnt = dpnt
        self.image = image
        #create a blank image
        self.image2=np.zeros(self.image.shape)
        
    def create_phantom(self, val1):
        
        for ii in range(len(self.dpnt)-1): #because of  (ii+1)
            #generating pixel values along each line segment
            x, y=self.pixel_values(self.dpnt[ii,:], self.dpnt[ii+1,:])
            
            #assigning those pixel values to a new image
            for jj in range(len(x)):
                self.image2[x[jj], y[jj]]=val1
                
    def smooth(self, val1=5):
        self.image2=nd.gaussian_filter(self.image2,val1)
            
    
    def plot(self):
        plt.imshow(self.image2, cmap=plt.cm.gray)

    
    @staticmethod
    def pixel_values(p1, p2):
        d=(((p2-p1)**2).sum())**0.5
        d=np.floor(d+0.5)
        xrange=np.linspace(p1[0], p2[0], num=d)
        yrange=np.linspace(p1[1], p2[1], num=d)
        xrange=np.floor(xrange+0.5) #all positive values so rounding to nearest integer
        yrange=np.floor(yrange+0.5)
        return xrange, yrange
        
        

class CheckFit:
    '''
    Compare the optimisation points to the manual annotation points
    '''
    
    def __init__(self, optpoints, manpoints):

        #initialise using the generated points and the fname being used
        self.optpoints=optpoints
        self.manpoints=manpoints

    def __call__(self):

        #returns the output measure
        print self.dist   
        
    def find_match(self):
        '''
        Finds the closest point on the linear interpolated manual point set
        '''
        
        self.dist=np.zeros((len(self.optpoints),1))
        self.perp=np.zeros((len(self.optpoints),2))
        # for each optimised silhouette point
        for ii in range(len(self.optpoints)): 
            dd=np.empty((len(self.manpoints)-1,1))
            vv=np.empty((len(self.manpoints)-1,2))
            
            # for each line segment
            for jj in range(len(self.manpoints)-1): 
                dd[jj], vv[jj,:]=CheckFit.point_to_line_segment(self.optpoints[ii,:], 
                                        self.manpoints[jj,:], self.manpoints[jj+1,:])                                              
            #closest point
            self.dist[ii]=dd.min()
            # stores the vector direction to the point on the closest line segment
            self.perp[ii,:]=vv[dd.argmin(),:] 
            #measure of fit (mean distance)
            self.dmeas=self.dist.sum()/self.dist.shape[0]
            #hausdorff distance (max distance)
            self.dmeasmax=self.dist.max()
        
        print "Fit: ", self.dmeas, self.dmeasmax
        
    def plot(self, Im, figure_num, type1='all', extra1=0):

        #Visualises the class
        #plot both sets of points with an arrow going between the points

        xaxis1=[(self.manpoints[:,0]).min(), (self.manpoints[:,0]).max()]
        yaxis1=[(self.manpoints[:,1]).min(), (self.manpoints[:,1]).max()]   
            
        plt.figure(figure_num)
        if Im!=0:
            plt.imshow(Im.I[xaxis1[0]:xaxis1[1], yaxis1[0]:yaxis1[1]], cmap=plt.cm.gray)
#        plt.ylim(ymin, ymax)
#        plt.xlim(xmin, xmax)
        
        
        if type1=='all':
            plt.plot(self.optpoints[:,1]-yaxis1[0], self.optpoints[:,0]-xaxis1[0], 'yo')
            plt.plot(self.manpoints[:,1]-yaxis1[0], self.manpoints[:,0]-xaxis1[0], 'r.')
            plt.plot(self.manpoints[:,1]-yaxis1[0], self.manpoints[:,0]-xaxis1[0], 'r')
            #plt.quiver(self.optpoints[:,1]-yaxis1[0], self.optpoints[:,0]-xaxis1[0], 
            #           self.perp[:,1], self.perp[:,0], angles='xy', scale=300)
            
        if type1=='comp':
            plt.plot(extra1[:,1]-yaxis1[0], extra1[:,0]-xaxis1[0], 'yo')
            plt.plot(self.optpoints[:,1]-yaxis1[0], self.optpoints[:,0]-xaxis1[0], 'go')
            plt.plot(self.manpoints[:,1]-yaxis1[0], self.manpoints[:,0]-xaxis1[0], 'r')
            
            
        plt.axis('image') # tight axis
        plt.axis('off')

        
    @staticmethod
    def point_to_line_segment(pnt, a, b):
        """
        Find the closest distance between a point and a line segment
        ref: mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
        http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
        
        pnt: point
        a and b: beginning and end of line
        """
        #For line y = a + u(b-a)
        #point of intersection where dot product is 0
        u=((pnt[0] - a[0])*(b[0]-a[0]) + (pnt[1]-a[1])*(b[1]-a[1]))/((b-a)**2).sum()
        
        if u >1:
            v=b-pnt
        elif u<0:
            v=a-pnt        
        else:
            pline=a+u*(b-a)
            v=pline-pnt
        d=np.sqrt((v**2).sum())
        
        return d, v
        # returns the distance and the vector orientation to the closest point on each line segment
                
def load_ij_text(fname):
    
    '''
    A small script to load imagej manual annotations
    
    Read a text file to a ndarray
    '''    
    with open(fname) as f:
        content=f.readlines()
    
    #Convert content to an integer array
    cc=0;
    pos_array=np.zeros((len(content),2))
    for ii in content:
        ii=ii.strip('\n')
        ii=ii.split('\t')
        # flipped around but otherwise the same
        pos_array[cc,:]=[int(ii[1]), int(ii[0])]
        cc += 1
        
    return pos_array
    