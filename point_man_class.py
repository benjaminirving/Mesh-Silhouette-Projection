# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 17:34:54 2012

@author: benjamin
"""

import numpy as np
import matplotlib.pyplot as plt

class PointMan:
    """
    Used to store and manipulate projected points and 2D manual annotations:
    
        point3d # projected 3D edge points found from model
        point2d_init
        point2d # 2D edge points aligned
        match3d # 3d bifurcation points for alignment
        match3d_align
        match2d # 2d bifurcation annotations
        
        dir1= outwards direction of each silhouette point
    """
    
    def __init__(self, match3d, match2d):
        
        self.match3d=match3d
        self.match2d=match2d
        
        self.first_time=True # checks if 3d inputs have been called
        self.skel=[0, 0, 0] # empty skel in case no skeleton is input
        
    def get_skel(self, skel):
        '''
        Getting the skeleton as input
        :param skel: skeleton
        '''
        self.skel = skel
    
    def similarity_transform(self):
        """ 
        translation, scaling and rotation of p1 onto p2
        (see Stegman2002bis - A brief introduction to statistical shape analysis
        as well as "kheng000rr - rigid registration - rigid - similarity - affine.pdf")            
        """
        
        p1=self.match3d 
        p2=self.match2d
        
        #Translation (translate both sets of vectors to their centroid)
        m1=p1.mean(axis=0)
        m2=p2.mean(axis=0)
        p1=p1-m1
        p2=p2-m2
        #plt.figure(1)
        #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
        
        #Scale p1 to p2 (Frobenius norm)
        p1s=np.sqrt(sum(sum(p1**2)))
        p2s=np.sqrt(sum(sum(p2**2)))
        self.s=p2s/p1s;
        p1=p1*self.s;
        #plt.figure(2)
        #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
        
        #Rotation (SVD of correlation to find optimal rotatin of p1)
        corr1=np.dot(p1.transpose(), p2)
        u, ss, v=np.linalg.svd(corr1)
        self.R=np.dot(v, u.transpose())
        p1=np.dot(self.R, p1.transpose())
        p1=p1.transpose()    
        #plt.figure(3)
        #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
        
        #Calculate transform T
        m1prime=self.s * np.dot(self.R, m1.transpose())
        m1prime=m1prime.transpose()
        self.T=m2-m1prime
        
        #p1orig2=s*(np.dot(R, self.match3d.transpose())).transpose() + T
        #plt.figure(14)
        #plt.plot(p1orig2[:,0], p1orig2[:,1],'ro', self.match2d[:,0], self.match2d[:,1], 'go')
        #plt.show()
        
        #print "Just testing without rotation"
        #self.R=np.eye(2)
        #rotating the match3D using these newly calculated transforms
        
        self.match3d_align=self.s*(np.dot(self.R, self.match3d.transpose())).transpose() + self.T
        
        Er1=np.sqrt(((self.match2d-self.match3d_align)**2).sum())
        return Er1        
                
    def get_p3D(self, point3d):
        """ Get 3d point input """
        self.point3d=point3d
                        
        
    def trans_p3D(self):
        """ 
        Updates the 2D point positions from the 3D points 
        Uses transforms generated from similarity_transform()        
        """
        try:
            self.point2d=self.s*(np.dot(self.R, self.point3d.transpose())).transpose() + self.T
            
            if self.first_time is True:
                # stores the first point2d as point3d init
                self.point2d_init=self.point2d
                self.first_time=False

        except:
            print "Could not transform. Has similarity_transform been run?"
            
    def trans_skel(self):
        '''
        Aligning the skeleton points. Similar way to trans_p3D.
        '''
        
        try:
            self.skel_align=self.s*(np.dot(self.R, self.skel.transpose())).transpose() + self.T

        except:
            print "Could not transform. Has similarity_transform been run?"
    
    def get_directions(self):
        '''
        Gets the outward pointing direction of each initial silhouette point 
        after 2d alignment
        '''
        self.dir1=np.zeros(self.point2d_init.shape)
        for ii in range(len(self.point2d_init)):
            self.dir1[ii,:]=self.pnt_dir(self.point2d_init[ii,:], self.skel_align)
            
    def plot_landmarks(self, figurenum, image):
        '''
        plot the manually input and aligned landmarks
        '''
        try:
            plt.figure(figurenum)
            ex1=60
            xaxis1=[(self.point2d[:,0]).min()-ex1, (self.point2d[:,0]).max()+ ex1]
            yaxis1=[(self.point2d[:,1]).min() - ex1, (self.point2d[:,1]).max() + ex1]       
            plt.imshow(image[xaxis1[0]:xaxis1[1], yaxis1[0]:yaxis1[1]], cmap=plt.cm.gray)
            
            plt.plot(self.match2d[:,1]-yaxis1[0], self.match2d[:,0]-xaxis1[0], 'go')
            plt.plot(self.match3d_align[:,1]-yaxis1[0], self.match3d_align[:,0]-xaxis1[0], 'rs')
            plt.axis('image')
            plt.show()
        except:
            print "Could not plot. Have landmarks been generated?"
            
    def plot_get_dir(self):
        '''
        plot the directional components
        '''
        plt.figure(1)
        plt.plot(self.point2d_init[:,0], self.point2d_init[:,1], 'go')
        plt.plot(self.match3d_align[:,0], self.match3d_align[:,1], 'ro')
        #adding the skeleton
        plt.plot(self.skel_align[:, 0], self.skel_align[:, 1], 'bo')
        plt.quiver(self.point2d_init[:,0], self.point2d_init[:,1], 
                            self.dir1[:,0], self.dir1[:,1], angles='xy')
        plt.axis('off')
        plt.show()
        
        pass 
        
        
    @staticmethod
    def pnt_dir(pnt, skel):
        
        '''
        Finding the outwards directions for each point
        Based on the four branch points Top trachea, centre, lmb end, rmb end
        1) Identifies closest two structure points to a silhouette point
        2) Finds outwards direction by choosing othogonal direction to centreline
        3) With cross product used to define sign
        '''
         
        from operator import itemgetter #used to select item from which to sort
        
        dist=np.sqrt(((skel-pnt)**2).sum(1))
        dist_tup=list(enumerate(dist))
        dist_tup.sort(key=itemgetter(1), reverse=False) #sorting
        # now use indices 0 and 1 to create the vector
        
        p_use=dist_tup[0][0]
        dir1=pnt-skel[p_use,:]
        dir1=dir1/np.sqrt((dir1**2).sum())
        return dir1
            
        
        #find nearest 2 points
        #cross product of trachea-mid, mid-lim, mid-rmb with point-trachea / point-mid
        #point outwards directions
    
#    def trans_p2D(self):
#        """ Updates the 3D point positions from the 2D points """
#        try:
#            x=1
#        except:
#            print "Could not transform. Do 2d points exist?"
            
#    def move_image(self, Im_edx, Im_edy):
#        hh=1.0; #step size
#        cp=np.int16(self.point2d.round())
#        self.dir1=np.zeros(self.point2d.shape)
#        self.dir1[:,0]=Im_edx[cp[:,0], cp[:,1]]
#        self.dir1[:,1]=Im_edy[cp[:,0], cp[:,1]]
#        
#        dir_n=((self.dir1**2).sum(axis=1))**0.5
#        dir_n=dir_n.reshape(dir_n.shape[0], 1)
#        self.dir1=self.dir1/dir_n
#        
#        self.point2d=self.point2d + self.dir1*hh
        
        
    