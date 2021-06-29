from itertools import groupby
import itertools
from operator import itemgetter
#from trackplot import *
import random 
from scipy import ndimage
#import gaussfitter as gf
import numpy as np
from src.bgo import *
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon
#from collections import Counter

#import time

#from functools import reduce
##############################################################################################################
def det(a, b):
    return a[0] * b[1] - a[1] * b[0]
##############################################################################################################
def line_2d(x,p):
    y = -p[0]/p[1]*x - p[2]/p[1]    
    return y  
##############################################################################################################    
def determine_vertex(line_params, bgo_max):
    #print "VERTEX RECONSTRUCTION: "
    
    #print line_params
    vertex_points = []
    
    for p in line_params:
        for p2 in line_params:
        
            if p == p2:
                continue
               
            x01 = 10
            x02 = -20                    
            y01 = line_2d(x01,p)
            y02 = line_2d(x02,p) 

            x11 = 10
            x12 = -20                    
            y11 = line_2d(x11,p2)
            y12 = line_2d(x12,p2) 
            
            line1 = [[x01,y01],[x02,y02]]
            line2 = [[x11,y11],[x12,y12]]
            
            #print "line1", line1
            #print "line2", line2
            
            xdiff = (x01 - x02, x11 - x12)
            ydiff = (y01 - y02, y11 - y12) 

            
            div = det(xdiff, ydiff)
            
            if div == 0: # no intersection point found
               return vertex_points  

            d = (det(*line1), det(*line2))
            x = det(d, xdiff) / div
            y = det(d, ydiff) / div
            
            #print "POINTS: ", x, y
            vertex_points.append([x,y])
            #return x, y
    vertex = calc_arithm_mean(np.array(vertex_points), bgo_max)
    #print "vertex:    ", vertex
        
    return vertex_points, vertex
##############################################################################################################   
##############################################################################################################
def determine_vertex_tpx(line_params, bgo_max):
    #print "VERTEX RECONSTRUCTION: "
    
    #print line_params
    vertex_points = []
    
    for p in line_params:
        for p2 in line_params:
        
            if p == p2:
                continue
               
            x01 = 10
            x02 = -20                    
            y01 = line_2d(x01,p)
            y02 = line_2d(x02,p) 

            x11 = 10
            x12 = -20                    
            y11 = line_2d(x11,p2)
            y12 = line_2d(x12,p2) 
            
            line1 = [[x01,y01],[x02,y02]]
            line2 = [[x11,y11],[x12,y12]]
            
            #print "line1", line1
            #print "line2", line2
            
            xdiff = (x01 - x02, x11 - x12)
            ydiff = (y01 - y02, y11 - y12) 

            
            div = det(xdiff, ydiff)
            
            if div == 0: # no intersection point found
               return vertex_points  

            d = (det(*line1), det(*line2))
            x = det(d, xdiff) / div
            y = det(d, ydiff) / div
            
            #print "POINTS: ", x, y
            vertex_points.append([x,y])
            #return x, y
    vertex = calc_arithm_mean_tpx(np.array(vertex_points), bgo_max)
    #print "vertex:    ", vertex
        
    return vertex_points, vertex
############################################################################################################## 
##############################################################################################################   
##############################################################################################################    
def calc_arithm_mean(vertex_points, bgo_max):
    #print "ARITHM MEAN"
    bound_minx = bgo_max[0]- 30
    bound_maxx = bgo_max[0]+ 30
    
    bound_miny = bgo_max[1]- 30
    bound_maxy = bgo_max[1]+ 30
    #bgo_bound_x = 50
    #bgo_bound_y = 50
    
    vertex_points = vertex_points[np.logical_and(vertex_points[:,0] < bound_maxx, vertex_points[:,0] >  bound_minx)]
    vertex_points = vertex_points[np.logical_and(vertex_points[:,1] < bound_maxy, vertex_points[:,1] >  bound_miny)]
    
    #print vertex_points
    
    nr_of_vert = len(vertex_points)
    vertex = np.zeros(2)
    vertex[0] = np.sum(vertex_points[:,0])/nr_of_vert
    vertex[1] = np.sum(vertex_points[:,1])/nr_of_vert
    return vertex  
##############################################################################################################    
def calc_arithm_mean_tpx(vertex_points, bgo_max):


    #vertex_points = vertex_points[np.logical_and(vertex_points[:,0] < bound_maxx, vertex_points[:,0] >  bound_minx)]
    #vertex_points = vertex_points[np.logical_and(vertex_points[:,1] < bound_maxy, vertex_points[:,1] >  bound_miny)]
    
    #print vertex_points
    
    nr_of_vert = len(vertex_points)
    
    print "Nr of Vertices", nr_of_vert
    
    vertex = np.zeros(2)
    vertex[0] = np.sum(vertex_points[:,0])/nr_of_vert
    vertex[1] = np.sum(vertex_points[:,1])/nr_of_vert
    
    print vertex
    
    return vertex  
            
##############################################################################################################
########################################################################################################    
def do_2D_gaussian_fit():
        print "copied from main script to save it... not functional yet, needs function parameters set"
        # do a fit of the hit position in BGO with a 2D gaussian

        """fit = gf.gaussfit(orig)
        data = gf.twodgaussian(fit, shape=orig.shape)
        xfit = coords[fit[2], fit[3], 0]
        yfit = coords[fit[2], fit[3], 1]
        print "center x", xfit
        print "center y", yfit"""
       
        # this is just the highest peak and a circle around
        #circle = Circle(xy=(xm,ym), radius=50.0,  color = 'white', lw=2, alpha = 0.5)
        
        
        #circle = Circle(xy=(xm,ym), radius=25.0,  color = 'white', lw=2, alpha = 0.5)
        #ellipse = Ellipse
        #(xy, width, height, angle=0.0, theta1=0.0, theta2=360.0, )
          
        """fit = gf.gaussfit(orig)        
        print "gaussfit params: ***********************"
        #inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)
        #(b is background height, a is peak amplitude)
        print fit
                
        data = gf.twodgaussian(fit, shape=orig.shape)
        
        xfit = coords[fit[2], fit[3], 0]
        yfit = coords[fit[2], fit[3], 1]
        print "center x", yfit
        print "center y", xfit
        #ypos= fit[3]               
        #xpos = fit[2]"""
    
    
    
########################################################################################################    
    
def get_list_of_polygons(max_cluster, coords, bgo_hit_bounds):

        coord_mask = max_cluster[1]
        

        coords_test = coords[coord_mask]
        
        #this is how you get the coordinates of the ports
        portAcoord = coords[0: 8,8:16]
        portBcoord = coords[0: 8,0: 8]
        portCcoord = coords[8:16,8:16] 
        portDcoord = coords[8:16,0: 8] 
        
        coord_maskA = coord_mask[0: 8,8:16]
        coord_maskB = coord_mask[0: 8,0: 8]
        coord_maskC = coord_mask[8:16,8:16] 
        coord_maskD = coord_mask[8:16,0: 8]
       
       # get the rects now 
        hpA = BGOPortPainter('A', coolwarm)
        hpB = BGOPortPainter('B', coolwarm)
        hpC = BGOPortPainter('C', coolwarm)
        hpD = BGOPortPainter('D', coolwarm)
        
        #ind = ndarray.tolist(coord_maskA.reshape((64,1)))
        #ind = [item for sublist in ind for item in sublist]
        ind = coord_maskA.reshape((64,))
        #print ind.shape
        #print np.array(hpA.patches).shape
        A_hitrects = np.array(hpA.patches)[ind]
        B_hitrects = np.array(hpB.patches)[coord_maskB.reshape((64,))]
        C_hitrects = np.array(hpC.patches)[coord_maskC.reshape((64,))]
        D_hitrects = np.array(hpD.patches)[coord_maskD.reshape((64,))]
                        
        all_hitrects = np.append(A_hitrects, B_hitrects)
        all_hitrects = np.append(all_hitrects, C_hitrects)
        all_hitrects = np.append(all_hitrects, D_hitrects)
        
        polygon_hitlist = []
        
        
        
        for hitr in all_hitrects:
            #print hitr.middle
            x1 = hitr.get_x()
            y1 = hitr.get_y()
            
            x2 = x1 + hitr.get_width()
            y2 = y1 + hitr.get_height()
            
            if x1 >= x2:
                if x2 <= bgo_hit_bounds[0][0]:
                    bgo_hit_bounds[0][0] = x2     
                if x1 > bgo_hit_bounds[0][1]:
                    bgo_hit_bounds[0][1] = x1
            else:
                if x1 <= bgo_hit_bounds[0][0]:
                    bgo_hit_bounds[0][0] = x1     
                if x2 > bgo_hit_bounds[0][1]:
                    bgo_hit_bounds[0][1] = x2 
                    
            ##################################################        
                    
            if y1 >= y2:
                if y2 <= bgo_hit_bounds[1][0]:
                    bgo_hit_bounds[1][0] = y2     
                if y1 > bgo_hit_bounds[1][1]:
                    bgo_hit_bounds[1][1] = y1
            else:
                if y1 <= bgo_hit_bounds[1][0]:
                    bgo_hit_bounds[1][0] = y1     
                if y2 > bgo_hit_bounds[1][1]:
                    bgo_hit_bounds[1][1] = y2         
                  
            ext = [(x1, y1), (x1, y2), (x2, y1), (x2,y2)]
            
            polygon = Polygon(ext)   
            polygon_hitlist.append(polygon)    
    
        return polygon_hitlist, all_hitrects    
    
    
    
############################################################################################################## 
def get_bgo_multi(bgo_energy_map, bgo_coords): # find hit clusters in the BGO map 
      
        maxhit_test = np.amax(bgo_energy_map)
        coord_mask_test = np.logical_and(bgo_energy_map <= maxhit_test, bgo_energy_map > maxhit_test*0.5)       

        #print "-----------------------------------------------------------------------------------------------------"
        #print coord_mask_test.astype(int)
        #print "-----------------------------------------------------------------------------------------------------"
        #struct2 = ndimage.generate_binary_structure(2, 2)
        struct1 = ndimage.generate_binary_structure(2, 1)
        # dilate mask to get rid of single pixel cluster and connect close clusters
        coord_mask_dilation = ndimage.binary_dilation(coord_mask_test, structure=struct1).astype(coord_mask_test.dtype)
        #print "-----------------------------------------------------------------------------------------------------"
        #print coord_mask_dilation.astype(int)
        #print "-----------------------------------------------------------------------------------------------------"       
        #        
        s = [[1,1,1],[1,1,1],[1,1,1]] # allow also diagonal connections        
        labelled, numfeatures = ndimage.label(coord_mask_dilation, structure = s)
        
        #print labelled
        #print "-----------------------------------------------------------------------------------------------------"
        #print "features: ", numfeatures             
        #print "-------------------------------------------------------------------------------------------------"

        bgo_coords[:,:,0] +=  6.08*0.5
        bgo_coords[:,:,1] +=  6.08*0.5
        
        nr_of_cluster = 0
        
        cluster = []
        
        cluster_points = []
        for n in range(1,numfeatures+1):
            #print "n", n
            f_bool = (labelled == n)
            nr_bars_in_cluster = len(labelled[f_bool])
            #print nr_bars_in_cluster
            if nr_bars_in_cluster <= 6:
                continue
                
            #nr_of_cluster += 1
            
            #print "-------------------------------------------------------------------------------------------------"
            #print f_bool.astype(int)
            #print "-------------------------------------------------------------------------------------------------"
            coords_f = bgo_coords[f_bool]
            bgo_map_f = bgo_energy_map[f_bool] 
                       
            coords_f = np.column_stack((coords_f, bgo_map_f)) # add pixel energy as 3 column
            
            #print coords_f
            cl_hit = [np.average(coords_f[:,0], weights = coords_f[:,2]), np.average(coords_f[:,1], weights = coords_f[:,2])]
            cluster_points.append(cl_hit) # weighted average, weights = pixel energy
            cluster.append([coords_f, f_bool, cl_hit])
        
        nr_of_cluster = len(cluster)    
        #print "------------------------------------------------------------------------------------------"
        #print "number of cluster: ", nr_of_cluster
        #print "------------------------------------------------------------------------------------------"
        
        if nr_of_cluster == 0:
            x,y =  np.unravel_index(bgo_energy_map.argmax(), bgo_energy_map.shape) 
            xm = bgo_coords[x,y,0]
            ym = bgo_coords[x,y,1] 
            return 0, [0.0,0.0], [0.0], [bgo_coords[coord_mask_test], coord_mask_test, [xm, ym]]
            
        distances = []
        if nr_of_cluster > 1:        
            # first, find cluster with largest E deposit - this will be the hit position of the event
            c_e_max = 0.0
            max_cluster = []
            for c in cluster:
                #print c
                #print "c02: "
                e = np.array(c[0])
                
                c_e = np.sum(e[:,2])
                #print c_e 
                
                if c_e > c_e_max:
                    c_e_max = c_e
                    max_cluster = c        
            #print "-------------------------------------------------------------------------------------------------"           
            for n, cl in enumerate(cluster_points):
                for m, cl2 in enumerate(cluster_points):
                    if cl != cl2 and n > m:
                        #print cl[0:2]
                        distances.append(np.linalg.norm(np.array(cl[0:2]) - np.array(cl2[0:2]))) # distance between clusters                        
        else:
            max_cluster = cluster[0]
                   
        return nr_of_cluster, cluster_points, distances, max_cluster

