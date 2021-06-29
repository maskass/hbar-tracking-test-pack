from itertools import groupby
import itertools
from operator import itemgetter
from trackplot_final import *
import random 

from collections import Counter
#from cycler import cycler

import time

import shapely

#from functools import reduce
from uncertainties import ufloat

from src.hodoscope_trackhelper import *
from src.vertexfit_3D import *

#from ROOT import TF3, Fit, Minuit2
import ROOT

import warnings


parameterisation = 'angle'
#parameterisation = 'p_d'

vertexfit_choice = 'manual' # 'manual' # python, root or manual


def set_list_intersection(set_list):
  if not set_list:
    return set()
  result = set_list[0]
  for s in set_list[1:]:
    result &= s
  return result

##############################################################################################################    
############################################################################################################## 
############################################################################################################## 

def find_consecutive_nr(numbers):    
    clusters = []
    yesno = False
    for k, g in groupby(enumerate(numbers), lambda (i, x): i-x): # total clusters... also single bars
        clusters.append(map(itemgetter(1), g))
    print "clusters::::::::::: ", clusters
    tempc = [] # save clusters with hits of 31 and 0
    nwclusters = list(clusters)  # needed for loop
        
    for c in clusters:
        if len(c) > 1:
            #print "c::::::::::", c
            yesno = True  
        if 31 in c:
            #print "yes inside!"
            tempc.append(c)
            nwclusters.remove(c)
        if 0 in c:
            tempc.append(c)
            nwclusters.remove(c)
    
    #print "clusters temp::::::::::: ", tempc                
    if len(tempc) >= 1: # if both clusters with 31 and 0 exist, merge them
        newcluster = [item for sublist in tempc for item in sublist]
        nwclusters.append(newcluster)
 
    #print "clusters::::::::::: ", nwclusters
    
    
    
    return yesno, nwclusters    # bool: found cluster: yes or no? 
    
##############################################################################################################      
    
def find_consecutive_nr_singletoo(numbers):    
    clusters = []
    yesno = False
    for k, g in groupby(enumerate(numbers), lambda (i, x): i-x): # total clusters... also single bars
        clusters.append(map(itemgetter(1), g))
    #print "clusters::::::::::: ", clusters
    tempc = [] # save clusters with hits of 31 and 0
    nwclusters = list(clusters)  # needed for loop
        
    for c in clusters:
        if len(c) > 1:
            #print "c::::::::::", c
            yesno = True  
        if 31 in c:
            #print "yes inside!"
            tempc.append(c)
            if c in nwclusters:
                nwclusters.remove(c)
        if 0 in c:
            tempc.append(c)
            if c in nwclusters:
                nwclusters.remove(c)
    
    #print "clusters temp::::::::::: ", tempc                
    if len(tempc) >= 1: # if both clusters with 31 and 0 exist, merge them
        newcluster = [item for sublist in tempc for item in sublist]
        nwclusters.append(newcluster)

    return yesno, nwclusters    # bool: found cluster: yes or no?     

##############################################################################################################  
"""def track_finding_circle(pmapI, cmap, currI, currO, circle, BGOXPos, BGOYPos):    # with density of lines added
    ### First: Find Clusters ###
    #print "####################"
    #print "inner bars hit: ", cmap[:nI,1] # corner map, all rects with corner points plus nr
    #print "outer bars hit: ", cmap[nI:,1]     
    #plt.show()
    
    #print "currI ", currI

    nI = len(currI)
    nO = len(currO)
   
    cmapI = cmap[:nI,:]
    cmapO = cmap[nI:,:] 
    
    hI = [] # to get it into the right shape for the big loop below...
    hO = []
    for i in currI:
        hI.append([i])
        
    for o in currO:
        hO.append([o])
    
    # use this for no unification of rects
    yesnoI, cI = find_consecutive_nr(currI)
    yesnoO, cO = find_consecutive_nr(currO)
    
    # find tracks with clusters: merge rects of clusters to a big one
    #yesnoI, cI = find_consecutive_nr(cmap[:nI,1])
    #print "innter: ####################"
    #print cI
    #yesnoO, cO = find_consecutive_nr(cmap[nI:,1])
    #print "outer : ####################"
    #print cO
    #union_cl_I = merge_rects_of_cluster(cmap[:nI], cI) # union contains the rects of the bars
    #union_cl_O = merge_rects_of_cluster(cmap[nI:], cO)   
    #print "union thing "
    #print union_cl_O
    
    union_cl_O = get_rects_of_hits(cmapO)
    union_cl_I = get_rects_of_hits(cmapI)
   
    hitcollection = []   # 2dim. lines : tracklines columns: hit bars
    hitcollection_count = []
    p0 = [0,0]
#    t1 = []
#    t2 = []
    p1 = [0,0]
    #pbgo = [BGOXPos, BGOYPos]
    line_points_collection = []

    
    for m,curr_pos in enumerate(pmapI): # loop over outer bars i[0],i[1] are x,y coord, i[2] is the bar number
            p0[0] = curr_pos[0]
            p0[1] = curr_pos[1]                                    
            #t1, t2 = Get_TangentP(p0)      
            for i in range(0,30):  # draw XX points from the circle / polygons aka all pmt pixel with a hit
                p1 = draw_pt_from_circle(circle.radius,[BGOXPos,BGOYPos])
                #p1 = draw_point_from_bgohit(polygon_hitlist)
                hitbars = []  # contains bars on one line
                           
                for count,(pol,nrI) in enumerate(zip(union_cl_I,hI)):  # for ever line, check if it intersects one of the hit bars (inner and outer)                    
                    yes = Polygon_contains_line(pol,p0,p1)
                    if yes == True:                   
                        if not [0,nrI] in hitbars:
                            hitbars.append([0,nrI])    
                            line_points_collection.append([p0[0],p0[1],p1[0],p1[1]]) #this is only used for drawing           
                        
                for count,(pol,nrO) in enumerate(zip(union_cl_O,hO)):  # for ever line, check if it intersects one of the hit bars (inner and outer)                    
                    yes = Polygon_contains_line(pol,p0,p1)
                    if yes == True:                   
                        if not [1,nrO] in hitbars:
                            hitbars.append([1,nrO])    
                            line_points_collection.append([p0[0],p0[1],p1[0],p1[1]]) #this is only used for drawing                           
                             
                sorted_track = sorttrack(hitbars)
                #print "sorted track ", sorted_track                
                
                
                if len(hitcollection) == 0 and len(sorted_track) > 0:
                    #np.insert(sorted_track, )
                    #hitcollection = np.array([sorted_track, 1])
                    hitcollection.append(sorted_track)
                    hitcollection_count.append(1)

                if not sorted_track in hitcollection and len(sorted_track) != 0:
                   #print "not inside!", sorted_track
                   hitcollection.append(sorted_track)
                   hitcollection_count.append(1)
                   
                elif sorted_track in hitcollection:
                   ind = hitcollection.index(sorted_track)#
                   hitcollection_count[ind] += 1
                   #print "inside!", sorted_track
                
    #sumup = 0               
    all_lines = sum(hitcollection_count)
    #merge_coll = []
    for h,c in zip(hitcollection, hitcollection_count):
        #print h,c
        h.append(float(c)/float(all_lines))
        #sumup += float(c)/float(all_lines)
        #print sumup
        #merge_coll.append([h,c])
        
    #print "sum "   , all_lines 

    hitcollection.sort(key=len,reverse=True)
    
    
    # TODO drawing: print hit collection on the side for presentations 

    
    #print len(line_points_collection)
    hitcollection2 = []
    #print "##################################"
    for h in hitcollection:
        if len(h[0])>1: 
            hitcollection2.append(h)
            #print h
    #print "##################################"    
    
    final_tracks = new_select_tracks_for_2D(hitcollection2,hI, hO, cI, cO)
    
    #plot_all_bars_2D(circle,cornermapI,cornermapO,union_cl_I, union_cl_O, line_points_collection)
    
    
    return final_tracks, line_points_collection"""

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
    
"""def get_list_of_polygons(max_cluster, coords, bgo_hit_bounds):

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
def get_bgo_multi(bgo_energy_map, bgo_coords):
      
        maxhit_test = np.amax(bgo_energy_map)
        coord_mask_test = np.logical_and(bgo_energy_map <= maxhit_test, bgo_energy_map > maxhit_test*0.5)       

        #print "-----------------------------------------------------------------------------------------------------"
        #print coord_mask_test.astype(int)
        #print "-----------------------------------------------------------------------------------------------------"
        #struct2 = ndimage.generate_binary_structure(2, 2)
        struct1 = ndimage.generate_binary_structure(2, 1) # dilate mask to get rid of single pixel cluster and connect close clusters
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
        #print "-----------------------------------------------------------------------------------------------------"

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
        print "------------------------------------------------------------------------------------------"
        print "number of cluster: ", nr_of_cluster
        print "------------------------------------------------------------------------------------------"
        
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
"""

##################################################################################################################################
##################################################################################################################################
# draw random points (x,y) from a circle around 0 with radius r
def points_on_circ(r, n=100):
    return [[math.sin(2*math.pi/n*x)*r, math.cos(2*math.pi/n*x)*r] for x in xrange(0,n+1)]
##################################################################################################################################
##################################################################################################################################
def Calc_xy_error(barnr, layer):

    tot_ang = []
    rep = [1.0]*4
    for i in range(8):
        for j in rep:
            angle = i*2*np.pi/8
            tot_ang.append(angle)
    
    #print "tot ang", tot_ang
    #print len(tot_ang)
    
    if layer == 1:
        b_ = 35.0 # outer
    else:
        b_ = 20.0
        
    sigmau2 = b_*b_/12.0
    sigmav2 = 5.0*5.0/12.0
    
    V1 = np.array([[sigmau2, 0.0],[0.0, sigmav2]])
    
    tot_ang = np.array(tot_ang)
    
    alpha = tot_ang[barnr]
    
    #print "barnr alpha", barnr, alpha
    
    #for i in tot_ang:
    #    alpha = i
    rota = np.array([[math.cos(alpha),-math.sin(alpha)],[math.sin(alpha), math.cos(alpha)]])
    rota_t = np.array([[math.cos(alpha),math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])
    mat = np.dot(rota,np.dot(V1,rota_t))
        #print math.sqrt(mat[0][0]), math.sqrt(mat[1][1])
   
    return  math.sqrt(mat[1][1]), math.sqrt(mat[0][0])  
     
##################################################################################################################################

def get_hit_error(rectlist): # x and y error... determined by geometry 
    xerr = []
    yerr = []
    for rect in rectlist:   
        maxy = rect[np.where(rect[:,1] == np.max(rect[:,1]))][0][1]
        miny = rect[np.where(rect[:,1] == np.min(rect[:,1]))][0][1]
        maxx = rect[np.where(rect[:,0] == np.max(rect[:,0]))][0][0]
        minx = rect[np.where(rect[:,0] == np.min(rect[:,0]))][0][0] 
        
        xerr.append(abs(maxx-minx)*0.5)  
        yerr.append(abs(maxy-miny)*0.5)

    return sum(xerr)/len(xerr), sum(yerr)/len(yerr)
    
##################################################################################################################################
def overlap_interval(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
##################################################################################################################################
########################################################################################################################################################
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
########################################################################################################################################################

# checks if there is overlap between bars and fibres z position and saves the matches (all of them, no removal if already added)
def bar_fibre_matching_2point0(pmapI, pmapO, zposII, zposOO, zerrI, zerrO, currI, currO, z_I_Fibre_hits, z_O_Fibre_hits, curr_fI, curr_fO):
   
    #print z_I_Fibre_hits
    print "curr_fI", curr_fI
    print "curr_fO", curr_fO
    
    #print " "
    
    z_curr_fI = z_I_Fibre_hits[curr_fI] 
    z_curr_fO = z_O_Fibre_hits[curr_fO]  
      
      
    combs_ch_nrs = cartesian((curr_fI, curr_fO))
    #print combs_ch_nrs
    #print " "
    combs_z_pos = cartesian((z_curr_fI, z_curr_fO))
    #print combs_z_pos
    #print " "
    diff_z_pos = combs_z_pos[:,1] - combs_z_pos[:,0]
    #print diff_z_pos
    
    x_i = combs_z_pos[:,0] + diff_z_pos*0.284 # list of z positions where the inner bar hits could be
    x_o = combs_z_pos[:,0] + diff_z_pos*1.484
    
    # create search intervals for hits with the error on the bars
    
    #print "  "
    #print x_i - zerrI
    #print x_i + zerrI

    search_int_i = np.vstack(([x_i - zerrI], [x_i + zerrI])).T
    search_int_o = np.vstack(([x_o - zerrO], [x_o + zerrO])).T
    
    print " "
    print search_int_i
    
    #print " "
    #print "x_i", x_i
    #print "x_o", x_o              
        
    print " "
    print "zposII",zposII
    print "zposOO", zposOO   
    
    print combs_ch_nrs.shape, combs_z_pos.shape, search_int_i.shape, search_int_o.shape
    
    collection_tracks = []
    
    for chnr_, zcom_, sI_, sO_ in zip(combs_ch_nrs, combs_z_pos, search_int_i, search_int_o):   
        track = []
        for innerbarnr, zposibar_, xy_ibar in zip(currI, zposII, pmapI):
            #print track
            track = []
            #print zposibar_, sI_
            if zposibar_ >= sI_[0] and zposibar_ <= sI_[1]:
                #print " in if! "
                ibar_xx = xy_ibar[0]
                ibar_yy = xy_ibar[1]
                err_ibar_xx, err_ibar_yy = Calc_xy_error(innerbarnr, 0) 
                # add inner bar hit
                track.append([zposibar_, zerrI/2.355, ibar_xx, err_ibar_xx, ibar_yy, err_ibar_yy, innerbarnr, 999, 1, 0])                  
                # inner fibre second hit of track
                track.append([zcom_[0], 4.0/math.sqrt(12.0), ibar_xx, err_ibar_xx, ibar_yy, err_ibar_yy, innerbarnr, chnr_[0], 1, 0])  
        
        
                for outerbarnr, zposobar_, xy_obar in zip(currO, zposOO, pmapO):        
                    if zposobar_ >= sO_[0] and zposobar_ <= sO_[1]:
                        obar_xx = xy_obar[0]
                        obar_yy = xy_obar[1]
                        err_obar_xx, err_obar_yy = Calc_xy_error(outerbarnr, 1) 
                        # add outer bar and fibre hit
                        to_add_ob = [zposobar_, zerrO/2.355, obar_xx, err_obar_xx, obar_yy, err_obar_yy, outerbarnr, 999, 1, 1]                        
                        to_add_of = [zcom_[1], 4.0/math.sqrt(12.0), obar_xx, err_obar_xx, obar_yy, err_obar_yy, outerbarnr, chnr_[1], 1, 1 ]
                        
                        track_temp = list(track)
                        track_temp.append(to_add_ob)                                         
                        
                        track_temp.append(to_add_of)   
               
                        collection_tracks.append(track_temp)  
                    
                    """if len(track) == 2: # if only inner bars have been added
                        obar_xx = xy_obar[0]
                        obar_yy = xy_obar[1]
                        err_obar_xx, err_obar_yy = Calc_xy_error(outerbarnr, 1) 
                        # add outer bar and fibre hit
                        to_add_ob = [zposobar_, zerrO/2.355, obar_xx, err_obar_xx, obar_yy, err_obar_yy, outerbarnr, 999, 1, 1]     
                        track_temp = list(track)
                        track_temp.append(to_add_ob)
                        collection_tracks.append(track_temp)   # """
    #fo_p_3D.append([oz[0], oz[1]/math.sqrt(12.0), point[0], errxx, point[1], erryy, oo, fo[0], len(fo), 1])
    
    
    
    
    """for tr in collection_tracks:
        for tt in tr:
            print tt
            
        print " "    #"""
    
    print "leeeeeeeeeeeeeeeeeeeeeen new", len(collection_tracks) 
    return collection_tracks       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                
        #z = np.mean(z)  
         
    """for i_point, zposI in zip(pmapI, zposII):   # loop over bars and their z    
            zerr = zerrI
            #print "zposI", zposI
            bar = int(i_point[2])
            if overlap_interval([z-2.0,z+2.0], [zposI - zerr, zposI + zerr]): # check for overlap of z 
                #print "z: ", z, "in: [", zposI-zerr, ", ", zposI+zerr, "]"
                bars_i_noleft.append(bar)                 
                matching_i.append([bar, j])                              
                                
    matching_find = []    
    if len(matching_i) != 0:
        #matching_find = np.array(matching_i)[:,1]   # get all the fibres that match a hodor bar
        for i in matching_i:
            #print i[1]
            matching_find.append(i[1])

    for i in curr_fI:
       
        if not i in matching_find:
            nomatch_i.append(i)
        
    #print "OUTER FIBRE: --------------------------------------------------"      
    for j in curr_fO:
        z = z_O_Fibre_hits[j]        
        z = np.mean(z)
        for o_point, zposO in zip(pmapO,zposOO):       
            zerr = zerrO   
            bar = int(o_point[2])         
            if overlap_interval([z-2.0,z+2.0], [zposO - zerr, zposO + zerr]):
                #print "z: ", z, "in: [", zposO-zerr, ", ", zposO+zerr, "]"  
                 
                bars_o_noleft.append(bar)            
                matching_o.append([bar, j])

    bars_i_left = list(set(bars_i_noleft).symmetric_difference(set(currI)))
    bars_o_left = list(set(bars_o_noleft).symmetric_difference(set(currO)))
        
    #print "inner match:", matching_i
    #print "inner nomatch: ", nomatch_i


    matching_find = []    
    if len(matching_o) != 0:
        #matching_find = np.array(matching_o)[:,1]   # get all the fibres that match a hodor bar
        for o in matching_o:
            #print i[1]
            matching_find.append(o[1])
    
    for o in curr_fO:
        if not o in matching_find:
            nomatch_o.append(o)

    #print "outer match:", matching_o
    #print "outer nomatch: ", nomatch_o  
        
    return matching_i, matching_o, nomatch_i, nomatch_o, bars_i_left, bars_o_left # """
    
    
########################################################################################################################################################
########################################################################################################################################################
    
# checks if there is overlap between bars and fibres z position and saves the matches (all of them, no removal if already added)
def bar_fibre_matching(pmapI, pmapO, zposII, zposOO, zerrI, zerrO, currI, currO, z_I_Fibre_hits, z_O_Fibre_hits, curr_fI, curr_fO):
   
    matching_i = [] # fibres with a bar match
    matching_o = []
    
    nomatch_i = [] # fibres that dont match any bar
    nomatch_o = []
    
    bars_i_noleft = [] # bars that have a fibre hit (I need this to find bars that DONT have a fibre hit)
    bars_o_noleft = []
    
    #print "INNER FIBRE: --------------------------------------------------"     
    for j in curr_fI: # loop over fibre hits
        z = z_I_Fibre_hits[j]        
        z = np.mean(z)   
        for i_point, zposI in zip(pmapI, zposII):   # loop over bars and their z    
            zerr = zerrI
            #print "zposI", zposI
            bar = int(i_point[2])
            if overlap_interval([z-2.0,z+2.0], [zposI - zerr, zposI + zerr]): # check for overlap of z 
                #print "z: ", z, "in: [", zposI-zerr, ", ", zposI+zerr, "]"
                bars_i_noleft.append(bar)                 
                matching_i.append([bar, j])                              
                                
    matching_find = []    
    if len(matching_i) != 0:
        #matching_find = np.array(matching_i)[:,1]   # get all the fibres that match a hodor bar
        for i in matching_i:
            #print i[1]
            matching_find.append(i[1])

    for i in curr_fI:
       
        if not i in matching_find:
            nomatch_i.append(i)
        
    #print "OUTER FIBRE: --------------------------------------------------"      
    for j in curr_fO:
        z = z_O_Fibre_hits[j]        
        z = np.mean(z)
        for o_point, zposO in zip(pmapO,zposOO):       
            zerr = zerrO   
            bar = int(o_point[2])         
            if overlap_interval([z-2.0,z+2.0], [zposO - zerr, zposO + zerr]):
                #print "z: ", z, "in: [", zposO-zerr, ", ", zposO+zerr, "]"  
                 
                bars_o_noleft.append(bar)            
                matching_o.append([bar, j])

    bars_i_left = list(set(bars_i_noleft).symmetric_difference(set(currI)))
    bars_o_left = list(set(bars_o_noleft).symmetric_difference(set(currO)))
        
    #print "inner match:", matching_i
    #print "inner nomatch: ", nomatch_i


    matching_find = []    
    if len(matching_o) != 0:
        #matching_find = np.array(matching_o)[:,1]   # get all the fibres that match a hodor bar
        for o in matching_o:
            #print i[1]
            matching_find.append(o[1])
    
    for o in curr_fO:
        if not o in matching_find:
            nomatch_o.append(o)

    #print "outer match:", matching_o
    #print "outer nomatch: ", nomatch_o  
        
    return matching_i, matching_o, nomatch_i, nomatch_o, bars_i_left, bars_o_left
    
########################################################################################################################################################
########################################################################################################################################################
def plot_event_3_views(fig, ax, h, list_i, list_o, pmapcI, pmapcO, zposII, zposOO, zerrI, zerrO, width_list_x, width_list_xo, width_list_y, width_list_yo,
                        currI, currO, match_i, match_o, nomatch_i, nomatch_o, if_z, of_z, curr_fI_c, curr_fO_c, circle, union_cl_I, union_cl_O, vertex_BGO, fI_proj, fO_proj):

    # just some line widths for plotting
    bluelw = 1.25
    bluelw2 = 1.25
    bluems = 2
    bluems2 = 1
    blacklw = 0.3
    
    zerrI =zerrI/2.355
    zerrO = zerrO/2.355


    
    # plot the gray outlines of the bar hodoscope -- this is for the side view. top view see further down
    ax[0].set_ylim(-250, 250)
    ax[0].set_xlim(-300, 300)
    
    # FIXME zoomed version:    
    #ax[0].set_ylim(-60, 60)
    #ax[0].set_xlim(-30, 30)
    
    ax[0].set_xlabel("z (mm)")
    ax[0].set_ylabel("y (mm)")
    #fig.suptitle("SIDE/TOP VIEW")
    
    ohodor_side = Polygon([(225.0, -175.0), (225.0, 175.0), (-225.0, 175.0), (-225.0, -175.0)]) 
    ax[0].add_patch(PolygonPatch(ohodor_side, fc='white', ec='grey', lw=1, alpha=1.0)) # 
    ihodor_side = Polygon([(150.0, -100.0), (150.0, 100.0), (-150.0, 100.0), (-150.0, -100.0)]) 
    ax[0].add_patch(PolygonPatch(ihodor_side, fc='white', ec='grey', lw=1, alpha=1.0)) #
    
    # this is the timepix/BGO: 
    bgo_patch = Polygon([(2.5, -45.0), (2.5, 45.0), (-2.5, 45.0), (-2.5, -45.0)]) 
    ax[0].add_patch(PolygonPatch(bgo_patch, fc='black', ec='black', lw=0.2, alpha=0.5)) #
    
    for n, point in enumerate(pmapcI):
            zposI = zposII[n]        
            zerr = zerrI
            # plot the z point
            ax[0].plot(zposI, point[1], linestyle='None', marker='.', c = 'cornflowerblue',linewidth=bluelw, markersize =bluems)
            
            #xerr, yerr = get_hit_error([rec])
            xerr, yerr = Calc_xy_error(list_i[n],0)
            # plot the error bars:
            ax[0].plot([zposI, zposI], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'cornflowerblue', markersize =bluems2)
            ax[0].plot([zposI, zposI], [point[1]+yerr, point[1]-yerr],  marker="_",linewidth=bluelw2, c = 'cornflowerblue', markersize =bluems2)
            ax[0].plot([zposI+zerr, zposI-zerr], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'cornflowerblue',  markersize =bluems2)
            
    # plot the outer hodoscope bar hits ########    
    for n, point in enumerate(pmapcO):
            zposO = zposOO[n]
            #print zposO
            zerr = zerrO
            ax[0].plot(zposO, point[1], linestyle="None", marker=".", c = "crimson",linewidth=bluelw, markersize =bluems)        
            #xerr, yerr = get_hit_error([rec])
            xerr, yerr = Calc_xy_error(list_o[n],1)
            ax[0].plot([zposO, zposO], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)
            ax[0].plot([zposO, zposO], [point[1]+yerr, point[1]-yerr],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)
            ax[0].plot([zposO+zerr, zposO-zerr], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)
     
        
    ax[1].set_ylim(-250, 250)
    ax[1].set_xlim(-300, 300)
    # FIXME zoomed version:    
    #ax[1].set_ylim(-60, 60)
    #ax[1].set_xlim(-30, 30)
    
    ax[1].set_xlabel("z (mm)")
    ax[1].set_ylabel("x (mm)")
    
    # plot the outlines for the hodoscope top view
    ohodor_top = Polygon([(225.0, -175.0), (225.0, 175.0), (-225.0, 175.0), (-225.0, -175.0)]) 
    ax[1].add_patch(PolygonPatch(ohodor_top, fc='white', ec='grey', lw=1, alpha=1.0)) # 
    ihodor_top = Polygon([(150.0, -100.0), (150.0, 100.0), (-150.0, 100.0), (-150.0, -100.0)]) 
    ax[1].add_patch(PolygonPatch(ihodor_top, fc='white', ec='grey', lw=1, alpha=1.0)) # 
    # bgo patch for the top view:
    bgo_patch = Polygon([(2.5, -45.0), (2.5, 45.0), (-2.5, 45.0), (-2.5, -45.0)]) 
    ax[1].add_patch(PolygonPatch(bgo_patch, fc='black', ec='black', lw=0.2, alpha=0.5)) # 
    
    
    # drawing the hodoscope bar hits for the side view
    for n, point in enumerate(pmapcI):
        zposI = zposII[n]        
        zerr = zerrI
        ax[1].plot(zposI, point[0], linestyle="None", marker=".", c = "cornflowerblue",linewidth=bluelw, markersize =bluems)
        #xerr, yerr = get_hit_error([rec])
        xerr, yerr = Calc_xy_error(list_i[n],0)
        ax[1].plot([zposI, zposI], [point[0], point[0]],  marker="_",linewidth=bluelw2, c = 'cornflowerblue', markersize =bluems2)
        ax[1].plot([zposI, zposI], [point[0]+xerr, point[0]-xerr],  marker="_",linewidth=bluelw2, c = 'cornflowerblue', markersize =bluems2)
        ax[1].plot([zposI+zerr, zposI-zerr], [point[0], point[0]],  marker="_",linewidth=bluelw2, c = 'cornflowerblue',  markersize =bluems2)
        
    for n, point in enumerate(pmapcO):
        zposO = zposOO[n]
        #print zposO
        zerr = zerrO
        ax[1].plot(zposO, point[0], linestyle="None", marker=".", c = 'crimson',linewidth=bluelw, markersize =bluems)        
        #xerr, yerr = get_hit_error([rec])
        xerr, yerr = Calc_xy_error(list_o[n],1)
        ax[1].plot([zposO, zposO], [point[0], point[0]],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)
        ax[1].plot([zposO, zposO], [point[0]+yerr, point[0]-yerr],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)
        ax[1].plot([zposO+zerr, zposO-zerr], [point[0], point[0]],  marker="_",linewidth=bluelw2, c = 'crimson', markersize =bluems2)#"""

   
   
    for iz, fi in zip(if_z, curr_fI_c):
        for point, oi in zip(fI_proj, currI):  

            if [oi,fi] in match_i:      
                xerr = width_list_x[oi]*0.5
                yerr = width_list_y[oi]*0.5
                #print "inner xerr, yerr: ", oi, xerr, yerr
                
                # those polis are just for plotting: (I also did some track fitting tests with the 2D projections)
                poli_side = Polygon([(iz[0] + iz[1]*0.5, point[1]-yerr), (iz[0] + iz[1]*0.5, point[1]+yerr), (iz[0] - iz[1]*0.5, point[1]+yerr), (iz[0] - iz[1]*0.5, point[1]-yerr)])
                poli_top = Polygon([(iz[0] + iz[1]*0.5, point[0]-xerr), (iz[0] + iz[1]*0.5, point[0]+xerr), (iz[0] - iz[1]*0.5, point[0]+xerr), (iz[0] - iz[1]*0.5, point[0]-xerr)])
                               
                ax[0].add_patch(PolygonPatch(poli_side, fc='cornflowerblue', ec='#555555', lw=0.2, alpha=0.7))                         
                ax[1].add_patch(PolygonPatch(poli_top, fc='cornflowerblue', ec='#555555', lw=0.2, alpha=0.7))       
                             
            if fi in nomatch_i:          
                xerr = width_list_x[oi]*0.5
                yerr = width_list_y[oi]*0.5
                #print "inner xerr, yerr: ", oi, xerr, yerr
                poli_side = Polygon([(iz[0] + iz[1]*0.5, point[1]-yerr), (iz[0] + iz[1]*0.5, point[1]+yerr), (iz[0] - iz[1]*0.5, point[1]+yerr), (iz[0] - iz[1]*0.5, point[1]-yerr)])
                poli_top = Polygon([(iz[0] + iz[1]*0.5, point[0]-xerr), (iz[0] + iz[1]*0.5, point[0]+xerr), (iz[0] - iz[1]*0.5, point[0]+xerr), (iz[0] - iz[1]*0.5, point[0]-xerr)])

       
                ax[0].add_patch(PolygonPatch(poli_side, fc='black', ec='#555555', lw=0.2, alpha=0.5))                         
                ax[1].add_patch(PolygonPatch(poli_top, fc='black', ec='#555555', lw=0.2, alpha=0.5))   #"""
                                                 
    for oz, fo in zip(of_z, curr_fO_c):
        for point, oo in zip(fO_proj, currO):  

            if [oo,fo] in match_o:
                xerr = width_list_xo[oo]*0.5
                yerr = width_list_yo[oo]*0.5
                polo_side = Polygon([(oz[0] + oz[1]*0.5, point[1]-yerr), (oz[0] + oz[1]*0.5, point[1]+yerr), (oz[0] - oz[1]*0.5, point[1]+yerr), (oz[0] - oz[1]*0.5, point[1]-yerr)])
                polo_top = Polygon([(oz[0] + oz[1]*0.5, point[0]-xerr), (oz[0] + oz[1]*0.5, point[0]+xerr), (oz[0] - oz[1]*0.5, point[0]+xerr), (oz[0] - oz[1]*0.5, point[0]-xerr)])
                              

                ax[0].add_patch(PolygonPatch(polo_side, fc='crimson', ec='#555555', lw=0.2, alpha=0.7))         
                ax[1].add_patch(PolygonPatch(polo_top, fc='crimson', ec='#555555', lw=0.2, alpha=0.7))   
                                     
            if fo in nomatch_o:
                xerr = width_list_xo[oo]*0.5
                yerr = width_list_yo[oo]*0.5
                polo_side = Polygon([(oz[0] + oz[1]*0.5, point[1]-yerr), (oz[0] + oz[1]*0.5, point[1]+yerr), (oz[0] - oz[1]*0.5, point[1]+yerr), (oz[0] - oz[1]*0.5, point[1]-yerr)])
                polo_top = Polygon([(oz[0] + oz[1]*0.5, point[0]-xerr), (oz[0] + oz[1]*0.5, point[0]+xerr), (oz[0] - oz[1]*0.5, point[0]+xerr), (oz[0] - oz[1]*0.5, point[0]-xerr)])
                               

                ax[0].add_patch(PolygonPatch(polo_side, fc='black', ec='#555555', lw=0.2, alpha=0.5))         
                ax[1].add_patch(PolygonPatch(polo_top, fc='black', ec='#555555', lw=0.2, alpha=0.5))   #"""           
                 

    fig.tight_layout()


    #plot_event_hodoscope(ax[2], h)
    plot_all_bars_2D(ax[2], circle, union_cl_I, union_cl_O, vertex_BGO)  
    
    
    fig.savefig("./event_plots_sims/event_3tr.pdf", dpi = 200) 
    #figt.savefig("./event_plots/top_{0}_{1}.pdf".format(cuspnr, eventnr), bbox_inches='tight', dpi = 200)      
    #plt.show()



########################################################################################################################################################
def flatten(l): return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]
########################################################################################################################################################
def track_finding_4layers(polygon_hitlist, h, nO, nI, I, O, IFC, OFC,
                           zerrI, zerrO, cuspnr, eventnr, circle, union_cl_I, union_cl_O,
                           sims_fibre_inner_hits, sims_fibre_outer_hits, sims_BGO_hit,
                           InnerPosX,InnerPosY,InnerPosZ,
                           OuterPosX,OuterPosY,OuterPosZ):

    print "----------------------------- NEW TRACKING FUNCTION, TESTING -----------------------------"
      
    InnerPosX = InnerPosX[I]
    OuterPosX = OuterPosX[O]
    InnerPosY = InnerPosY[I]
    OuterPosY = OuterPosY[O]
    InnerPosZ = InnerPosZ[I]
    OuterPosZ = OuterPosZ[O]
    
    zposII = InnerPosZ
    zposOO = OuterPosZ
   
    plot = False
    ax = 999
    if plot == True:
        fig, ax = plt.subplots(3,1, figsize = (7,13))
    
    
    pmapcI, pmapcO = getPositionMap()
    pmapcI, pmapcO = pmapcI[I], pmapcO[O] 
    
    numbers = np.arange(32)   
    currI = numbers[I]
    currO = numbers[O] 
    
    if len(currI) == 0 and len(currO) == 0:
        plt.close('all')
        #print "first if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        return ax, 999, 999, [999,999,999], [999,999,999], 999, [999,999,999,999,999], [999], [999,999,999], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]], 999, 999,999, [999,999,999], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]]
  
        
    
    list_i = [i for i,x in enumerate(I) if x == True]   
    list_o = [o for o,x in enumerate(O) if x == True]
    
    #print polygon_hitlist[0].bounds
    #print polygon_hitlist[0].bounds[1::2] # access y coordniates of bounderies of polygons
  
    ####
    # now I will create 'polygons' of all the hits in the hodoscope. I named it polygons because I am plotting the 2D projections
    # I will at the same time also write out the 3D case and name it as well polygons   
        
    ################################################################################
    #######################################################
    
    
    #if BGO is in use, the following lines get the boundaries of the hit cluster in the pixel map:
    
    if len(polygon_hitlist) > 0:
        # get the Y coordinates of the cluster in the BGO pixel map
        max_y = -100.0
        min_y = 100.0
        for p in polygon_hitlist:
            current = max(p.bounds[1::2])
            current2 = min(p.bounds[1::2])
            if current > max_y:
                max_y = current
            if current2 < min_y:
                min_y = current2
                
        print min_y, max_y 
        ###################################################
        # get the X coordinates of the cluster in the BGO pixel map
        max_x = -100.0
        min_x = 100.0
        for p in polygon_hitlist:
            current = max(p.bounds[::2])
            current2 = min(p.bounds[::2])
            if current > max_x:
                max_x = current
            if current2 < min_x:
                min_x = current2

     
    
    # Now we need the polygons of the hodoscope hits in the zy plane    
    numbers_fi = np.arange(63) # number array of the inner fibre channels
    numbers_fo = np.arange(100) # number array of the outer fibre channels 
    
    
    fO = OFC > 0 # active fibres outer layer (boolian array)
    fI = IFC > 0 # active finbres inner layer (boolian array)
    
    fI_not = np.logical_not(fI==True) # (inner fibres without a hit)
    fO_not = np.logical_not(fO==True)
    
    curr_fI = numbers_fi[fI] # channels numbers of hits inner fibres
    curr_fO = numbers_fo[fO] # channel numbers of hits outer fibres

    #print curr_fI
    curr_fI_c = []
    curr_fO_c = []
    for i in curr_fI:
        curr_fI_c.append([i])
        
    for i in curr_fO:
        curr_fO_c.append([i])    
       

    
    #ynfI, curr_fI_c = find_consecutive_nr_singletoo(curr_fI) # FIXME find clusters (neighbouring bars with hits) in the inner layer of fibres
    #ynfO, curr_fO_c = find_consecutive_nr_singletoo(curr_fO)

    curr_fI_not = numbers_fi[fI_not]
    curr_fO_not = numbers_fo[fO_not]
        
    #print "CLUSTER fI", curr_fI_c
    #print "CLUSTER fO", curr_fO_c   
    
    fi_dia = h.fi.diameter # get the diameters
    fo_dia = h.fo.diameter
    
    z_I_Fibre_hits = h.fi.getPositionMap()# get the position map of the inner fibre layer
    z_O_Fibre_hits = h.fo.getPositionMap()#[fO]
    
    of_z = [] # merge z positions together if there is a cluster
    if_z = []
    
    fI_proj, fO_proj = Get_Fibre_XYBar_Proj()

    for j in curr_fI_c: # loop over array where neighbouring channels are clustered. merge them (their z positions) together
        z = z_I_Fibre_hits[j]                
        err = len(z)*4.0 # len of fibres in cluster times 4 mm (width)     
        z = np.mean(z)     
        #print z
        if_z.append([z,err])#"""
    #if_z = np.hstack(z_I_Fibre_hits[curr_fI], np.array)
    for j in curr_fO_c: # loop over array where neighbouring channels are clustered. merge them (their z positions) together
        z = z_O_Fibre_hits[j]              
        err = len(z)*4.0 # mm
        z = np.mean(z) 
        #print z
        of_z.append([z,err])#"""
     
    #of_z = z_O_Fibre_hits[curr_fO]
   
         
    pmapI, pmapO = prepare_PosMap_nottracking(I,O,nI,nO) # get the x,y positions of the hodoscope bars (centers of bars)  
    
    # find matching bars nad fibres 
    """track_collection_new = bar_fibre_matching_2point0(pmapI, pmapO, zposII, zposOO, zerrI, zerrO, currI, currO, z_I_Fibre_hits, z_O_Fibre_hits, curr_fI, curr_fO)
    
    
    
    clean_track_coll = []    
    for i in track_collection_new:
        #print "cleaning..", len(i[0]), len(i[1])
        #
        #print "i", i
        #print "i0", i[0] 
        #print "i2", i[2]
        #print i[0][6], i[2][6]
        
        if i[2][6] == 0:
            if i[0][6] in [31,1]:# [30,31,1,2]: # look only at combinations where the inner hit is in range [outer-1, outer+1]
                clean_track_coll.append(i) 

        if i[2][6] == 31:
            if i[0][6] in [30,0] :#[29,30,0,1]: # look only at combinations where the inner hit is in range [outer-1, outer+1]
                clean_track_coll.append(i)             
                
        if i[0][6] in range(int(i[2][6] - 1), int(i[2][6] + 1+1)): # look only at combinations where the inner hit is in range [outer-1, outer+1]
            clean_track_coll.append(i) """
    
    #exit(0)
    # Old matching function
    match_i, match_o, nomatch_i, nomatch_o, barsi_left, barso_left = bar_fibre_matching(pmapI, pmapO, zposII, zposOO, zerrI, zerrO, currI, currO,
                                                                z_I_Fibre_hits, z_O_Fibre_hits, curr_fI_c, curr_fO_c)
                   
                   
    #print "machting inner hits:"
    #print match_i
    #print " "
    #print "machting outer hits:"
    #print match_o
    #print " "               
                                                                

    # x and y error of the bar hodoscope ####################
    x = 14.14   # cos(45grad) = x / 20.0      
    width_list_x = [5.0,5.0,5.0,5.0,x,x,x,x,20.0,20.0,20.0,20.0,x,x,x,x,5.0,5.0,5.0,5.0, x,x,x,x,20.0,20.0,20.0,20.0,x,x,x,x]
    width_list_y = [20.0,20.0,20.0,20.0, x,x,x,x,5.0,5.0,5.0,5.0,x,x,x,x,20.0,20.0,20.0,20.0,x,x,x,x,5.0,5.0,5.0,5.0,x,x,x,x]
    
    x = 24.75  
    width_list_xo = [5.0,5.0,5.0,5.0,x,x,x,x,35.0,35.0,35.0,35.0,x,x,x,x,5.0,5.0,5.0,5.0, x,x,x,x,35.0,35.0,35.0,35.0,x,x,x,x]
    width_list_yo = [35.0,35.0,35.0,35.0, x,x,x,x,5.0,5.0,5.0,5.0,x,x,x,x,35.0,35.0,35.0,35.0,x,x,x,x,5.0,5.0,5.0,5.0,x,x,x,x]
    #########################################################

    bi_p_3D = [] # 3D points of the inner bar layer
    bo_p_3D = [] # 3D points of the outer bar layer
    fi_p_3D = [] # 3D points of the inner fibre layer
    fo_p_3D = [] # 3D points of the outer fibre layer        
    # structure:
    #[z,zerr, x, xerr, y, yerr, barnr, fibrenr[0], length of fibre cluster]  
       
    #####
    
    # if you want to add also all bar combinations inner/outer 
    barsi_left = currI
    barso_left = currO
 
    ##### 
       
    if len(barsi_left) != 0:   # if there are bars with a hit but no matching fibre hit
        # get the indices of the bars left in the currO array - so we can get the corresponding x,y,z vals
        inter =  np.in1d(currI, barsi_left)        
        indices=np.array(range(len(currI)))[inter]
        
        xy = pmapcI[indices]
        z = zposII[indices]
        #widx = np.array(width_list_x)[indices]
        #widy = np.array(width_list_y)[indices]
        for barsi, xy_, z_ in zip(barsi_left, xy, z):
            # x and y error devided by root(12). assuming uniform distribution in x and y 
            errxx, erryy = Calc_xy_error(barsi, 0)  
            #print "calc func",barsi, errxx, erryy
            #print "by width", widx_/3.46, widy_/3.46    
            bi_p_3D.append([z_, zerrI/2.355, xy_[0], errxx, xy_[1], erryy, barsi, 999, 0, 0])

    
    if len(barso_left) != 0:
        # get the indices of the bars left in the currO array - so we can get the corresponding x,y,z vals
        inter =  np.in1d(currO, barso_left)        
        indices = np.array(range(len(currO)))[inter]

        xy = pmapcO[indices]
        z = zposOO[indices]
        #widx = np.array(width_list_xo)[indices]
        #widy = np.array(width_list_yo)[indices]
        for barso, xy_, z_ in zip(barso_left, xy, z):
            errxx, erryy = Calc_xy_error(barso, 1)
            #print "calc func",barso, errxx, erryy
            #print "by width", widx_/3.46, widy_/3.46
            # x and y error devided by root(12). assuming uniform distribution in x and y       
            bo_p_3D.append([z_, zerrO/2.355, xy_[0], errxx, xy_[1], erryy, barso, 999, 0, 1])   
            
    if_z = np.array(if_z) # just convert to a numpy array    
    of_z = np.array(of_z) 

    ###############################################################################################################################################        
    fI_proj = fI_proj[I]
    fO_proj = fO_proj[O]
    
    for iz, fi in zip(if_z, curr_fI_c):
        for point, oi in zip(fI_proj, currI):  

            if [oi,fi] in match_i:      
                #xerr = width_list_x[oi]/3.46
                #yerr = width_list_y[oi]/3.46
                errxx, erryy = Calc_xy_error(oi, 0)
                #print "calc func",oi, errxx, erryy
                # this is an array cotaining all the 3D hits in the inner layer:
                #print "sqrt erro 12222222222222222", iz[1]/math.sqrt(12.0)
                fi_p_3D.append([iz[0], iz[1]/math.sqrt(12.0), point[0], errxx, point[1], erryy, oi, fi[0], len(fi), 0])
      
                                          
    for oz, fo in zip(of_z, curr_fO_c):
        for point, oo in zip(fO_proj, currO):  

            if [oo,fo] in match_o:
                errxx, erryy = Calc_xy_error(oo, 1)
                #print "calc func",oi, errxx, erryy
                # x and y error devided by root(12) = 3.46.... assuming uniform distribution in x and y                     
                fo_p_3D.append([oz[0], oz[1]/math.sqrt(12.0), point[0], errxx, point[1], erryy, oo, fo[0], len(fo), 1])

    ################################################################################################################################################
    ################################################################################################################################################
    #
    ################################################################################################################################################
    ################################################################################################################################################    
    
    # combine inner fibre and inner bar hits to point pairs #################################
    col_fbi_hits = [] 
    col_fbo_hits = [] 
    
    # entry with index 6 is the bar number - if they are equal, add     
    for fip in fi_p_3D:
        for bip in bi_p_3D:            
            if fip[6] == bip[6]:
                col_fbi_hits.append([bip, fip]) # combination of inner fibre and inner bar points that will be part of the same track candidate
            
    for fop in fo_p_3D:
        for bop in bo_p_3D:            
            if fop[6] == bop[6]:         
                col_fbo_hits.append([bop, fop])
                
                
    #for i in col_fbi_hits:
    #    print "haha", i             
                
    # I also add hits in the bar hodoscope that dont have a fibre hit: ####################        
    for bip in bi_p_3D:
        check_bool = False
        for fbi in col_fbi_hits:    
            if len(fbi) == 2:
                # if none of the fibre-bar combinations (inner) has the bar included
                if bip[6] == fbi[0][6] or bip[6] == fbi[1][6]:
                    check_bool = True           
                   
        if check_bool == False:
            col_fbi_hits.append([bip])
            #print "BINGOOOOOOOO bip", bip

    for bop in bo_p_3D:
        check_bool = False
        for fbo in col_fbo_hits:          
            if len(fbo) == 2:
                # if none of the fibre-bar combinations (outer) has the bar included
                if bop[6] == fbo[0][6] or bop[6] == fbo[1][6]:
                    check_bool = True 
                 
        if check_bool == False:
            col_fbo_hits.append([bop])
            #print "BINGOOOOOOOO bop", bop
    ############################################################################################                    
    
    
                     
    # all point combinations of inner and outer ################################################
    all_combs = [[x,y] for x in col_fbi_hits for y in col_fbo_hits]
    
    
    #print "leeeeeeeeeeeeeen", len(all_combs)  
    
    #for c in all_combs:
    #    for j in c:
    #        print j
    #    print " "
    
    # remove those that dont have bars close to each other #####################################
    clean_combs = []    
    for i in all_combs:
        #print "cleaning..", len(i[0]), len(i[1])
        #print i[0][0][6], i[1][0][6]
        #print "range", range(int(i[1][0][6] - 1), int(i[1][0][6] + 1+1))        
        if i[1][0][6] == 0:
            if i[0][0][6] in [31,1]:#[30,31,1,2]: # look only at combinations where the inner hit is in range [outer-1, outer+1]
                clean_combs.append(i) 

        if i[1][0][6] == 31:
            if i[0][0][6] in [30,0]:#[29,30,0,1]: # look only at combinations where the inner hit is in range [outer-1, outer+1]
                clean_combs.append(i)             
                
        if i[0][0][6] in range(int(i[1][0][6] - 1), int(i[1][0][6] + 1+1)): # look only at combinations where the inner hit is in range [outer-1, outer+1]
            clean_combs.append(i) 
    
    # the next bit just brings the track candidates in a nicer form, the above creates nested lists... in principle I can make the above nicer
    # programming style wise, but  this will only happen if I have time in the end... for now this works and thats the main point

    new_clean_combs = []    
    for hhh in clean_combs:
        
        curr_comb = []

        for hui in hhh:
            flat = flatten(hui)

            if len(flat) > 10:
                curr_comb.append(flat[0:10])
                curr_comb.append(flat[10:20])
            else:
                curr_comb.append(flat)
                
        new_clean_combs.append(curr_comb) 


    new_clean_combs = [s for s in new_clean_combs if len(s) > 2]
        
    #######################################################    
            
    # start the fitting ###############################################################################        
    
            
            
    print "leeeeeeeeeeeeeen", len(new_clean_combs)         
    print "-------------------------- FITTING THE TRACK CANDIDATES ----------------------------------"  
    # now the tracks will be fitted. different functions are defined in this file using different parameterisations of the line. 
    # The "angles" version uses a point given by (px, py, pz), where we can set one of them zero (py = 0 in this case) and a 
    # direction vector defined by 2 angles phi and theta
       
    #new_clean_combs = clean_track_coll

    
    if parameterisation == 'angle': 
        lines_fit, chisq, fit_status, cov_mats, cov_stats = track_selection_fitting_3D_4layers_root_angles(new_clean_combs) #FIXME
    elif parameterisation == 'p_d':
        lines_fit, chisq, fit_status, cov_mats, cov_stats = track_selection_fitting_3D_4layers_root_new(new_clean_combs)


    clean_combs = np.array(new_clean_combs)
   
    lines_fit = np.array(lines_fit)
    
    ####################################################  
    if len(chisq) == 0: # if no lines at all have been found
        #plt.show()
        
        print "second if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        if plot == True:
            plot_event_3_views(fig, ax, h, list_i, list_o, pmapcI, pmapcO, zposII, zposOO, zerrI, zerrO, width_list_x, width_list_xo, width_list_y, width_list_yo,
                                    currI, currO, match_i, match_o, nomatch_i, nomatch_o, if_z, of_z, curr_fI_c, curr_fO_c, circle, union_cl_I, union_cl_O, sims_BGO_hit,
                                    fI_proj, fO_proj)
            plt.show()                        
        plt.close('all')                         
        return ax, 999, 999, [999,999,999], [999,999,999], 999, [999,999,999,999,999], [999], [999,999,999], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]], 999, 999, 999, [999,999,999], [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]]
  
    
    """print "check normalisation direction vector..."
       
    for line in lines_fit:
        print "line: ", line
        px = line[0]
        pz = line[1]
        dx = line[2]
        dz= line[3]
        norm_line = np.linalg.norm(np.array([dx, 1.0, dz ]))
        
        
        print norm_line #"""
    
    
    
    ####################################################################################################
    # sort according to least sum of squares of fit
    
    
    #sorted_clean_combs, sorted_chisq, sorted_lines_fit, sorted_cov_mats, sorted_cov_stats, sorted_fit_status = smart_sort_fitted_tracks(clean_combs, chisq, lines_fit, cov_mats, cov_stats, fit_status)

    sorted_indices = [b[0] for b in sorted(enumerate(chisq),key=lambda i:i[1])]
    
    sorted_chisq = chisq[sorted_indices]
            
    sorted_clean_combs = clean_combs[sorted_indices]
    
    sorted_lines_fit = lines_fit[sorted_indices] 
    
    sorted_cov_mats = cov_mats[sorted_indices]
    
    sorted_cov_stats = cov_stats[sorted_indices]
    
    sorted_fit_status = fit_status[sorted_indices] #"""
   
    """mask0 = fit_status == 0
    mask1 = np.invert(mask0)
    
    sorted_chisq0 = chisq[mask0]
    sorted_chisq1 = chisq[mask1]
    
    print chisq
    print sorted_chisq0
    
    sorted_indices0 = [b[0] for b in sorted(enumerate(sorted_chisq0),key=lambda i:i[1])]
    sorted_indices1 = [b[0] for b in sorted(enumerate(sorted_chisq1),key=lambda i:i[1])]
    
    print "sorted inds 0", sorted_indices0
    print "sorted inds 1", sorted_indices1
    
    sorted_chisq0 = sorted_chisq0[sorted_indices0]
    sorted_chisq1 = sorted_chisq1[sorted_indices1]
    sorted_chisq = np.append(sorted_chisq0, sorted_chisq1)
            
    sorted_clean_combs0 = clean_combs[mask0][sorted_indices0]
    sorted_clean_combs1 = clean_combs[mask1][sorted_indices1]
    sorted_clean_combs = np.vstack((sorted_clean_combs0,sorted_clean_combs1))
    
    sorted_lines_fit0 = lines_fit[mask0][sorted_indices0]
    sorted_lines_fit1 = lines_fit[mask1][sorted_indices1]
    sorted_lines_fit = np.vstack((sorted_lines_fit0,sorted_lines_fit1))
    
    sorted_fit_status0 = fit_status[mask0][sorted_indices0]
    sorted_fit_status1 = fit_status[mask1][sorted_indices1]
    sorted_fit_status = np.append(sorted_fit_status0,sorted_fit_status1)
    
    sorted_cov_mats0 = cov_mats[mask0][sorted_indices0]
    sorted_cov_mats1 = cov_mats[mask1][sorted_indices1]
    sorted_cov_mats = np.vstack((sorted_cov_mats0,sorted_cov_mats1))
    
    sorted_cov_stats0 = cov_stats[mask0][sorted_indices0]
    sorted_cov_stats1 = cov_stats[mask1][sorted_indices1]
    sorted_cov_stats = np.append(sorted_cov_stats0,sorted_cov_stats1) # """
      
    chisq = sorted_chisq
    lines_fit = sorted_lines_fit
    clean_combs = sorted_clean_combs
    fit_status = sorted_fit_status
    cov_stats = sorted_cov_stats
    cov_mats = sorted_cov_mats
    
    print "chisq, ", chisq        
          
    
    if len(chisq) != 0:

        # colour for plotting lines
        colours = ['gold', 'sandybrown', 'orangered', 'hotpink', 'mediumorchid','purple', 'thistle', 'skyblue', 'blue', 'darkturquoise', 'limegreen', 'darkgreen', 'greenyellow','grey','black']
        
        print " "
        print "number of lines fitted: ", len(chisq)
        print " "
                   
        print "  "         
        #############################################################################################################################################################
        #############################################################################################################################################################
        # now we need to select tracks out of the track collection 
        angle_cut = False          # cut that removes tracks that are parallel
        fibre_comb_cut = True      # cut that removes all tracks that contain a fibre-bar combination that has already been used
        # INFO: definetly use the fibre_comb_cut - the angle cut was just for a test, I only left it in this python file for testing
        #############################################################################################################################################################
        # cut away lines that dont go through the detector (close to the center (depends on the bounds))
        # need to go at least though the inner hodoscope:
        #xbound_vertex = 150.0 # mm
        #ybound_vertex = 150.0 # mm
        
        clean_combs_1 = []  # save selected tracks
        chisq_1 = []   # save chi2 of selected tracks
        lines_fit_after_select_1 = []  # save fit parameter of selected tracks
        cov_mats_1 = [] # save covariance matrices of selected tracks
        status_fit = []       
               
        combs_copy = clean_combs #[:,:,6:10].tolist()  # copy bar/fibre hit combinations (track point collection)
                
        fibres_already_added = []
        bars_already_added = []
        
        # do this check only if the line is not close to parallel to the y-axis (calculation fails in that case)
        counter = -1 # this is just for plotting the fitted lines in colour
        counter2 = -1 # this is just for plotting the fitted and selected lines in colour
        
        for n, (paramsl, cq, track, fstat, covm, covs) in enumerate(zip(lines_fit, chisq, clean_combs, fit_status, cov_mats, cov_stats)):
            ###################################
            # get x,y of the line at z = 0
            #x, y = get_xy_of_line_new(0.0,paramsl)
            #pz0 = [x, y, 0.0]
            #print x, y
            ##################################
            # calculate angle with y-axis:
            if parameterisation == 'angle':
                x0,y0,z0 = line_param_angles(1000.0,paramsl)    # get a point on the current line (this is for plotting)
                x1,y1,z1 = line_param_angles(-1000.0,paramsl)   # get another point
            elif parameterisation == 'p_d':
                x0,y0,z0 = line_param_new(1000.0,paramsl)    # get a point on the current line (this is for plotting)
                x1,y1,z1 = line_param_new(-1000.0,paramsl)   # get another point
                            
            p0 = np.array([x0,y0,z0])
            p1 = np.array([x1,y1,z1])
           
            #d1 = p1 - p0
            #print "d1", d1
            #angle_with_y_axis = angle_3Dtracks(d1, np.array([0.0,1.0,0.0]))
            #angle_with_my_axis = angle_3Dtracks(d1, np.array([0.0,-1.0,0.0]))
            #angle_with_y_axis =  angle_with_y_axis/(2.0*np.pi)*360.0
            #angle_with_my_axis =  angle_with_my_axis/(2.0*np.pi)*360.0
            
            #cutting_angle = min(angle_with_y_axis, angle_with_my_axis)
            #print "angle_with_y_axis",angle_with_y_axis
            #print "angle_with_my_axis",angle_with_my_axis
            #print  "abs(y)", abs(y),  "abs(x)", abs(x) #, color     """
                  
            # plot all found lines, with low alpha:      
            if plot == True:
                #ind = n
                counter+=1
                ind = counter
                if counter >= len(colours): # if there are more lines than the colour array is long, just start from the beginning again
                    ind = int(math.fmod(n, len(colours)))
                        
                color = colours[ind] # select a colour to plot

                #ax[0].plot([z0, z1],[y0,y1], '--', c = color, alpha = 1)
                #ax[1].plot([z0, z1],[x0,x1], '--', c = color, alpha = 1)
                #ax[2].plot([x0, x1],[y0,y1], '--', c = color, alpha = 1)

                #print " "
                #print "color ", color, "chi2", cq
                #for fff in track:
                #    print np.array(fff)[[6,7,8,9]]
                
                
            if fstat > 0:
                #print "SKIP!!!!!!!!!!!!!!!! fit status"
                continue

            if covs < 2 and covs < 0:
            #    #print "SKIP!!!!!!!!!!!!!!!! cov matrix"
                continue                
            #if cutting_angle > 15:
            #    if (abs(y) > ybound_vertex or abs(x) > xbound_vertex):
            #            print "continue"
            #            continue
                               
            ###################################################################################   
            # in case of angle cut... 
            if  angle_cut == True and fibre_comb_cut == False: 
                lines_fit_after_select_1.append(paramsl)
                clean_combs_1.append(track)
                chisq_1.append(cq)
                cov_mats_1.append(covm)
                
                
            ###################################################################################
            # here the selection accoring to fibre/bar combinations happens
            #print "--------"  
            #print "track", track
            #print " "
            # here the selection accoring to fibre/bar combinations happens
            if fibre_comb_cut == True:
                add = True
                for blub in track:                 
                    fibre_comb = list( blub[i] for i in [6,7,9] )
                    #print fibre_comb
                    if fibre_comb in fibres_already_added:# and fibre_comb[0] != 999:# and 
                        add = False
                        break
                    
                if add == True:
                    for blub in track:
                        fibre_comb = list( blub[i] for i in [6,7,9] )
                        #if fibre_comb[0] != 999:                   
                        fibres_already_added.append(fibre_comb)#"""
                            
            ##############################################################################################
            
                
                if add == True:
                    #print "add track: ", np.array(track)[:,[6,7,9]]
                    lines_fit_after_select_1.append(paramsl)
                    clean_combs_1.append(track)
                    chisq_1.append(cq)
                    cov_mats_1.append(covm)
                    status_fit.append([fstat, covs])
                    #color_after_select.append(color)
                    if plot == True:
                        counter2+=1
                        ind = counter2
                        if n >= len(colours): # if there are more lines than the colour array is long, just start from the beginning again
                            ind = int(math.fmod(n, len(colours)))
                                
                        color = colours[ind] # select a colour to plot
                        
                        #ax[0].plot([z0, z1],[y0,y1], '-', c = color)
                        #ax[1].plot([z0, z1],[x0,x1], '-', c = color)
                        #ax[2].plot([x0, x1],[y0,y1], '-', c = color)
                        #print "------------------------------------------------------------------------------------------"
                        #print "added! ", color, np.array(track)[:,[6,7,9]]
                        #print " " #"""

            ###################################################################################
            ###################################################################################
            
        ############################################################################################################################################################       
        #############################################################################################################################################################

          
        if angle_cut == True:
            # get rid of parallel lines that go through similar bar combinations  
            lines_fit_after_select = select_tracks_angle(lines_fit_after_select_1, chisq_1, clean_combs_1)
        else:
            lines_fit_after_select = lines_fit_after_select_1

        

        
        
        print "lines after select ............................................................."
        for i in lines_fit_after_select:
            print i 
        
        print "......................................................................................"
        print "total lines: ", len(lines_fit_after_select)         
        print "......................................................................................"
        
        if plot == True:
            for n, (params, cq, cov, stat, track_) in enumerate(zip(lines_fit_after_select, chisq_1, cov_mats_1, status_fit, clean_combs_1)):           
                ind = n
                if n >= len(colours): # if there are more lines than the colour array is long, just start from the beginning again
                    ind = int(math.fmod(n, len(colours)))
                    
                color = colours[ind] # select a colour to plot
                if angle_cut == True:
                    pars = params[0]
                else:
                    pars = params
                    
                    
                if parameterisation == 'angle':
                    x0,y0,z0 = line_param_angles(4000.0,pars)    # get a point on the current line (this is for plotting)
                    x1,y1,z1 = line_param_angles(-4000.0,pars)   # get another point on the line
                elif parameterisation == 'p_d':
                    x0,y0,z0 = line_param_new(4000.0,pars)    # get a point on the current line (this is for plotting)
                    x1,y1,z1 = line_param_new(-4000.0,pars)   # get another point

                ax[0].plot([z0, z1],[y0,y1], '-', c = color, lw = 2)
                ax[1].plot([z0, z1],[x0,x1], '-', c = color, lw = 2)
                ax[2].plot([x0, x1],[y0,y1], '-', c = color, lw = 2) #"""
                
                #print "line: ", color 
                #print "track ", np.array(track_)[:,[6,7,8,9]], ".... fit status:", stat[0]
                #print cov
                #print "..............................."
               
        #############################################################################################################################################################        
        #############################################################################################################################################################            
        #############################################################################################################################################################            
                    
        dist_points = []
        m_cluster_arr = []
        mean_d_pt = [999,999,999] # create a dummy, in case it cant find a vertex
        min_point = [999,999,999]
        av_d_pt = [999,999,999]
        vertex_fitted = [999,999,999]
        dist_points = 999
        vertex_fitted_covar = [[999,999,999],[999,999,999],[999,999,999]]
        covmat_largest_cl = [[999,999,999],[999,999,999],[999,999,999]]
        total_covar_av = [[999,999,999],[999,999,999],[999,999,999]]
        mean_clus_xyz_ =[999,999,999]
        total_covar_mean_ =[[999,999,999],[999,999,999],[999,999,999]]
        covar_largest_cl_mn_ = [[999,999,999],[999,999,999],[999,999,999]]
        dists = [999]
        cl_pts = []
        n_clust = 0
        av_dist_l_cl = [999]
        nr_lines_fit_after_select_unique = 999
        
        hesse_bool = 999       
        if len(clean_combs_1) > 1:
            #print "IN! 1"
        
            ######################################
            # this function creates a grid / 3D matrix around the center of the detector
            # then it calculates for all points on that grid the distance to the lines
            # it returns the point with the minimum sum of distances
            # it works but is very, very computing power intense! 
                      
            #min_point = dist_lines_volume(lines_fit_after_select)
            ######################################
        
            if len(lines_fit_after_select) > 1:
                # cov_mats_1 ... covariance matrices for all fits of the selected tracks
                # angle_cut... boolian
                # lines_fit_after_select... collection of all the parameter values of all track fits
                    
                # return values:
                # dist_points ... array of mid points between lines, [x,y,z, ind1, ind2]
                # n_cluster ... number of lcuster found
                # m_cluster_arr ... [x,y,z, color, number of points, density] (x,y,z) of mean point
                # dists ...array containing the distances of line pairs
                # cl_pts ... array with the closests points on lines of line pairs
                # largest_cl ... largest cluster, entry of m_cluster_arr
                # covmat_largest_cl ... covariance matrix of propagated error: fit errors -> mean of midpoints
                # total_covar_av ... same as above but not only for the largest cluster but for all tracks
                # inds_largest_cluster ... pairs of indices of tracks of the largest cluster
               
                dist_points, n_clust, m_cluster_arr, dists, cl_pts, largest_cl, covmat_largest_cl, total_covar_av, inds_largest_cl, av_dist_l_cl, mean_clus_xyz_, total_covar_mean_,  covar_largest_cl_mn_ = cluster_finding_midpoints(lines_fit_after_select, angle_cut, cov_mats_1) #, chisq, clean_combs)
                unique_l_cl_inds = np.unique(inds_largest_cl).astype(int)
                
               
                print "unique inidces of tracks in cluster", unique_l_cl_inds
                #if n_clust== 0:
                    #plt.close('all') 
                    #return len(lines_fit_after_select), [999,999,999], [999,999,999], 999, [999,999,999,999,999]
                    
                print " "
                print "..........................." 
                print "av dist largest cl",  np.mean(av_dist_l_cl)  
                           
                lines_fit_after_select_unique =  [lines_fit_after_select[i] for i in unique_l_cl_inds]  
                nr_lines_fit_after_select_unique = len(lines_fit_after_select_unique)
                #nr_lines_fit_after_select_uniqu
                #cov_mats_1_unique =  [cov_mats_1[i] for i in unique_l_cl_inds]  
                
                mean_d_pt = np.average(dist_points[:,[0,1,2]], weights = 1.0/(dists), axis = 0)                 
                vertex_estimate = largest_cl[[2,3,4]]
                
                
                # give:
                # vertex_estimate ... first estimate of vertex (use mean of all midpoints of lines or mean of largest cluster)
                # fit parameters ... use those of all selected tracks, or those that are part of largest cluster of midpoints
                # covariance matrices of fit results               

                """if vertexfit_choice == 'manual':
                    if parameterisation == 'angle':
                        vertex_fitted, vertex_fitted_covar = do_the_vertexfit(vertex_estimate, lines_fit_after_select, cov_mats_1)
                        vertex_fitted, vertex_fitted_covar = do_the_vertexfit(vertex_estimate, lines_fit_after_select, cov_mats_1)                
                    elif parameterisation == 'p_d':
                        vertex_fitted, vertex_fitted_covar, hesse_bool = do_the_vertexfit_2(vertex_estimate, lines_fit_after_select, cov_mats_1)      
                elif vertexfit_choice == 'root':
                    vertex_fitted, vertex_fitted_covar = do_the_vertexfit_root(vertex_estimate, lines_fit_after_select, cov_mats_1, parameterisation)        
                elif vertexfit_choice == 'python':
                    vertex_fitted, vertex_fitted_covar = do_the_vertexfit_pythonopt(vertex_estimate, lines_fit_after_select, cov_mats_1, parameterisation)    # """     
                    
                    
                                   
                """for pars,tracks in zip(lines_fit_after_select, clean_combs_1):
                    print len(lines_fit_after_select), len(clean_combs_1)
                    #print tracks.shape
                    print vertex_fitted
                    print tracks[-1][2]
                    print " "
                                        
                    x1 = vertex_fitted[0]
                    x2 = tracks[-1][2]

                    y1, dum = get_zy_of_line_angles(x1, pars)
                    y2, dum = get_zy_of_line_angles(x2, pars)
                    p0 = [x1,y1]
                    p1 = [x2,y2]
                    slope = slope_from_points(p0, p1)
                    intercept = p0[1] - slope*p0[0] 
               
                    x = np.zeros(2)
                    x[0] = x1
                    x[1] = x2   
                    data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
                    line = mpl_lines.Line2D(x, data_y, color='red', lw = 1)
                    ax[2].add_line(line)  # """   
             
                print "vertex estimate: ", vertex_estimate
                print " "
                #print covmat_largest_cl
                if len(covmat_largest_cl) > 0:
                    if plot == True: 
                        print "vertex fitted:", vertex_fitted, "error:", np.sqrt(np.diag(vertex_fitted_covar))
                        print " "
                        print ".............. errors on x,y,z according to covariance matrix (using clustering)", np.sqrt(np.diag(covmat_largest_cl))
                        print "..... errors on x,y,z according to covariance matrix (using total covar. matrix)", np.sqrt(np.diag(total_covar_av))
                else:
                    covmat_largest_cl = [[999,999,999],[999,999,999],[999,999,999]]    
                # just the mean of all mid points
                #if n_clust >= 0:
                    
                
                    #av_d_pt = np.average(dist_points[:,[0,1,2]], weights=1.0/np.array(dists), axis =  0) #"""
                
                    if plot == True:
                        print 
                        for f in dist_points:
                            ax[2].plot(f[0], f[1], 'x',c = 'black',  markersize = 5)
                            ax[1].plot(f[2], f[0], 'x', c = 'black',markersize = 5) 
                            ax[0].plot(f[2], f[1], 'x', c = 'black',markersize = 5)  #"""      
                            #print f
                        """colours = ["green", "red", "blue", "yellow"]    
                        for n,f in enumerate(cl_pts):
                            c = colours[n]
                            ax[2].plot(f[0][0], f[0][1], 'o', c = c, markersize = 5)
                            ax[1].plot(f[0][2], f[0][0], 'o',c = c, markersize = 5) 
                            ax[0].plot(f[0][2], f[0][1], 'o', c = c,markersize = 5)
                            ax[2].plot(f[1][0], f[1][1], 'o', c = c,markersize = 5)
                            ax[1].plot(f[1][2], f[1][0], 'o', c = c,markersize = 5) 
                            ax[0].plot(f[1][2], f[1][1], 'o', c = c, markersize = 5) #"""                                           
                        # mean of all the mid points
                        #mp = mean_d_pt                                
                        #ax[2].plot(mp[0], mp[1], 'go', markersize = 20)
                        #ax[1].plot(mp[2], mp[0], 'go', markersize = 20) 
                        #ax[0].plot(mp[2], mp[1], 'go', markersize = 20) 

                        # weighted mean of all the mid points
                        #mp = av_d_pt
                        #ax[2].plot(mp[0], mp[1], 'bo', markersize = 8)
                        #ax[1].plot(mp[2], mp[0], 'bo', markersize = 8) 
                        #ax[0].plot(mp[2], mp[1], 'bo', markersize = 8)  #"""
            
        if plot == True:    
        
            

            ax[2].plot(sims_BGO_hit[0], sims_BGO_hit[1], 'd', markerfacecolor = 'tomato',  markersize=10, alpha = 1)
            ax[1].plot(sims_BGO_hit[2], sims_BGO_hit[0], 'd', markerfacecolor = 'tomato',  markersize=10, alpha = 1)    
            ax[0].plot(sims_BGO_hit[2], sims_BGO_hit[1], 'd', markerfacecolor = 'tomato',  markersize=10, alpha = 1)
            
            if vertex_fitted[0] != 999:
                mp = vertex_fitted

                               
                ax[2].plot(mp[0], mp[1], 'o', color = 'lightgreen', markersize = 7)
                ax[1].plot(mp[2], mp[0], 'o', color = 'lightgreen', markersize = 7) 
                ax[0].plot(mp[2], mp[1], 'o', color = 'lightgreen', markersize = 7)   


            """for ix, iy, iz in zip(InnerPosX,InnerPosY,InnerPosZ):
                        ax[2].plot(ix, iy, 'o', markerfacecolor = 'black',  markersize=3)
                        ax[1].plot(iz, ix, 'o', markerfacecolor = 'black',  markersize=3)    
                        ax[0].plot(iz, iy, 'o', markerfacecolor = 'black',  markersize=3)             

            for ox, oy, oz in zip(OuterPosX,OuterPosY,OuterPosZ):
                        ax[2].plot(ox, oy, 'o', markerfacecolor = 'black',  markersize=3)
                        ax[1].plot(oz, ox, 'o', markerfacecolor = 'black',  markersize=3)    
                        ax[0].plot(oz, oy, 'o', markerfacecolor = 'black',  markersize=3)     

            sims_fibre_inner_hits = sims_fibre_inner_hits.T[curr_fI]
            sims_fibre_outer_hits = sims_fibre_outer_hits.T[curr_fO]
            
            #print " SIMS INNER ", sims_fibre_inner_hits
            #print " SIMS OUTER ", sims_fibre_outer_hits
            
            for sims_if in sims_fibre_inner_hits:
                ax[2].plot(sims_if[0], sims_if[1], 'o', markerfacecolor = 'black',  markersize=3)
                ax[1].plot(sims_if[2], sims_if[0], 'o', markerfacecolor = 'black',  markersize=3)    
                ax[0].plot(sims_if[2], sims_if[1], 'o', markerfacecolor = 'black',  markersize=3) 
                
            for sims_if in sims_fibre_outer_hits:
                ax[2].plot(sims_if[0], sims_if[1], 'o', markerfacecolor = 'black',  markersize=3)
                ax[1].plot(sims_if[2], sims_if[0], 'o', markerfacecolor = 'black',  markersize=3)    
                ax[0].plot(sims_if[2], sims_if[1], 'o', markerfacecolor = 'black',  markersize=3) #"""
            #sims_BGO_hit 
                               
            plot_event_3_views(fig, ax, h, list_i, list_o, pmapcI, pmapcO, zposII, zposOO, zerrI, zerrO, width_list_x, width_list_xo, width_list_y, width_list_yo,
                                currI, currO, match_i, match_o, nomatch_i, nomatch_o, if_z, of_z, curr_fI_c, curr_fO_c, circle, union_cl_I, union_cl_O, sims_BGO_hit,
                                fI_proj, fO_proj)
            if len(m_cluster_arr) > 0:             
                        """for mc in m_cluster_arr:                     
                            ax[2].plot(mc[0][0], mc[0][1], 'o', markerfacecolor = mc[1],  markersize=10)
                            ax[1].plot(mc[0][2], mc[0][0], 'o', markerfacecolor = mc[1],  markersize=10)    
                            ax[0].plot(mc[0][2], mc[0][1], 'o', markerfacecolor = mc[1],  markersize=10)     
                            
                            p1 = np.array(mc[0])
                            p2 = np.array(vertex_BGO)
                            #print p1, p2
                            squared_dist = np.sum(p1**2 + p2**2, axis=0)
                            dist = np.sqrt(squared_dist)        
                            print "dist to BGO vertex", dist, "point:", mc[0], mc[2], mc[3]#"""
                        mc = largest_cl
                        ax[2].plot(mc[2], mc[3], 'o', markerfacecolor = 'r',  markersize=8, alpha = 0.5)
                        ax[1].plot(mc[4], mc[2], 'o', markerfacecolor = 'r',  markersize=8, alpha = 0.5)   
                        ax[0].plot(mc[4], mc[3], 'o', markerfacecolor = 'r',  markersize=8, alpha = 0.5)          #"""
                         

    if len(m_cluster_arr) > 0:    
        
        biggest_cl = largest_cl
    else:
        biggest_cl = [999,999,999,999,999]    
        
    #if n_clust > 0:
    #plt.show()          
    if plot == True:   # FIXME
        #if (av_d_pt[2] < -40 or av_d_pt[2] > 40):
            plt.show()
            plt.close('all') 
        #else:
       #     plt.close('all') 
            
    return ax, dist_points, len(lines_fit_after_select), mean_d_pt, av_d_pt, n_clust, biggest_cl, dists, vertex_fitted, vertex_fitted_covar, covmat_largest_cl, total_covar_av, np.mean(av_dist_l_cl), hesse_bool, nr_lines_fit_after_select_unique, mean_clus_xyz_, total_covar_mean_,  covar_largest_cl_mn_ # biggest cluster 
    
######################################################################################################################################################################### 
#########################################################################################################################################################################
######################################################################################################################################################################### 
#########################################################################################################################################################################     
def cluster_finding_midpoints(lines_fit_, angle_cut_, cov_mats_): #, chisq_, clean_combs_):
        # cov_mats_ ... covariance matrices for all fits of the selected tracks
        # angle_cut_... boolian
        # lines_fit_... collection of all the parameter values of all track fits
        
        # return values:
        # dist_points ... array of mid points between lines, [x,y,z, ind1, ind2]
        # n_cluster ... number of lcuster found
        # m_cluster_arr ... [x,y,z, color, number of points, density] (x,y,z) of mean point
        # dists ...array containing the distances of line pairs
        # cl_pts ... array with the closests points on lines of line pairs
        # largest_cl ... largest cluster, entry of m_cluster_arr
        # covmat_largest_cl ... covariance matrix of propagated error: fit errors -> mean of midpoints
        # total_covar_av ... same as above but not only for the largest cluster but for all tracks
        # inds_largest_cluster ... pairs of indices of tracks of the largest cluster        
        print "--------------------------------- VERTEX FINDIND -----------------------------------------"
        print "vertex finding via cluster finding in mid point cloud of closest points of line pairs"
        
        mean_nps = []
        dists = []
        closest_points =  []
        covar_arr = []
        total_covar_mean = [[999,999,999],[999,999,999],[999,999,999]]
        mean_clus_xyz = [999,999,999]
        covar_largest_cl_mn  = [[999,999,999],[999,999,999],[999,999,999]]
        
        largest_cluster = np.zeros(8) # first entry: nr of points in cluster, second entry: density
       
        # TODO is there a fast way to first get all combinations and then loop only once? nesting loops is slow!
        for n, (params, comat) in enumerate(zip(lines_fit_, cov_mats_)):
            for m, (params1, comat1) in enumerate(zip(lines_fit_, cov_mats_)):
            
                if n == m:
                    continue
                    
                if n > m:
                    continue  
                    
                if angle_cut_ == True:
                    pars = params[0]
                    pars1 = params1[0] 
                else:
                    pars = params
                    pars1 = params1     
                    
                    
                if parameterisation == "angle":                                          
                    x0,y0,z0 = line_param_angles(1.0,pars)    # get a point on the current line
                    x1,y1,z1 = line_param_angles(5.0,pars)   # get another point on the line

                    x01,y01,z01 = line_param_angles(2.0,pars1)    # get a point on the current line 
                    x11,y11,z11 = line_param_angles(4.0,pars1)   # get another point on the line  
                    
                elif parameterisation == "p_d": 
                    x0,y0,z0 = line_param_new(1.0,pars)    # get a point on the current line
                    x1,y1,z1 = line_param_new(5.0,pars)   # get another point on the line

                    x01,y01,z01 = line_param_new(2.0,pars1)    # get a point on the current line 
                    x11,y11,z11 = line_param_new(4.0,pars1)   # get another point on the line               
                ####################################################
                # line one points and direction
                p0 = np.array([x0,y0,z0])
                p1 = np.array([x1,y1,z1])
                d1 = p1 - p0
                
                #print "p0", p0
                #print "p1", p1
                
                # line 2 points and direction
                p2 = np.array([x01,y01,z01])
                p3 = np.array([x11,y11,z11])
                d2 = p3 - p2

                #print "p2", p2
                #print "p3", p3
                
                # angle between tracks
                #angle = angle_3Dtracks(d1, d2)                
                #angle =  angle/(2.0*np.pi)*360.0

                # test for skewness: ###############################
                # the below stuff only works for skew lines
                # there shouldnt be any lines that are completely parallel though, so I dont do a cut on it
                c1 = p0 - p1
                c2 = p1 - p2
                c3 = p2 - p3
                skewness = 1.0/6.0*np.linalg.norm(np.linalg.det((c1,c2,c3)))
                #print "skewness", skewness
                if skewness < 1.0e-05:
                    continue
                
                ####################################################                                
                # nearstes points: ( wikipedia , articel for skew lines)
                # https://en.wikipedia.org/wiki/Skew_lines
                # np1 and np0 are the points of line0 and line1 that are closest to the other line
                n2_temp = np.cross(d1, d2) 
                n2 = np.cross(d2, n2_temp)
                
                n0_temp = np.cross(d2, d1) 
                n0 = np.cross(d1, n0_temp)                
                
                np0 = np.dot((p2 - p0),n2) 
                np0 = np0/(np.dot(d1,n2))
                np0 = np0 * d1 + p0
                
                np1 = np.dot((p0 - p2),n0)
                np1 = np1 / (np.dot(d2,n0))
                np1 = np1 * d2 + p2
                
                # mean of the two nearest points:           
                m_np = 0.5*(np0 + np1)
                
                # distance betw lines ###############################               
                n_vec = np.cross(d1, d2)
                n_vec = n_vec/np.linalg.norm(n_vec)                             
                dist = np.dot(n_vec, (p2 - p0))
                
                #help_vec = (np0 - np1)
                #dist2 = np.sqrt(help_vec[0]**2 + help_vec[1]**2 + help_vec[2]**2)
                if parameterisation == "angle":                                          
                    cov = get_midpoint_error_prop_covs_angles(pars, pars1, comat, comat1)                
                    
                elif parameterisation == "p_d":                        
                    cov = get_midpoint_error_prop_covs_2(pars, pars1, comat, comat1) 
                
                closest_points.append([np0,np1])
                dists.append(abs(dist)) # save the distances in an array
                mean_nps.append([m_np[0], m_np[1], m_np[2], n, m]) # save the mid point of the 2 points,
                                                                   # and the indices to find the corresponding lines later,
                covar_arr.append(cov) # plus the covariance matrix of error propagation
        
        mean_nps = np.array(mean_nps) 
        dists = np.array(dists)
        covar_arr = np.array(covar_arr)
        
        ############# 
        #total_covar_average = 999
        total_covar_average = np.average(covar_arr, axis = 0, weights = 1.0/(dists))
        total_covar_mean = np.mean(covar_arr, axis = 0)
        print "TOTAL COVARIANCE MATRIX, AVERAGE........."
        #print np.sqrt(np.diagonal(total_covar_average))
        print " "
        #############
        
        if len(mean_nps) > 0:
            print 'more than one line...'

        else:
            #mean_nps =[[999,999,999]]
            return mean_nps, 0, [], dists, 999, 999, [[999,999,999],[999,999,999],[999,999,999]], [[999,999,999],[999,999,999],[999,999,999]], [], [999], mean_clus_xyz, total_covar_mean,  covar_largest_cl_mn
        
        #######################################################################
        # look for clusters in the cloud of points
        print "Looking for clusters in the cloud of mid points of clostest points of line pairs..."
        #######################################################################
        
        plot_clusters = False
        
        if plot_clusters == True:
        
            fig_3d2 = plt.figure()
            ax_3d2 = fig_3d2.add_subplot(111, projection='3d')
            
            ax_3d2.set_xlim(-100, 100)
            ax_3d2.set_ylim(-100, 100)        
            ax_3d2.set_zlim(-100, 100) 
        
        xyz = mean_nps[:,0:3]
               
        #f.flush()
        
        """for xxx in xyz:
            for xx in xxx:
                f.write("%s " % xx)
            f.write("\n") 
            
        #"""
        
        #for i,j in zip(xyz, dists):
        #    print i,j
               
        from sklearn.cluster import DBSCAN
        from sklearn import metrics        
        # find clusters with high density
        # eps... max distance between points to belong to a cluster (I chose 2cm)
        # min_samples... minimum number of points in a cluster (I set it to half the number of tracks
        
        eps_ = 90.0
        
        if len(lines_fit_) <= 2:
            min_samples_ = 1.0
        else:
            min_samples_ = int(len(lines_fit_)*0.5) # math.ceil(len(lines_fit_)*0.5)
        
        
        min_samples_ = 1.0
        
        
        
        db = DBSCAN(eps=eps_, min_samples=min_samples_).fit(xyz)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        
        print "parameters of DBSCAN: eps:", eps_, " min_samples: ", min_samples_

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        print('Estimated number of clusters: %d' % n_clusters_)

        # array storing the mean point of the cluster
        mean_cluster_array = []
        
        if n_clusters_ < 1:
            return mean_nps, n_clusters_, mean_cluster_array, dists, closest_points, largest_cluster, [[999,999,999],[999,999,999],[999,999,999]], total_covar_average, [], [999], mean_clus_xyz, total_covar_mean,  covar_largest_cl_mn
                
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        counter_cl = 0          
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]

            class_member_mask = (labels == k)
            
            xyz2 = xyz[class_member_mask & core_samples_mask] # points in the cluster
            xyz_inds = mean_nps[class_member_mask & core_samples_mask][:,3:5]
            
            print "CLUSTER INDS", xyz_inds
            
            if len(xyz2) == 0:
                continue

            ####################################################################################               
            # distance between the 2 points that correspond to the mean point            
            #print "dist",dists
            #print "shape dist", dists.shape 
            
            curr_dist = dists[class_member_mask & core_samples_mask]

            # calculate weighted mean of point could on cluster
            av_clus_xyz = np.average(xyz2, weights = 1.0/(curr_dist), axis = 0)   
            mean_clus_xyz = np.mean(xyz2, axis = 0)          
            ################################################
            """print " ---------------------------------------------------", av_clus_xyz
            if (av_clus_xyz[0] > -40 and av_clus_xyz[0] < 40) and (av_clus_xyz[2] > -10 and av_clus_xyz[2] < 10):
                if counter_cl == 0:
                    f = open('asdaf_cluster_3.dat','a')
                    
                    for d in curr_dist:
                        f.write("%s " % d)
                        f.write("\n")
                    f.flush()"""

            counter_cl += 1
            
            # number of points
            nr_cl_p = len(xyz2)
            # density: number of points in a cluster devided by the number of clusters
            density_cl = np.sum(curr_dist)/(1.0*nr_cl_p)
            
            # covar matrices of cluster pts:            
            curr_covar_arr = covar_arr[class_member_mask & core_samples_mask]        
            #print curr_covar_arr.shape
      
            average_covar = np.average(curr_covar_arr, axis = 0, weights = 1.0/(curr_dist))
            mean_covar = np.mean(curr_covar_arr, axis = 0) #TODO
            #print "av", average_covar
            #print "sqrt of diagonal .........................................................", np.sqrt(np.diag(average_covar))
            
            print "number of points in cluster: ", len(xyz2), av_clus_xyz, "density:", density_cl
            #print "covariance matrix, error prop:", average_covar            
            #print "col:", tuple(col)
            
            if largest_cluster[0] == 0:
                largest_cluster[0] = nr_cl_p
                largest_cluster[1] = density_cl
                largest_cluster[2:5] = av_clus_xyz 
                largest_cluster[5:8] = mean_clus_xyz                
                covar_largest_cl =  average_covar 
                covar_largest_cl_mn = mean_covar
                largest_cluster_inds = xyz_inds
                largest_cl_dists = curr_dist
                
            elif nr_cl_p > largest_cluster[0]:
                largest_cluster[0] = nr_cl_p
                largest_cluster[1] = density_cl
                largest_cluster[2:5] = av_clus_xyz  
                largest_cluster[5:8] = mean_clus_xyz                 
                covar_largest_cl =  average_covar 
                covar_largest_cl_mn = mean_covar
                largest_cluster_inds = xyz_inds
                largest_cl_dists = curr_dist
                
            elif largest_cluster[0] == nr_cl_p and density_cl > largest_cluster[1]:
                largest_cluster[0] = nr_cl_p
                largest_cluster[1] = density_cl 
                largest_cluster[2:5] = av_clus_xyz 
                largest_cluster[5:8] = mean_clus_xyz                   
                covar_largest_cl =  average_covar
                covar_largest_cl_mn = mean_covar
                largest_cluster_inds = xyz_inds        
                largest_cl_dists = curr_dist
                
            # save, weighted mean of points in cluster, colour, number of points in cluster, density
            mean_cluster_array.append([av_clus_xyz, tuple(col), nr_cl_p, density_cl])
            ######################################################################################
            
            if plot_clusters == True:
                ax_3d2.plot(xyz2[:, 0], xyz2[:, 1],  xyz2[:,2],'o', markerfacecolor=tuple(col),
                        markersize=6)

                xyz2 = xyz[class_member_mask & ~core_samples_mask] # "noise" hits, not belonging to a cluster
                ax_3d2.plot(xyz2[:, 0], xyz2[:, 1],  xyz2[:,2], 'o', markerfacecolor=tuple(col),
                         markeredgecolor='k', markersize=6)

                plt.title('Estimated number of clusters: %d' % n_clusters_)
                
        return mean_nps, n_clusters_, mean_cluster_array, dists, closest_points, largest_cluster, covar_largest_cl, total_covar_average, largest_cluster_inds, largest_cl_dists, mean_clus_xyz, total_covar_mean,  covar_largest_cl_mn #p0,p1,p2,p3, np0, np1
######################################################################################################################################################################### 
#########################################################################################################################################################################  
#########################################################################################################################################################################
def angle_3Dtracks(v1, v2):
# v1 is your firsr vector
# v2 is your second vector
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

    angle2 =  2 * np.pi - angle
    #print "angle", angle/(2.0*np.pi)*360.0, angle2/(2.0*np.pi)*360.0

    return min(angle,angle2)

######################################################################################################################################################################### 
######################################################################################################################################################################### 
######################################################################################################################################################################### 

def ndmesh(*xi,**kwargs):
    if len(xi) < 2:
        msg = 'meshgrid() takes 2 or more arguments (%d given)' % int(len(xi) > 0)
        raise ValueError(msg)

    args = np.atleast_1d(*xi)
    ndim = len(args)
    copy_ = kwargs.get('copy', True)

    s0 = (1,) * ndim
    output = [x.reshape(s0[:i] + (-1,) + s0[i + 1::]) for i, x in enumerate(args)]

    shape = [x.size for x in output]

    # Return the full N-D matrix (not only the 1-D vector)
    if copy_:
        mult_fact = np.ones(shape, dtype=int)
        return [x * mult_fact for x in output]
    else:
        return np.broadcast_arrays(*output)


######################################################################################################################################################################### 
######################################################################################################################################################################### 
def get_xy_of_line_new(z,p):
    
    #x = p[0] + p[1]*z
    #y = p[1] + p[3]*z
    
    y = (z - p[2])/p[3]
    x = p[0] + p[1]*y
        
    return x,y
######################################################################################################################################################################### 
def get_xy_of_line(z,p):
    t = (z - p[2]) / p[5]
    
    x = p[0] + p[3]*t
    y = p[1] + p[4]*t
        
    return x,y
#########################################################################################################################################################################
def line_param_new(t,p): # parameter form of a straight line 
    #x = p[0] + p[1]*t
    #y = p[2] + p[3]*t
    #z = t
    ####
    x = p[0] + p[2]*t    
    y = t    
    z = p[1] + p[3]*t
    return x,y,z 
#########################################################################################################################################################################
def line_param_angles2(t,p): # parameter form of a straight line 
    #x = p[0] + p[1]*t
    #y = p[2] + p[3]*t
    #z = t
    ####
    x = p[0] + math.cos(p[2])*math.sin(p[3])*t    
    y = p[4] +  math.sin(p[2])*math.sin(p[3])*t    
    z = p[1] + math.cos(p[3])*t
    return x,y,z     
#########################################################################################################################################################################
def line_param_angles(t,p): # parameter form of a straight line 
    #x = p[0] + p[1]*t
    #y = p[2] + p[3]*t
    #z = t
    ####
    x = p[0] + math.cos(p[2])*math.sin(p[3])*t    
    y =        math.sin(p[2])*math.sin(p[3])*t    
    z = p[1] + math.cos(p[3])*t
    return x,y,z 
#########################################################################################################################################################################
def get_xy_of_line_angles(z,p):
    
    t = (z - p[1])/math.cos(p[2])*math.sin(p[3])
    
    y =        math.sin(p[2])*math.sin(p[3])*t    
    z = p[1] + math.cos(p[3])*t
        
    return y,z
#########################################################################################################################################################################
def get_zy_of_line_angles(x,p):
    
    t = (x - p[0])/math.cos(p[3])
    x = p[0] + math.cos(p[2])*math.sin(p[3])*t    
    y =        math.sin(p[2])*math.sin(p[3])*t        
    
        
    return x,y
#########################################################################################################################################################################
def get_xy_of_line_angles2(z,p):
    
    t = (z - p[1])/math.cos(p[3])
    x = p[0] + math.cos(p[2])*math.sin(p[3])*t    
    y = p[4] + math.sin(p[2])*math.sin(p[3])*t        
    
        
    return x,y    
#########################################################################################################################################################################    
def dist_p_line_angles2(x,y,z,p): # function to calculate the distance in 3D of point [x,y,z] to the line with parameters p. returns distance^2
    # distance line point is D= | (xp-x0) cross  ux | 
    # where ux is direction of line and x0 is a point on the line (like t = 0) 
    xp = np.array([x,y,z])

    x0 = np.array([p[0], p[4], p[1]])

    x1 = np.array([p[0] + math.cos(p[2])*math.sin(p[3]), p[4] + math.sin(p[2])*math.sin(p[3]),  p[1] + math.cos(p[3]) ])

    uabs = math.sqrt( (x1[0]-x0[0])*(x1[0]-x0[0]) + (x1[1]-x0[1])*(x1[1]-x0[1]) + (x1[2]-x0[2])*(x1[2]-x0[2])) #np.absolute(x1 - x0) #
    u = (x1-x0)/uabs
    tempvv = np.array([xp[0]-x0[0], xp[1]-x0[1], xp[2]-x0[2]])
    crossp = np.cross(tempvv,u)
    #d2 = np.dot(crossp,crossp)    
   
    return crossp # distance squared!!
#########################################################################################################################################################################    
def dist_p_line_angles(x,y,z,paras): # function to calculate the distance in 3D of point [x,y,z] to the line with parameters p. 
    # distance line point is D= | (xp-x0) cross  ux | 
    # where ux is direction of line and x0 is a point on the line (like t = 0) 
 
    xp = np.array([x,y,z])
    
    
    theta = paras[3]
    phi = paras[2]
    px = paras[0]
    pz = paras[1]
    
    #spherical_matrix = np.array([[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi), math.cos(theta)], 
    #                             [math.cos(theta)*math.cos(phi), math.cos(theta)*math.sin(phi), - math.sin(theta)],
    #                             [- math.sin(phi), math.cos(phi), 0 ]   ])
                                     
    #spherical_matrix_inv = np.array([[math.sin(theta)*math.cos(phi),math.cos(theta)*math.cos(phi), - math.sin(phi)], 
    #                             [math.sin(theta)*math.sin(phi), math.cos(theta)*math.sin(phi), math.cos(phi)],
    #                             [math.cos(theta), - math.sin(theta), 0 ]   ])     
    

    x0 = np.array([px, 0.0, pz])
    
    #unit_cart = np.dot(spherical_matrix_inv,unit_spher)

    x1 = np.array([px + math.cos(phi)*math.sin(theta), math.sin(phi)*math.sin(theta),  pz + math.cos(theta) ])


    uabs = math.sqrt( (x1[0]-x0[0])*(x1[0]-x0[0]) + (x1[1]-x0[1])*(x1[1]-x0[1]) + (x1[2]-x0[2])*(x1[2]-x0[2])) #np.absolute(x1 - x0) #
    u = (x1-x0)/uabs
    tempvv = np.array([xp[0]-x0[0], xp[1]-x0[1], xp[2]-x0[2]])
    crossp = np.cross(tempvv,u)
    #d2 = np.dot(crossp,crossp) 
    
    #print "--------------------------------------", np.sqrt(d2)   
   
    return crossp # distance squared!!
#########################################################################################################################################################################
#########################################################################################################################################################################    
def dist_p_line_new(x,y,z,p): # function to calculate the distance in 3D of point [x,y,z] to the line with parameters p. returns distance^2
    # distance line point is D= | (xp-x0) cross  ux | 
    # where ux is direction of line and x0 is a point on the line (like t = 0) 
    xp = np.array([x,y,z])

    x0 = np.array([p[0], 0.0, p[1] ])

    x1 = np.array([p[0] + p[2], 1.0,  p[1] + p[3] ])

    uabs = math.sqrt( (x1[0]-x0[0])*(x1[0]-x0[0]) + (x1[1]-x0[1])*(x1[1]-x0[1]) + (x1[2]-x0[2])*(x1[2]-x0[2])) #np.absolute(x1 - x0) #
    u = (x1-x0)/uabs
    tempvv = np.array([xp[0]-x0[0], xp[1]-x0[1], xp[2]-x0[2]])
    crossp = np.cross(tempvv,u)
    #d2 = np.dot(crossp,crossp)       
    return crossp # distance squared!!
#########################################################################################################################################################################
def line_param(t,p): # parameter form of a straight line 
    x = p[0] + p[3]*t
    y = p[1] + p[4]*t
    z = p[2] + p[5]*t
    return x,y,z
######################################################################################################################################################################### 
#########################################################################################################################################################################
#########################################################################################################################################################################
def track_selection_fitting_3D_4layers_root_angles(point_collections):

    print " "
    print "-------------------------------- 3D TRACK SELECTION START --------------------------------"
    final_lines_params = [] # store parameter of fitted lines, size Nr_tracks x Nr_params
    chisq = [] # list to save chi2 values
    conv_status = [] # integer, did the fit converge?
    cov_matrices = [] # array containing the covariance matrices
    covar_mat_status = [] # status of covariance matrices

    for n, ps in enumerate(point_collections):
                #print " "
                points = []
                errors = []
                
                if len(ps) < 3:
                    continue
                
                for ui in ps:      
                    points.append([ui[2], ui[4], ui[0]]) # [x, y, z]
                    errors.append([ui[3], ui[5], ui[1]]) # errors of x,y,z
                              
                samplesize = len(points) # number of points to fit
                                
                #####################################################################################
                # finding good start values:
                #
                # it can be that depending on the event, different start values lead to better convergence of the fit
                # below i calculate from two hits of the track candidate, the direction vector of the straight line,
                # and from that the two parameter angles theta and phi.
                # for the start parameters px and pz, zero is often a good start value
                # I also calculate the values of x and z of the point (x,0,z), this can also lead to good convergence
                # for certain events. all in all, I think zero is a good start, if there are not many events like in antihydrogen
                # runs, the candidates can be also looked at individually.
                # the finding of start values can cerntainly be improved.

                # track, inner fibre point, if available. if not, then the inner bar hit is used
                if ps[1][7] != 999: 
                    start_val_pt = np.array([ps[1][2], ps[1][4], ps[1][0]]) 
                else:
                    start_val_pt = np.array([ps[0][2], ps[0][4], ps[0][0]]) 
                #print "start vals:", start_val_pt

                # track, outer most point of track 
                start_val_pt2 = np.array([ps[len(ps)-1][2], ps[len(ps)-1][4], ps[len(ps)-1][0]])

                start_val_dir = start_val_pt - start_val_pt2                
                #x = p[0] + math.cos(p[2])*math.sin(p[3])*t    
                #y =        math.sin(p[2])*math.sin(p[3])*t    
                #z = p[1] + math.cos(p[3])*t
                
                rho_start = math.sqrt(start_val_dir[0]**2 + start_val_dir[1]**2 + start_val_dir[2]**2)
                theta_start = math.acos(start_val_dir[2]/rho_start) 
                phi_start = math.atan(start_val_dir[1]/start_val_dir[0])     
                #print "dir ", start_val_dir
                #print "dir check",  rho_start*math.cos(phi_start)*math.sin(theta_start), rho_start*math.sin(phi_start)*math.sin(theta_start), rho_start*math.cos(theta_start)
                                  
                # find point with y = 0 for start values of fit:
                #x = p[0] + p[3]*t
                #y = p[1] + p[4]*t
                #z = p[2] + p[5]*t
                t_ = - start_val_pt[1]/start_val_dir[1]
                x_start = start_val_pt[0] + t_*start_val_dir[0]
                z_start = start_val_pt[2] + t_*start_val_dir[2]
                      
                #spherical_matrix = np.array([[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi), math.cos(theta)], 
                #                 [math.cos(theta)*math.cos(phi), math.cos(theta)*math.sin(phi), - math.sin(theta)],
                #                 [- math.sin(phi), math.cos(phi), 0 ]   ])
                
                # set the startvalues for the fit                                               
                pars = np.array([x_start, z_start, phi_start, theta_start]) # start values
                                
                pars = np.array([0.0, 0.0, phi_start, theta_start]) # start values                             
                #pars = np.array([x_start,z_start,0.0,0.0]) # start values
                #pars = np.array([0.0,0.0,0.0,0.0]) # start values
                
                step = [0.001, 0.001, 0.001, 0.001] # step size for the minimiser 
                minimizer = ROOT.Minuit2.Minuit2Minimizer(0) # here are several options possible!
                #EMinimizerType algoType = kMigrad;
                #if (algoname == "simplex")   algoType = kSimplex;
                #if (algoname == "minimize" ) algoType = kCombined;
                #if (algoname == "scan" )     algoType = kScan;
                #if (algoname == "fumili" ) algoType = kFumili; 
                               
                # print errors etc during the minimisation
                #minimizer.SetPrintLevel(3) 
               
                errorfct = chi2fct3d_angles() # define the function to minimise
                
                errorfct.SetVals(samplesize, points, errors)
                minimizer.SetFunction(errorfct)

                minimizer.SetMaxFunctionCalls(500)
                minimizer.SetMaxIterations(500)
                minimizer.SetTolerance(0.001)
                        
                # for historical reasons (compatibility with the old Fortran version of Minuit)
                # the Migrad minimisation will stop when the
                # computed EDM is less than 0.002 * tolerance * UPERROR.
                # UPERROR = 1 (68% errors in a chi2 fit)
                # edm ..... expected distance reached from the minimum
                                
                #adding variables without limited range:
                minimizer.SetVariable(0,"px", pars[0], step[0])
                minimizer.SetVariable(1,"pz", pars[1], step[1])
          
                # adding them with a limited range:
                # the limited range sometimes doesnt let the fit converge...
                # so maybe it's better to use the unlimted case.
                #minimizer.SetLimitedVariable(0, "px", pars[0], step[0], -2000.0, 2000.0)
                #minimizer.SetLimitedVariable(1, "pz", pars[1], step[1], -2000.0, 2000.0)#
               
                minimizer.SetVariable(2,"phi", pars[2], step[2])
                minimizer.SetVariable(3,"theta", pars[3], step[3]) 
                #minimizer.SetLimitedVariable(2, "phi", pars[2], step[2], - 2.0*math.pi, 2.0*math.pi)
                #minimizer.SetLimitedVariable(3, "theta", pars[3], step[3], - math.pi, math.pi) 
                         
                retval_min = minimizer.Minimize() # start the minimiser
                
                minval = minimizer.MinValue() # minvalue of the function to minimise (chi2fct3d_angles() in this case)

                # print the fit results after the fit - see whether it converges etc.
                #minimizer.PrintResults()
                
                stat = minimizer.Status()
                #status 0 means that the errors of the fit are ok.
                #status = 1    : Covariance was made pos defined.  status 1 means, you cant trust the errorbars (internet info)
                #status = 2    : Hesse is invalid
                #status = 3    : Edm is above max
                #status = 4    : Reached call limit
                #status = 5    : Any other failure
                
                covar_stat = minimizer.CovMatrixStatus()
                #print "-------------------------------------------- status of cov matrix: ", covar_stat, "---------------------------------------------"                               
                #return the status of the covariance matrix status = -1 : not available (inversion failed or Hesse failed)
                #status = 0 : available but not positive defined
                #status = 1 : covariance only approximate
                #status = 2 : full matrix but forced pos def
                #status = 3 : full accurate matrix
                temp_cov = []
                
                """if stat == 0:
                    calls_file = open("saved_number_calles.dat","a+")
                    calls_file.write("%d " % minimizer.NIterations()) 
                    calls_file.write("%d " % minimizer.NCalls()) 
                    calls_file.write("\n") """
                               
                ret_params = minimizer.X() # parameters of the fit
                ret_errors = minimizer.Errors() # errors of the fit
                final_params = [] # list to save all the parameters of the track
                       
                for i in range(0,pars.size):
                    final_params.append(ret_params[i])
                    for j in range(0,pars.size):
                        #print minimizer.CovMatrix(i,j)
                        temp_cov.append(minimizer.CovMatrix(i,j))
                        
                cov_matrices.append(temp_cov)
                covar_mat_status.append(covar_stat)
                
                final_lines_params.append(final_params)                   
                chisq.append(minval)
                conv_status.append(stat)
                
    return final_lines_params, np.array(chisq), np.array(conv_status), np.array(cov_matrices), np.array(covar_mat_status)       
##############################################################################################################
######################################################################################################################################################################### 
######################################################################################################################################################################### 
def track_selection_fitting_3D_4layers_root_new(point_collections):

    print " "
    print "-------------------------------- 3D TRACK SELECTION START --------------------------------"
    final_lines_params = [] # store parameter of fitted lines, size Nr_tracks x Nr_params
    chisq = [] # list to save chi2 values
    conv_status = [] # integer, did the fit converge?
    cov_matrices = [] # array containing the covariance matrices
    covar_mat_status = [] # status of covariance matrices
    for n, ps in enumerate(point_collections):
                #print " "
                points = []
                errors = []
                
                if len(ps) < 3:
                    continue
                
                for ui in ps:      
                    #print "ui", ui
                    points.append([ui[2], ui[4],ui[0]]) # [x, y, z]
                    errors.append([ui[3], ui[5], ui[1]]) # errors of x,y,z
                    #points.append([i[2], i[4],i[0]]) # [x, y, z]
                #points.append([0, 0, 0]) # [x, y, z]
                #errors.append([50, 50, 50])
                #print "inner:", i[0][2], i[0][4],i[0][0]
                #print "i err:", i[0][3], i[0][5], i[0][1]
                #print "outer:", i[0][2], i[0][4],i[0][0]
                #print "o err:", i[1][3], i[1][5], i[1][1]
                               
                samplesize = len(points) # number of points to fit
                
                if ps[1][7] != 999: 
                    start_val_pt = np.array([ps[1][2], ps[1][4], ps[1][0]]) 
                else:
                    start_val_pt = np.array([ps[0][2], ps[0][4], ps[0][0]]) 
                #print "start vals:", start_val_pt

                # track, outer most point of track 
                start_val_pt2 = np.array([ps[len(ps)-1][2], ps[len(ps)-1][4], ps[len(ps)-1][0]])

                start_val_dir = start_val_pt - start_val_pt2               
                
                
                #t_for_pt_start = - start_val_pt[1] / start_val_dir[1]
                x_start = start_val_pt[0] #+ t_for_pt_start*start_val_dir[0]
                z_start = start_val_pt[2] #+ t_for_pt_start*start_val_dir[2]
                
                
                #norm_dir = np.linalg.norm(start_val_dir)
                #norm_pt = np.linalg.norm(start_val_pt)
                start_val_dir = start_val_dir/start_val_dir[1]
                
                pars = np.array([x_start, z_start, start_val_dir[0], start_val_dir[2]]) #"""
                pars = np.array([0, 0, start_val_dir[0], start_val_dir[2]]) 
                
                #pars = np.array([-1.0,-1.0,-1.0,-1.0]) # start values
                #pars = np.array([0.0,0.0,0.0,0.0,0.0,0.0]) # start values
                #pars = np.array([100.0,100.0,100.0,100.0,100.0,100.0]) # start values
                
                step = [0.01, 0.01, 0.01, 0.01] # step size for the minimiser 
                minimizer = ROOT.Minuit2.Minuit2Minimizer(1) # here are several options possible!
                
                #minimizer.SetPrintLevel(10)
                
                #EMinimizerType algoType = kMigrad;
                #if (algoname == "simplex")   algoType = kSimplex;
                #if (algoname == "minimize" ) algoType = kCombined;
                #if (algoname == "scan" )     algoType = kScan;
                #if (algoname == "fumili" ) algoType = kFumili;                
                
                errorfct = chi2fct3d_new() # define the function to minimise
                
                errorfct.SetVals(samplesize, points, errors)
                minimizer.SetFunction(errorfct)
                #minimizer.SetDefaultPrintLevel(2)

                minimizer.SetMaxFunctionCalls(500)
                minimizer.SetMaxIterations(700)
                minimizer.SetTolerance(0.001)
                # for historical reasons (compatibility with the old Fortran version of Minuit)
                # the Migrad minimisation will stop when the
                # computed EDM is less than 0.002 * tolerance * UPERROR.
                # UPERROR = 1 (68% errors in a chi2 fit)
                # edm ..... expected distance reached from the minimum
                
                

                minimizer.SetVariable(0,"par0", pars[0], step[0])
                minimizer.SetVariable(1,"par1", pars[1], step[1])
                minimizer.SetVariable(2,"par2", pars[2], step[2])
                minimizer.SetVariable(3,"par3", pars[3], step[3])
                
              
                retval_min = minimizer.Minimize()
                
                minval = minimizer.MinValue()

                # print the fit results after the fit - see whether it converges etc.
                #minimizer.PrintResults()
                
                stat = minimizer.Status() # status 0 means that the errors of the fit are ok.
                #status = 1    : Covariance was made pos defined.  status 1 means, you cant trust the errorbars (internet info)
                #status = 2    : Hesse is invalid
                #status = 3    : Edm is above max
                #status = 4    : Reached call limit
                #status = 5    : Any other failure
                
                covar_stat = minimizer.CovMatrixStatus()
                #print "-------------------------------------------- status of cov matrix: ", covar_stat, "-------------------------"                
                temp_cov = []
                #return the status of the covariance matrix status = -1 : not available (inversion failed or Hesse failed)
                #status = 0 : available but not positive defined
                #status = 1 : covariance only approximate
                #status = 2 : full matrix but forced pos def
                #status = 3 : full accurate matrix
                
                
                ret_params = minimizer.X() # parameters of the fit

                ret_errors = minimizer.Errors() # errors of the fit
                final_params = [] # list to save all the parameters of the track
                       
                for i in range(0,pars.size):
                    final_params.append(ret_params[i])
                    for j in range(0,pars.size):
                        temp_cov.append(minimizer.CovMatrix(i,j))
                        
                cov_matrices.append(temp_cov)
                covar_mat_status.append(covar_stat)
                
                final_lines_params.append(final_params)                   
                chisq.append(minval)
                conv_status.append(stat)
                
    return final_lines_params, np.array(chisq), np.array(conv_status), np.array(cov_matrices), np.array(covar_mat_status)       
##############################################################################################################
######################################################################################################################################################################### 
##############################################################################################################
##############################################################################################################
from scipy.optimize import fmin, leastsq
def line_points( s, p, t ):
    return [ s * t[0] + p[0], s * t[1] + p[1], s * t[2] + p[2] ]

def weighted_dist( s, p, t, xVec, sigmaVec ):
    q = line_points( s, p, t )
    d  = ( q[0] - xVec[0] )**2 / sigmaVec[0]**2
    d += ( q[1] - xVec[1] )**2 / sigmaVec[1]**2
    d += ( q[2] - xVec[2] )**2 / sigmaVec[2]**2
    return np.sqrt( d )

def weighted_od( p, t, xVec, sigmaVec ):
    f = lambda s: weighted_dist( s, p, t, xVec, sigmaVec )
    sol = fmin( f, 0, disp=False )
    d = weighted_dist( sol[0], p, t, xVec, sigmaVec )
    return d

def residuals( params, data, sigmas ): ###data of type [ allx, ally, allz], sigma of type [allsx, allsy, allsz]
    px, py, pz, tx, ty, tz = params
    out = list()
    for x0, y0, z0, sx, sy, sz in zip( *( data + sigmas ) ):
        out += [weighted_od( [ py, py, pz ], [ tx, ty, tz ], [ x0, y0, z0 ], [ sx, sy, sz ] ) ]
    print sum(out)
    return out
##############################################################################################################   
# FIXME this function doesnt work... it was just a try to implement the fitting purely in python
def track_selection_fitting_3D_4layers_python(point_collections):
    print " "
    print "-------------------------------- 3D TRACK SELECTION START --------------------------------"
    final_lines_params = [] # store parameter of fitted lines, size Nr_tracks x Nr_params
    chisq = [] # list to save chi2 values
    conv_status = [] # integer, did the fit converge?
    cov_matrices = [] # array containing the covariance matrices
    covar_mat_status = [] # status of covariance matrices
    
    for n, ps in enumerate(point_collections):
                #print " "
                points = []
                errors = []
                
                if len(ps) < 3:
                    continue
                
                for ui in ps:      
                    #print ui
                    points.append([ui[2], ui[4],ui[0]]) # [x, y, z]
                    errors.append([ui[3], ui[5], ui[1]]) # errors of x,y,z
                    #points.append([i[2], i[4],i[0]]) # [x, y, z]
                
                samplesize = len(points) # number of points to fit
                
                points = np.array(points)
                errors = np.array(errors)
                
                # track one, fibre point
                start_val_pt = np.array([ps[1][2], ps[1][4], ps[1][0]])                 
                # track one, another point
                start_val_pt2 = np.array([ps[2][2], ps[2][4], ps[2][0]])

                start_val_dir = start_val_pt2 - start_val_pt
              
                pars = np.array([[start_val_pt[0], start_val_pt[1], start_val_pt[2]],  [start_val_dir[0], start_val_dir[1], start_val_dir[2]]]).reshape((2,3)) # start values
                                
                pts_and_errs = np.append(points, errors)
                print pts_and_errs.shape 
                
                data = None
                #residual(pars, pts_and_errs, data = None)
                
                ###################################################
                
                from random import random
              
                print points
                xData =  points[:,0]
                yData =  points[:,1]
                zData =  points[:,2]
                
                xErr =  errors[:,0]
                yErr =  errors[:,1]
                zErr =  errors[:,2]                
                
                print xData
                print yData

                xyzData = [ xData, yData, zData ]
                sssData = [ xErr, yErr, zErr ]

                residuals( [ 1, 1, 3, -1,-3,.8 ],  xyzData, sssData )
                myFit, err = leastsq(residuals, [ 1, 1, 2 , -1, -2, -1 ], args=( xyzData, sssData ) )
                print myFit

                fitP = myFit[:3]
                fitT = myFit[3:]
                fitTN= np.linalg.norm( fitT )
                fitT = [ fitT[0] / fitTN, fitT[1] / fitTN, fitT[2] / fitTN ]
                fitLineList = [ line_points( s, fitP, fitT ) for s in sList ] 

                ax = m3d.Axes3D(plt.figure() )
                ax.plot( *zip(*lineList) )
                ax.plot( *zip(*fitLineList) )
                ax.scatter3D( xData, yData, zData )
                plt.show()
                
                
                
                
                ######################################################
                print "start paras", pars 
                result = scipy.optimize.leastsq(residual, pars, pts_and_errs)
                
                fittedParams = result[0]
                
                print fittedParams
                
                exit(0)
                #cov_x : ndarray
                #Uses the fjac and ipvt optional outputs to construct an estimate of the jacobian around the solution.
                #None if a singular matrix encountered (indicates very flat curvature in some direction).
                #This matrix must be multiplied by the residual variance to get the covariance of the parameter estimates - see curve_fit. """
                
                """final_params = [] # list to save all the parameters of the track
                       
                for i in range(0,pars.size):
                    final_params.append(ret_params[i])
                    for j in range(0,pars.size):
                        temp_cov.append(minimizer.CovMatrix(i,j))
                        
                cov_matrices.append(temp_cov)
                covar_mat_status.append(covar_stat)
                
                final_lines_params.append(final_params)                   
                chisq.append(minval)
                conv_status.append(stat)"""
                
    return final_lines_params, np.array(chisq), np.array(conv_status), np.array(cov_matrices), np.array(covar_mat_status)      
###################################################################################################################################################
###################################################################################################################################################        
##############################################################################################################
##############################################################################################################    
def dist_p_line(x,y,z,p): # function to calculate the distance in 3D of point [x,y,z] to the line with parameters p. returns distance^2
    # distance line point is D= | (xp-x0) cross  ux | 
    # where ux is direction of line and x0 is a point on the line (like t = 0) 
    xp = np.array([x,y,z])
    
    #print "xp", xp

    #4 -> 2
    #1 -> 3
    #2 -> 1
    #3 -> 4
    x0 = np.array([p[0], p[1], p[2] ])
    
    #print "x0", x0
    
    
    
    x1 = np.array([p[0] + p[3], p[1] + p[4], p[2] + p[5] ])
    
    
    #print "x1", x1
    
    
    uabs = math.sqrt( (x1[0]-x0[0])*(x1[0]-x0[0]) + (x1[1]-x0[1])*(x1[1]-x0[1]) + (x1[2]-x0[2])*(x1[2]-x0[2])) #np.absolute(x1 - x0) #
    u = (x1-x0)/uabs
    tempvv = np.array([xp[0]-x0[0], xp[1]-x0[1], xp[2]-x0[2]])
    crossp = np.cross(tempvv,u)
    d2 = np.dot(crossp,crossp)    
   
    return d2 # distance squared!!
############################################################################################################## 
##############################################################################################################     
class chi2fct3d( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(12)
        self.errors = np.zeros(12)

    def NDim( self ): # total number of variables
        #print 'PYTHON NDim called'
        return 6.0

    def SetVals(self, samplesizet, pointst, errorst):
        self.samplesize = samplesizet
        self.points = pointst
        self.errors = errorst
        #print "pts", self.points
        #print "err", self.errors
    
    def DoEval( self, pars):
        #print "hai!"   
        chisq = []   
       
        for i in range(0,self.samplesize):
                 
            errs = 0.0
            
            errsx = self.errors[i][0]**2
            errsy = self.errors[i][1]**2
            errsz = self.errors[i][2]**2

            xvalue = self.points[i][0]
            yvalue = self.points[i][1]
            zvalue = self.points[i][2]
            distx, disty, distz = dist_p_line_new(xvalue, yvalue, zvalue, pars) # distance squared of line to point i
            distx = distx**2
            disty = disty**2
            distz = distz**2            
            chisq.append(distx/errsx + disty/errsy + distz/errsz) # add chisq of point i
            

        chisq = np.sum(np.array(chisq))

        return chisq 
###########################################################################################################################################
        
class chi2fct3d_angles2( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(12)
        self.errors = np.zeros(12)

    def NDim( self ): # total number of variables
        #print 'PYTHON NDim called'
        return 5.0

    def SetVals(self, samplesizet, pointst, errorst):
        self.samplesize = samplesizet
        self.points = pointst
        self.errors = errorst
        #print "pts", self.points
        #print "err", self.errors
    
    def DoEval( self, pars):
        #print "hai!"   
        chisq = []   
       
        for i in range(0,self.samplesize):
                 
            # sum of errors of point i     
            errs = 0.0
            
            errsx = self.errors[i][0]**2
            errsy = self.errors[i][1]**2
            errsz = self.errors[i][2]**2

            xvalue = self.points[i][0]
            yvalue = self.points[i][1]
            zvalue = self.points[i][2]
            distx, disty, distz = dist_p_line_angles2(xvalue, yvalue, zvalue, pars) # distance squared of line to point i
            distx = distx**2
            disty = disty**2
            distz = distz**2            
            chisq.append(distx/errsx + disty/errsy + distz/errsz) # add chisq of point i
            #chisq.append((distx + disty + distz)/(errsx + errsy + errsz)) # add chisq of point i
            

        chisq = np.sum(np.array(chisq))

        return chisq 
##############################################################################################################  
##############################################################################################################        
        
##############################################################################################################
##############################################################################################################
class chi2fct3d_angles( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(12)
        self.errors = np.zeros(12)

    def NDim( self ): # total number of variables
        #print 'PYTHON NDim called'
        return 4.0

    def SetVals(self, samplesizet, pointst, errorst):
        self.samplesize = samplesizet
        self.points = pointst
        self.errors = errorst
        #print "pts", self.points
        #print "err", self.errors
    
    def DoEval( self, pars):
        #print "hai!"   
        chisq = []   
        
        #print "points:"#, self.points
        
        #for j in self.points:
        #    print j
      
        for i in range(0,self.samplesize):
        
            #print "..................",i 
            
            errsx = self.errors[i][0]**2
            errsy = self.errors[i][1]**2
            errsz = self.errors[i][2]**2
            
            #print errsx, errsy, errsz

            xvalue = self.points[i][0]
            yvalue = self.points[i][1]
            zvalue = self.points[i][2]
            distx, disty, distz = dist_p_line_angles(xvalue, yvalue, zvalue, pars) # distance squared of line to point i
            
            #print distx, disty, distz
            
            distx = distx**2
            disty = disty**2
            distz = distz**2       
            #print distx, disty, distz
            
            #print "_______________________________________________________________________________________", math.sqrt(distx + disty + distz)     
            #chisq.append(np.sqrt(distx/errsx + disty/errsy + distz/errsz)) # add chisq of point i
            chisq.append((distx + disty + distz)/(errsx + errsy + errsz)) # add chisq of point i
            #print (distx + disty + distz)/(errsx + errsy + errsz)

        chisq = np.sum(np.array(chisq))
        
        #print chisq

        return chisq 
##############################################################################################################  
##############################################################################################################
class chi2fct3d_new( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(12)
        self.errors = np.zeros(12)

    def NDim( self ): # total number of variables
        #print 'PYTHON NDim called'
        return 4.0

    def SetVals(self, samplesizet, pointst, errorst):
        self.samplesize = samplesizet
        self.points = pointst
        self.errors = errorst
        #print "pts", self.points
        #print "err", self.errors
    
    def DoEval( self, pars):
        #print "hai!"   
        chisq = []   
       
        for i in range(0,self.samplesize):
                 
            # sum of errors of point i     
            errs = 0.0
            
            errsx = self.errors[i][0]**2
            errsy = self.errors[i][1]**2
            errsz = self.errors[i][2]**2

            xvalue = self.points[i][0]
            yvalue = self.points[i][1]
            zvalue = self.points[i][2]
            distx, disty, distz = dist_p_line_new(xvalue, yvalue, zvalue, pars) # distance squared of line to point i
            distx = distx**2
            disty = disty**2
            distz = distz**2            
            chisq.append((distx + disty + distz)/(errsx + errsy + errsz)) # add chisq of point i
            

        chisq = np.sum(np.array(chisq))

        return chisq 
############################################################################################################## 
##############################################################################################################
def get_z_from_line(p,x,y):
    z1 = p[2] + (x - p[0])/p[3]*p[5]
    z2 = p[2] + (y - p[1])/p[4]*p[5]
    #4 -> 2
    #1 -> 3
    #2 -> 1
    #3 -> 4
    print "************** ", z1, "*********************", z2, "**********************"
    #return (z2+z1)/(2.0)
    return z2
######################################################################################################################################################################### 
######################################################################################################################################################################### 
##############################################################################################################
class chi2fct2d( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(3)
        self.errors = np.zeros(3)

    def NDim( self ):
        #print 'PYTHON NDim called'
        return 2

    def SetVals(self, samplesizet, pointst, errorst):
        self.samplesize = samplesizet
        self.points = pointst
        self.errors = errorst
    
    def DoEval( self, pars):
        #print "hai!"   
        chisq = 0.0    
        errsum = 0.0
        #print  self.errors
        #print self.points
        #for i in range(0,samplesize):
        #    errsum = errsum + errors[i][0] + errors[i][1] + errors[i][2] 
        errsum = np.sum(self.errors[:])
        
        #y = -p[0]/p[1]*x - p[2]/p[1]  
        for i in range(0,self.samplesize):
            chisq = chisq + (self.points[i][1] - pars[0] - pars[1]*self.points[i][0])**2/(self.errors[i][1]**2 + self.errors[i][0]**2*pars[1]**2)

        return chisq  
############################################################################################################## 
##############################################################################################################
def line_2d(x,p):
    y = p[0]*x + p[1]    
    return y  
############################################################################################################################################################################ 
############################################################################################################################################################################ 

def track_finding_circle(pmapI, cmap, currI, currO, circle, BGOXPos,BGOYPos, plot_trackfinding, innerzPos, outerzPos): #, BGOXPos, BGOYPos):  # with density of lines added
        ### First: Find Clusters ###
        #print "####################"
        #print "inner bars hit: ", cmap[:nI,1] # corner map, all rects with corner points plus nr
        #print "outer bars hit: ", cmap[nI:,1]     
        #plt.show()        
        #print "pmapI ", pmapI
        
        
        
        #print "Outer CFs "
        #print outerCFs

        nI = len(currI)
        nO = len(currO)
       
        cmapI = cmap[:nI,:]
        cmapO = cmap[nI:,:] 
        
        # use this for no unification of rects
        yesnoI, cI = find_consecutive_nr(currI)  
        yesnoO, cO = find_consecutive_nr(currO)
        
        merge_cl_i = []
        merge_cl_o = []
        #####################################################################################################################
        #####################################################################################################################         
        ynI, clusI = find_consecutive_nr_singletoo(currI)
        ynO, clusO = find_consecutive_nr_singletoo(currO)
               
        #print "Inner cluster: ", clusI
        #print "Outer cluster: ", clusO
        
        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
        #FWHM_it = 0.8 # ns
        FWHM_ip = 59 # mm
        FWHM_ip = 30 # mm
        print "FWHM inner z-pos: ", FWHM_ip
        #print "FWHM inner mtd: ", FWHM_it
        cl_i_merge  = []   
        added_bars = []
        for cluI in clusI:

            cl_i_zs = innerzPos[np.nonzero(np.in1d(np.arange(32),np.array(cluI)))]

            if len(cluI) == 1:
                cl_i_merge.append([0, [int(cluI[0])]])
            else:     
                all_combs = []          
                for n,(i, ii) in enumerate(zip(cl_i_zs, cluI)):
                    for m, (j, jj) in enumerate(zip(cl_i_zs, cluI)):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = ii # bar number 1
                            bar2 = jj # bar number 2
                            
                            zposdiff = i - j

                            meantimediff = 999 #meantime1 - meantime2

                            all_combs.append([0,[bar1,bar2], zposdiff, meantimediff])


                
                    
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_ip and abs(ac[3]) < 1000:
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        
                        if not bar1 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar1] in all_pot_cls:
                                pot_new_cl.append(bar1)
                            elif not bar1 in pot_new_cl and bar1-1 in pot_new_cl:
                                pot_new_cl.append(bar1)
                                
                        if not bar2 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar2] in all_pot_cls:
                                pot_new_cl.append(bar2)
                            elif not bar2 in pot_new_cl and bar2-1 in pot_new_cl:
                                pot_new_cl.append(bar2)
                    else:
                        if len(pot_new_cl) != 0 and not pot_new_cl in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)
                        
                if len(pot_new_cl) != 0 and pot_new_cl not in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)  
                                          

                for i in all_pot_cls:
                    ccc = [0,i]
                    if not ccc in cl_i_merge:
                        cl_i_merge.append(ccc)



        print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"  
        
                  
        FWHM_op = 73 # mm   
        FWHM_op = 50 # mm        

        #print "FWHM outer z-pos: ", FWHM_op    
        #print "FWHM outer mtd: ", FWHM_ot          
        cl_o_merge = []
        added_bars = []
        
        for cluO in clusO:

            cl_o_zs = outerzPos[np.nonzero(np.in1d(np.arange(32),np.array(cluO)))]

            if len(cluO) == 1:
                cl_o_merge.append([1, [int(cluO[0])]])
            else:     
                all_combs = []          
                for n,(i, ii) in enumerate(zip(cl_o_zs, cluO)):
                    for m, (j, jj) in enumerate(zip(cl_o_zs, cluO)):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = ii # bar number 1
                            bar2 = jj # bar number 2
                            
                            zposdiff = i - j

                            meantimediff = 999 #meantime1 - meantime2

                            all_combs.append([1,[bar1,bar2], zposdiff, meantimediff])

                
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    #print "ac", ac
                    if abs(ac[2]) < FWHM_op and abs(ac[3]) < 1000:
                        #print "ac3", ac[3]
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        #print "yas", ac
                        #print pot_new_cl
                        if not bar1 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar1] in all_pot_cls:
                                pot_new_cl.append(bar1)
                                #print "here 1!"
                            elif not bar1 in pot_new_cl and bar1-1 in pot_new_cl:
                                pot_new_cl.append(bar1)
                                #print "here 2!"
                                
                        if not bar2 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar2] in all_pot_cls:
                                pot_new_cl.append(bar2)
                                #print "here 3!"
                            elif not bar2 in pot_new_cl and bar2-1 in pot_new_cl:
                                pot_new_cl.append(bar2)
                                #print "here 4!"
                    else:
                        if len(pot_new_cl) != 0 and not pot_new_cl in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)
                            #print "here 5!"
                        
                if len(pot_new_cl) != 0 and pot_new_cl not in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)  
                            #print "here 6!"
                                          
                #print "pot all : "
                for i in all_pot_cls:
                    ccc = [1,i]
                    if not ccc in cl_o_merge:
                        cl_o_merge.append(ccc)
                    #print i
                                                                                                          

        # if merge_cl_i still contains neighbouring bars, check again: ... while? recursive?
        # (cluster with more than two bars hit)


        # if a cluster is recognized as two hits, they are not yet added.... so add them now:        
        cl_i_list = [item for sublist in [c[1] for c in cl_i_merge] for item in sublist] # list of all bars in clusters (flatten cluster list)
        all_i_list = [item for sublist in  clusI for item in sublist] # list of ALL bars
        
        
        cl_o_list = [item for sublist in [c[1] for c in cl_o_merge] for item in sublist]
        all_o_list = [item for sublist in clusO for item in sublist] # list of all bars (flatten cluster list)

        for b in all_i_list:
            if b not in cl_i_list:
                cl_i_merge.append([0,[b]])
                
                
        for b in all_o_list:
            if b not in cl_o_list:
                cl_o_merge.append([1,[b]])
        
        #print "cl_i_merge ", cl_i_merge      
        #print "cl_o_merge ", cl_o_merge        
                             

        #####################################################################################################################"""
        ##################################################################################################################### 

        union_cl_I = merge_rects_of_cluster(cmap[:nI], cl_i_merge) # union contains the rects of the bars
        union_cl_O = merge_rects_of_cluster(cmap[nI:], cl_o_merge)  
        
        #print "uuuuuuuuuuuuuuunion", union_cl_I, union_cl_O
        
        # MODIFY PMAP
        pmapI_new = []
        #print "pmapI ", pmapI[:,2]
        
        for cl in cl_i_merge:
            cl_pmap = []
            for cb in cl[1]:
                cl_pmap.extend(pmapI[pmapI[:,2]==cb])
            #print cl_pmap
            pmapI_new.append(cl_pmap)
                    
        #print  pmapI_new
        #####################################################################################################################                 
        #exit(0)
        #union_cl_O = get_rects_of_hits(cmapO)
        #union_cl_I = get_rects_of_hits(cmapI)
       
        hitcollection = []   # 2dim. lines : tracklines columns: hit bars
        hitcollection_count = []
        p0 = [0,0]
    #    t1 = []
    #    t2 = []
        p1 = [0,0]
        #pbgo = [BGOXPos, BGOYPos]
        line_points_collection = []
        bgo_points_coll = []
        len_bgo_list = 20 # len(polygon_hitlist)
        
        start_time = time.time()
        
        ##################################################################################################################### 
        ##################################################################################################################### 
        #final_tracks = []
        #line_points_collection = []
        #bgo_points_coll = []
        #return final_tracks, line_points_collection, bgo_points_coll
        
        ##################################################################################################################### 
        ##################################################################################################################### 
         
        hO = cl_o_merge
        hI = cl_i_merge
        range_accept1 = range(0,3)
        range_accept2 = range(31,29)
        
        #print len(cl_i_merge)
        #print len(pmapI_new)
               
        #fig2, ax2 = plt.subplots(frameon=False, figsize = (10,10))
        
        for m_cl, (curr_cl, pmap_cl) in enumerate(zip(cl_i_merge, pmapI_new)):
            #print "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            #print "curr cluster", curr_cl
            for m,curr_pos in enumerate(pmap_cl): # loop over inner bars i[0],i[1] are x,y coord, i[2] is the bar number
                    p0[0] = curr_pos[0]
                    p0[1] = curr_pos[1]         
                    for i in range(0,int(len_bgo_list*1.1)):  # draw XX points from the circle / polygons aka all pmt pixel with a hit. make it dependent of the size of the hit shape
                    
                        p1 = draw_pt_from_circle(circle.radius,[BGOXPos,BGOYPos])
                        #p1 = draw_point_from_bgohit(polygon_hitlist, bgo_hit_bounds)
                        p1.append(i)
                        
                        if plot_trackfinding[0] == True:
                            bgo_points_coll.append(p1)
                        
                        hitbars = []  # contains bars on one line
                        hitbars.append(curr_cl)
       
                        for count,(pol,nrcO) in enumerate(zip(union_cl_O,hO)):  # for ever line, check if it intersects one of the hit bars (outer)  
                            nrO = nrcO[1]  # bar or cluster         
                            if not round(abs(sum(curr_cl[1])/len(curr_cl[1]) - sum(nrO)/len(nrO))) in range_accept1 and not round(abs(sum(curr_cl[1])/len(curr_cl[1]) - sum(nrO)/len(nrO))) in range_accept2:
                            # dont look on the other side of hodor
                                continue   
                             
                            yes = Polygon_contains_line(pol,p0,p1[:-1])
                            
                            if yes == True:
                                #x,y = pol.exterior.xy       
                                #ax2.plot(x, y, color='#6699cc', alpha=0.7,
                                #linewidth=3, solid_capstyle='round', zorder=2)   
   
                                if not [1,nrO] in hitbars and [1,[nrO[0]]] not in hitbars:                                
                                    hitbars.append([1,nrO]) 
                                    
                                    if plot_trackfinding[1]== True:    
                                        line_points_collection.append([p0[0],p0[1],p1[0],p1[1]]) #this is only used for drawing                           
                                     
                        sorted_track = sorttrack(hitbars)
                        
                        if len(hitcollection) == 0 and len(sorted_track) == 2: # only  tracks of length 2
                            hitcollection.append(sorted_track)
                            hitcollection_count.append(1)

                        if not sorted_track in hitcollection and len(sorted_track) == 2: # only tracks of length 2
                           hitcollection.append(sorted_track)
                           hitcollection_count.append(1)
                           
                        elif sorted_track in hitcollection:
                           ind = hitcollection.index(sorted_track) 
                           hitcollection_count[ind] += 1
                                
                                
                                
        all_lines = sum(hitcollection_count)
        
        #plt.show()
 
        
        #plot_all_bars_2D(cmapI,cmapO, union_cl_I, union_cl_O, line_points_collection)
        
        #merge_coll = []
        for h,c in zip(hitcollection, hitcollection_count):
            h.append(float(c)/float(all_lines))

        hitcollection.sort(key=len,reverse=True)
        
        hitcollection2 = []
        #print "################################## END OF ROTATION ALG, time (s): ", (time.time() - start_time)
        for h in hitcollection:
            
            if len(h[:-1])>=2: # and len(): 
            
                hitcollection2.append(h)
                #print h
        #print "##################################"    
         
        if len(hitcollection2) == 0:
            return 999, 999, 999, union_cl_I, union_cl_O   
            
        final_tracks = new_select_tracks_for_2D(hitcollection2, currI, currO) 
        
        if final_tracks == 999:
            return 999, 999, 999,union_cl_I, union_cl_O
        
        return final_tracks, line_points_collection, bgo_points_coll, union_cl_I, union_cl_O#"""
    
##############################################################################################################   
############################################################################################################################################################################################################################  
############################################################################################################################################################################################################################  


    
def track_finding_polygons_tpx(pmapI, cmap, currI, currO, polygon_hitlist, plot_trackfinding, bgo_hit_bounds, innerzPos, outerzPos): #, BGOXPos, BGOYPos):  # with density of lines added

        ### First: Find Clusters ###
        #print "####################"
        #print "inner bars hit: ", cmap[:nI,1] # corner map, all rects with corner points plus nr
        #print "outer bars hit: ", cmap[nI:,1]     
        #plt.show()        
        #print "pmapI ", pmapI
        
        
        
        #print "Outer CFs "
        #print outerCFs

        nI = len(currI)
        nO = len(currO)
       
        cmapI = cmap[:nI,:]
        cmapO = cmap[nI:,:] 
        
        # use this for no unification of rects
        yesnoI, cI = find_consecutive_nr(currI)  
        yesnoO, cO = find_consecutive_nr(currO)
        
        merge_cl_i = []
        merge_cl_o = []
        #####################################################################################################################
        #####################################################################################################################         
        ynI, clusI = find_consecutive_nr_singletoo(currI)
        ynO, clusO = find_consecutive_nr_singletoo(currO)
               
        #print "Inner cluster: ", clusI
        #print "Outer cluster: ", clusO
        
        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
        #FWHM_it = 0.8 # ns
        FWHM_ip = 59 # mm
        FWHM_ip = 30 # mm
        print "FWHM inner z-pos: ", FWHM_ip
        #print "FWHM inner mtd: ", FWHM_it
        cl_i_merge  = []   
        added_bars = []
        for cluI in clusI:

            cl_i_zs = innerzPos[np.nonzero(np.in1d(np.arange(32),np.array(cluI)))]

            if len(cluI) == 1:
                cl_i_merge.append([0, [int(cluI[0])]])
            else:     
                all_combs = []          
                for n,(i, ii) in enumerate(zip(cl_i_zs, cluI)):
                    for m, (j, jj) in enumerate(zip(cl_i_zs, cluI)):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = ii # bar number 1
                            bar2 = jj # bar number 2
                            
                            zposdiff = i - j

                            meantimediff = 999 #meantime1 - meantime2

                            all_combs.append([0,[bar1,bar2], zposdiff, meantimediff])


                
                    
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_ip and abs(ac[3]) < 1000:
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        
                        if not bar1 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar1] in all_pot_cls:
                                pot_new_cl.append(bar1)
                            elif not bar1 in pot_new_cl and bar1-1 in pot_new_cl:
                                pot_new_cl.append(bar1)
                                
                        if not bar2 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar2] in all_pot_cls:
                                pot_new_cl.append(bar2)
                            elif not bar2 in pot_new_cl and bar2-1 in pot_new_cl:
                                pot_new_cl.append(bar2)
                    else:
                        if len(pot_new_cl) != 0 and not pot_new_cl in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)
                        
                if len(pot_new_cl) != 0 and pot_new_cl not in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)  
                                          

                for i in all_pot_cls:
                    ccc = [0,i]
                    if not ccc in cl_i_merge:
                        cl_i_merge.append(ccc)



        print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"  
        
                  
        FWHM_op = 73 # mm   
        FWHM_op = 50 # mm        

        #print "FWHM outer z-pos: ", FWHM_op    
        #print "FWHM outer mtd: ", FWHM_ot          
        cl_o_merge = []
        added_bars = []
        
        for cluO in clusO:

            cl_o_zs = outerzPos[np.nonzero(np.in1d(np.arange(32),np.array(cluO)))]

            if len(cluO) == 1:
                cl_o_merge.append([1, [int(cluO[0])]])
            else:     
                all_combs = []          
                for n,(i, ii) in enumerate(zip(cl_o_zs, cluO)):
                    for m, (j, jj) in enumerate(zip(cl_o_zs, cluO)):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = ii # bar number 1
                            bar2 = jj # bar number 2
                            
                            zposdiff = i - j

                            meantimediff = 999 #meantime1 - meantime2

                            all_combs.append([1,[bar1,bar2], zposdiff, meantimediff])

                
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    #print "ac", ac
                    if abs(ac[2]) < FWHM_op and abs(ac[3]) < 1000:
                        #print "ac3", ac[3]
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        #print "yas", ac
                        #print pot_new_cl
                        if not bar1 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar1] in all_pot_cls:
                                pot_new_cl.append(bar1)
                                #print "here 1!"
                            elif not bar1 in pot_new_cl and bar1-1 in pot_new_cl:
                                pot_new_cl.append(bar1)
                                #print "here 2!"
                                
                        if not bar2 in pot_new_cl:
                            if len(pot_new_cl) == 0 and not [bar2] in all_pot_cls:
                                pot_new_cl.append(bar2)
                                #print "here 3!"
                            elif not bar2 in pot_new_cl and bar2-1 in pot_new_cl:
                                pot_new_cl.append(bar2)
                                #print "here 4!"
                    else:
                        if len(pot_new_cl) != 0 and not pot_new_cl in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)
                            #print "here 5!"
                        
                if len(pot_new_cl) != 0 and pot_new_cl not in all_pot_cls:
                            all_pot_cls.append(pot_new_cl)  
                            #print "here 6!"
                                          
                #print "pot all : "
                for i in all_pot_cls:
                    ccc = [1,i]
                    if not ccc in cl_o_merge:
                        cl_o_merge.append(ccc)
                    #print i
                                                                                                          

        # if merge_cl_i still contains neighbouring bars, check again: ... while? recursive?
        # (cluster with more than two bars hit)


        # if a cluster is recognized as two hits, they are not yet added.... so add them now:        
        cl_i_list = [item for sublist in [c[1] for c in cl_i_merge] for item in sublist] # list of all bars in clusters (flatten cluster list)
        all_i_list = [item for sublist in  clusI for item in sublist] # list of ALL bars
        
        
        cl_o_list = [item for sublist in [c[1] for c in cl_o_merge] for item in sublist]
        all_o_list = [item for sublist in clusO for item in sublist] # list of all bars (flatten cluster list)
        
        print "ALL inner: ", cl_i_list
        print "ALL outer: ", cl_o_list
        
        for b in all_i_list:
            if b not in cl_i_list:
                cl_i_merge.append([0,[b]])
                
                
        for b in all_o_list:
            if b not in cl_o_list:
                cl_o_merge.append([1,[b]])
        
        print "cl_i_merge ", cl_i_merge      
        print "cl_o_merge ", cl_o_merge                             
        print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
        #####################################################################################################################"""
        ##################################################################################################################### 

        union_cl_I = merge_rects_of_cluster(cmap[:nI], cl_i_merge) # union contains the rects of the bars
        union_cl_O = merge_rects_of_cluster(cmap[nI:], cl_o_merge)  
        
        # MODIFY PMAP
        pmapI_new = []
        #print "pmapI ", pmapI[:,2]
        
        for cl in cl_i_merge:
            cl_pmap = []
            for cb in cl[1]:
                cl_pmap.extend(pmapI[pmapI[:,2]==cb])
            #print cl_pmap
            pmapI_new.append(cl_pmap)
                    
        #print  pmapI_new
        #####################################################################################################################                 
        #exit(0)
        #union_cl_O = get_rects_of_hits(cmapO)
        #union_cl_I = get_rects_of_hits(cmapI)
       
        hitcollection = []   # 2dim. lines : tracklines columns: hit bars
        hitcollection_count = []
        p0 = [0,0]
    #    t1 = []
    #    t2 = []
        p1 = [0,0]
        #pbgo = [BGOXPos, BGOYPos]
        line_points_collection = []
        bgo_points_coll = []
        len_bgo_list = 50
        
        start_time = time.time()
        
        ##################################################################################################################### 
        ##################################################################################################################### 
        #final_tracks = []
        #line_points_collection = []
        #bgo_points_coll = []
        #return final_tracks, line_points_collection, bgo_points_coll
        
        ##################################################################################################################### 
        ##################################################################################################################### 
         
        hO = cl_o_merge
        hI = cl_i_merge
        range_accept1 = range(0,3)
        range_accept2 = range(31,29)
        
        #print len(cl_i_merge)
        #print len(pmapI_new)
               
        #fig2, ax2 = plt.subplots(frameon=False, figsize = (10,10))
        
        for m_cl, (curr_cl, pmap_cl) in enumerate(zip(cl_i_merge, pmapI_new)):
            #print "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            #print "curr cluster", curr_cl
            for m,curr_pos in enumerate(pmap_cl): # loop over inner bars i[0],i[1] are x,y coord, i[2] is the bar number
                    p0[0] = curr_pos[0]
                    p0[1] = curr_pos[1]         
                    for i in range(0,int(len_bgo_list*1.1)):  # draw XX points from the circle / polygons aka all pmt pixel with a hit. make it dependent of the size of the hit shape
                    
                        #p1 = draw_pt_from_circle(circle.radius,[BGOXPos,BGOYPos])
                        p1 = draw_point_from_tpx(polygon_hitlist, bgo_hit_bounds)
                        p1.append(i)
                        
                        if plot_trackfinding[0] == True:
                            bgo_points_coll.append(p1)
                        
                        hitbars = []  # contains bars on one line
                        hitbars.append(curr_cl)
       
                        for count,(pol,nrcO) in enumerate(zip(union_cl_O,hO)):  # for ever line, check if it intersects one of the hit bars (outer)  
                            nrO = nrcO[1]  # bar or cluster         
                            if not round(abs(sum(curr_cl[1])/len(curr_cl[1]) - sum(nrO)/len(nrO))) in range_accept1 and not round(abs(sum(curr_cl[1])/len(curr_cl[1]) - sum(nrO)/len(nrO))) in range_accept2:
                            # dont look on the other side of hodor
                                continue   
                             
                            yes = Polygon_contains_line(pol,p0,p1[:-1])
                            
                            if yes == True:
                                #x,y = pol.exterior.xy       
                                #ax2.plot(x, y, color='#6699cc', alpha=0.7,
                                #linewidth=3, solid_capstyle='round', zorder=2)   
   
                                if not [1,nrO] in hitbars and [1,[nrO[0]]] not in hitbars:                              
                                    hitbars.append([1,nrO]) 
                                    
                                    if plot_trackfinding[1]== True:    
                                        line_points_collection.append([p0[0],p0[1],p1[0],p1[1]]) #this is only used for drawing                           
                                     
                        sorted_track = sorttrack(hitbars)
                        
                        if len(hitcollection) == 0 and len(sorted_track) == 2: # only  tracks of length 2
                            hitcollection.append(sorted_track)
                            hitcollection_count.append(1)

                        if not sorted_track in hitcollection and len(sorted_track) == 2: # only tracks of length 2
                           hitcollection.append(sorted_track)
                           hitcollection_count.append(1)
                           
                        elif sorted_track in hitcollection:
                           ind = hitcollection.index(sorted_track) 
                           hitcollection_count[ind] += 1
                                
                                
                                
        all_lines = sum(hitcollection_count)
        
        #plt.show()
        
        
        #merge_coll = []
        for h,c in zip(hitcollection, hitcollection_count):
            h.append(float(c)/float(all_lines))

        hitcollection.sort(key=len,reverse=True)
        
        hitcollection2 = []
        #print "################################## END OF ROTATION ALG, time (s): ", (time.time() - start_time)
        for h in hitcollection:
            
            if len(h[:-1])>=2: # and len(): 
            
                hitcollection2.append(h)
                print h
        #print "##################################"    
         
        if len(hitcollection2) == 0:
            return 999, 999, 999, union_cl_I, union_cl_O
            
        final_tracks = new_select_tracks_for_2D(hitcollection2, currI, currO) 
        
        if final_tracks == 999:
            return 999, 999, 999, union_cl_I, union_cl_O
        
        return final_tracks, line_points_collection, bgo_points_coll, union_cl_I, union_cl_O
############################################################################################################################################################################################################################ 
############################################################################################################################################################################################################################      

def new_select_tracks_for_2D(hitcollection, currI, currO): # has the density of lines added!
    print "select tracks ############################################################################"

    #exit(0)

    # sort according to probability        
    hitcollection = list(reversed(sorted(hitcollection, key=itemgetter(-1))))
    
    """print "------------------------------------------------------------------------------------------"  
                         
    for i in hitcollection:
        print i                           

    print "------------------------------------------------------------------------------------------" # """
    
    
    finalhitcoll = []
    # add the first track with the highest prob
    finalhitcoll.append(hitcollection[0])
    
    #print "currI ", currI
    #print "currO ", currO
        
    #print hitcollection[0][0][1]
    #print hitcollection[0][1][1]
    for i in hitcollection[0][0][1]:
        #print "i", i
        currI = currI[currI != i]
        #print "currI", currI
        
    for i in hitcollection[0][1][1]:
    
        currO = currO[currO != i]    
        
    #print "currI after ", currI
    #print "currO after ", currO    
       
    for tr in hitcollection: # loop over all possible found tracks
            #chosen = False
            track = tr[:-1] #not the line count
            #inner = 0
            #outer = 1
            if tr not in finalhitcoll:
                #print "tr ", tr
                #print "tr if ", tr[0][1]
                """for trh in tr[0][1]: # in case there is a cluster
                    print "trh i", trh
                    if trh in currI and len(currI) != 0: 
                        finalhitcoll.append(tr) # but here we append with line count 
                        currI = currI[currI != trh]
                    elif len(currI) == 0: 
                        break
                        
                        
                for trh in tr[1][1]: # in case there is a cluster
                    print "trh o", trh
                    if trh in currO and len(currO) != 0: 
                        finalhitcoll.append(tr) # but here we append with line count
                        currO = currO[currO != trh]
                    elif len(currO) == 0:
                        break"""
                
                for trh in tr[0][1]: # in case there is a cluster
                    #print "trh i", trh
                    for trho in tr[1][1]:
                        if (trh in currI and len(currI)!=0) and (trho in currO and len(currO) != 0): 
                            #print "in 1 "
                            finalhitcoll.append(tr) # but here we append with line count"""
                            for i in tr[0][1]:
                                currI = currI[currI != i]
                            for o in tr[1][1]:
                                currO = currO[currO != o]
                        else:
                            break
                            


    finalhitcoll = [s for s in finalhitcoll if len(s[0]) in [2] ]
    finalhitcoll.sort()
    finalhitcoll = list(finalhitcoll for finalhitcoll,_ in itertools.groupby(finalhitcoll))
    
    #print "print after probability sort: -----------------------------------------------------------"  
                         
    #for i in finalhitcoll:
    #    print i                           

    #print "------------------------------------------------------------------------------------------" # """
    #finalhitcoll_set = list(set(finalhitcoll))
    #finalhitcoll = list(set(finalhitcoll) - set(finalhitcoll_set))

    return finalhitcoll 
    #return ultimatehitcoll

##############################################################################################################
##############################################################################################################
def draw_pt_from_circle(radius, midpoint): # get a point inside circle - uniform distributed
    point = [0,0]
    t = 2*math.pi*np.random.uniform(0,1,1)
    u = np.random.uniform(0,1,1)+np.random.uniform(0,1,1)
    if u>1:
        r = 2-u
    else:
        r = u
    point[0] = midpoint[0] + r*radius*math.cos(t)
    point[1] = midpoint[1] + r*radius*math.sin(t)
    
    #ellipse:
    # phi is angle between X ais and majorj axis of ellipse
    #point[0] = midpoint[0] + a*math.cos(t)*math.cos(phi) - b*math.sin(t)*math.sin(phi)
    #point[1] = midpoint[1] + a*math.cos(t)*math.sin(phi) + b*math.sin(t)*math.cos(phi)
    
    return point
    
    
    
##############################################################################################################
##############################################################################################################  
  
def draw_point_from_tpx(rects, bgo_hit_bounds): 
    # rects is a list of rects... merging them does not work, if they are seperated..
    #bgo_bound = 50 # 5 cm
    picked = False
       
    r1  = random.uniform(bgo_hit_bounds[0][0], bgo_hit_bounds[0][1]) 
    r2 =  random.uniform(bgo_hit_bounds[1][0], bgo_hit_bounds[1][1])            
          
          
    return [r1,r2]
    
##############################################################################################################   
##############################################################################################################

def draw_point_from_bgohit(rects, bgo_hit_bounds): 
    # rects is a list of rects... merging them does not work, if they are seperated..
    #bgo_bound = 50 # 5 cm
    picked = False
       
    while picked == False:
    
        #r1 = np.random.uniform(-50,50,1)
        #r2 = np.random.uniform(-50,50,1)
        r1  = random.uniform(bgo_hit_bounds[0][0], bgo_hit_bounds[0][1]) 
        r2 =  random.uniform(bgo_hit_bounds[1][0], bgo_hit_bounds[1][1])            
        for r in rects:
            #print r
            if r.contains(Point(r1,r2)) == True:
                #print  "jaaaaaa ", r1,r2
                randompoint = [r1,r2]
                picked = True
                
                break
               
    return randompoint
    
##############################################################################################################   
##############################################################################################################
def Polygon_contains_line(polygon,p0,p1): # funktioniert aber mit allem auch ellipse oder so # ax uebergeben
    #get max and min values of polygon
    bounds = polygon.bounds
    x = np.arange(bounds[0],bounds[2],0.7)  # xmin and xmax
    #print "xxxxxxxxxxxxxxxxxxxxxxxxx"
    #print "x len", len(x)
    #print x
    #print "xxxxxxxxxxxxxxxxxxxxxxxxx"
    y = np.zeros(len(x))    
    y = (p1[1] - p0[1])/(p1[0] - p0[0])*(x - p0[0]) + p0[1]
    yes = np.logical_and(y[:]>bounds[1], y[:]<bounds[3])
    y = y[yes]
    x = x[yes]
    
    if polygon.geom_type == 'Polygon':
        for n,(i,j) in enumerate(zip(x,y)):
            p = Point(i,j)
            if polygon.contains(p) == True: 
                return True
                
    elif polygon.geom_type == 'MultiPolygon':
        for pol in polygon:    
            for n,(i,j) in enumerate(zip(x,y)):
                p = Point(i,j) 
                if pol.contains(p) == True: 
                    return True

    return False
##############################################################################################################
def merge_rects_of_cluster(cmap, cluster): # cmap for inner or outer   
    #print "Corner Map:", cmap # Corner map of inner or outer respectively
    # First: get rects of cluster   
    
    polygon_cmap = []
    for cl in cluster:
        temp_pol = []
        #print cl
        for b in cl[1]:
            temp_pol.append(cmap[cmap[:,1]==b])           
        polygon_cmap.append(temp_pol)
        
    new_polygon_cmap = [[i[0] for i in l] for l in polygon_cmap]
    
    #print "new "
    #print [ [i[1] for i in l] for l in new_polygon_cmap]
   
    union_pol_list = []    
    for cllust in new_polygon_cmap:
        polygon_list = []
        for i in range(0,len(cllust)):
            x0 = cllust[i][0][2][0]
            y0 = cllust[i][0][2][1]
            x1 = cllust[i][0][0][0]
            y1 = cllust[i][0][0][1]
            x2 = cllust[i][0][1][0]
            y2 = cllust[i][0][1][1]
            x3 = cllust[i][0][3][0]
            y3 = cllust[i][0][3][1]
            ext = [(x0, y0), (x1, y1), (x2, y2), (x3,y3)]
            polygon = Polygon(ext)
            polygon_list.append(polygon)
    
        up = cascaded_union(polygon_list)
        union_pol_list.append(up)
           
    return union_pol_list
##############################################################################################################
def get_rects_of_hits(cmap): # no unification of cluster hits to bigger rects
    polygon_list = []
    #print cmap
    for i in cmap:
            x0 = i[0][2][0]
            y0 = i[0][2][1]
            x1 = i[0][0][0]
            y1 = i[0][0][1]
            x2 = i[0][1][0]
            y2 = i[0][1][1]
            x3 = i[0][3][0]
            y3 = i[0][3][1]
            ext = [(x0, y0), (x1, y1), (x2, y2), (x3,y3)]
            polygon = Polygon(ext)
            polygon_list.append(polygon)
            
    return polygon_list
##############################################################################################################
##############################################################################################################
def sorttrack(hb):
    sorted_by_second = sorted(hb, key=lambda tup: tup[1])
    return sorted(sorted_by_second, key=lambda tup: tup[0])
##############################################################################################################

