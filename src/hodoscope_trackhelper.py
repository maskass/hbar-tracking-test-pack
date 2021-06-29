import numpy as np
import math 

##############################################################################################################    
def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]
##############################################################################################################      
def calc_angle_for_1track(trackcoll, vertex):    
    
        positionmap = getPositionMap()
        #print "Trackcoll before: "
        #for temp in trackcoll:
        #    print temp

        pos_inner = positionmap[0]
        pos_outer = positionmap[1]

        add = np.arange(32)

        pos_inner = np.column_stack((add, pos_inner))
        pos_outer = np.column_stack((add, pos_outer))
        #trackcoll_temp = list(trackcoll)
        tr = trackcoll[0]
        anglesY = []
        upwards = 0
        downwards = 0
        horizont = 0

       
        #print "n, m ", n, m
        merge_o = tr[1][1]
        
        ro_merge = pos_outer[merge_o]
        ro_new = np.mean(ro_merge[:,0:], axis = 0)
        
        vA = [(vertex[0]-ro_new[1]), (vertex[1]-ro_new[2])]
        magA = dot(vA, vA)**0.5   
                 
        vY = [0,1]
        vmY = [0,-1]
                    
        # Get dot prod
        dot_prodY = dot(vA/magA, vY)
        #print "dot_prodY ", dot_prodY
        if abs(dot_prodY) < 0.30:
            horizont += 1
            #print "horiz "
        elif dot_prodY >= 0.30:
            downwards += 1
            #print "down "
        elif dot_prodY < 0:
            upwards += 1
            #print "up "
                   
        dot_prodmY = dot(vA/magA, vmY)
        # Get magnitudes
        #magY = 1.0
        # Get cosine value
        cos_Y = dot_prodY
        cos_mY = dot_prodmY
                    
        angleY = math.acos(cos_Y)
        angleY = math.degrees(angleY)%360
        anglemY = math.acos(cos_mY)
        anglemY = math.degrees(anglemY)%360
        anglesY.append(min(angleY, anglemY))
                    
        #print "ANGLES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        #print anglesY
        #print "ANGLES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        anglesY = np.array(anglesY)
        
        return anglesY, [upwards, downwards, horizont]

###################################################################################################################################################################################

def calc_mtds_of_tracks_corr_std(trackcoll, cfs_i, cfs_o, zinner, zouter, vertex): #, rinner, router):
    
    mtds = []
    #print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    print cfs_i
    positionmap = getPositionMap()
    #print "Trackcoll before: "
    #for temp in trackcoll:
    #    print temp

    pos_inner = positionmap[0]
    pos_outer = positionmap[1]

    add = np.arange(32)

    pos_inner = np.column_stack((add, pos_inner))
    pos_outer = np.column_stack((add, pos_outer))

    angles = [999,999,999]
    #print "#################################"
    if len(trackcoll) >= 2:
    
        trackcoll_temp = list(trackcoll)
        angles = []
        anglesY = []
        upwards = 0
        downwards = 0
        horizont = 0
        for n, tr in enumerate(trackcoll_temp):            
            for m, tr2 in enumerate(trackcoll_temp):                              
                #if tr != tr2 and n < m: a
                if m == n+1 or (n == len(trackcoll)-1 and m == 0):                   
                    #print "n, m ", n, m
                    #print "tr ", tr
                    #print "tr2 ", tr2
                    merge_o = tr[1][1]
                    merge_o2 = tr2[1][1]
                    ro_merge = pos_outer[merge_o]
                    ro_merge2 = pos_outer[merge_o2]

                    ro_new = np.mean(ro_merge[:,0:], axis = 0)
                    ro_new2 = np.mean(ro_merge2[:,0:], axis = 0) 

                    vA = [(vertex[0]-ro_new[1]), (vertex[1]-ro_new[2])]
                    vB = [(vertex[0]-ro_new2[1]), (vertex[1]-ro_new2[2])]                    
                    vY = [0,1]
                    vmY = [0,-1]
                    
                    # Get dot prod
                    dot_prod = dot(vA, vB)
                    # Get magnitudes
                    magA = dot(vA, vA)**0.5
                    magB = dot(vB, vB)**0.5
                    # Get cosine value
                    cos_ = dot_prod/magA/magB
                    
                    # Get dot prod
                    dot_prodY = dot(vA/magA, vY)
                    #print "dot_prodY ", dot_prodY
                    if abs(dot_prodY) < 0.30:
                        horizont += 1
                        #print "horiz "
                    elif dot_prodY >= 0.30:
                        downwards += 1
                        #print "down "
                    elif dot_prodY < 0:
                        upwards += 1
                        #print "up "
                    
                    dot_prodmY = dot(vA/magA, vmY)
                    # Get magnitudes
                    #magY = 1.0
                    # Get cosine value
                    cos_Y = dot_prodY
                    cos_mY = dot_prodmY
                    
                    angleY = math.acos(cos_Y)
                    angleY = math.degrees(angleY)%360
                    anglemY = math.acos(cos_mY)
                    anglemY = math.degrees(anglemY)%360
                    #print "Y ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    #print angleY, anglemY
                    #print "Y ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    anglesY.append(min(angleY, anglemY))
                    
                    if abs(cos_) > 1:
                        angles.append(999)
                    else:                    
                        # Get angle in radians and then convert to degrees
                        angle = math.acos(cos_)
                        # Basically doing angle <- angle mod 360
                        angle = math.degrees(angle)%360
                        angles.append(angle)
                    #print "angle ", angle 
                    
                    
        #print "ANGLES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        #print angles 
        #print anglesY
        #print "ANGLES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        angles = np.array(angles)
        #print len(angles[angles < 30])
        #trackcoll_temp = list(trackcoll)
        for n,j in enumerate(angles):            
       
            if j < 30 and len(trackcoll) > 2:
                #print "NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW!!! ANGLE", j
                
                n1 = n
                n2 = n + 1
                if n == len(angles)-1:
                    n2 = 0
                    
                track1 = trackcoll[n1]
                track2 = trackcoll[n2]
                                
                """if n == 0:
                    track1 = trackcoll[0]
                    track2 = trackcoll[1]
                elif n == 1:
                    track1 = trackcoll[1]               # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! attention to the order!!!
                    track2 = trackcoll[2]
                elif n == 2:
                    track1 = trackcoll[2]
                    track2 = trackcoll[0]"""

                #print "track1 ", track1
                #print "track2 ", track2
                
                
                merge_i = track1[0][1]
                merge_i.append(track2[0][1][0])
                merge_o = track1[1][1]
                merge_o.append(track2[1][1][0])
                merge_i = list(set(merge_i))
                merge_o = list(set(merge_o))
                ##########################################################################
                merge_cfi = cfs_i[np.nonzero(np.in1d(cfs_i[:,0],np.array(merge_i)))[0]] 
                cfi_new = np.mean(merge_cfi[:,0:], axis = 0) 
                cfi_new[0] = cfi_new[0] + 999 

                cfs_i  = np.vstack([cfs_i, cfi_new])
                ##########################################################################
                merge_cfo = cfs_o[np.nonzero(np.in1d(cfs_o[:,0],np.array(merge_o)))[0]] 
                cfo_new = np.mean(merge_cfo[:,0:], axis = 0) 
                cfo_new[0] = cfo_new[0] + 999  
                cfs_o  = np.vstack([cfs_o, cfo_new])


                ##########################################################################
                merge_zi = zinner[np.nonzero(np.in1d(zinner[:,0],np.array(merge_i)))[0]] #[0][1]     
                merge_zo = zouter[np.nonzero(np.in1d(zouter[:,0],np.array(merge_o)))[0]] #[0][1]

                zi_track = np.mean(merge_zi[:,0:], axis = 0)
                zo_track = np.mean(merge_zo[:,0:], axis = 0)
                ##########################################################################
                zinner_new = np.mean(merge_zi[:,0:], axis = 0)  
                zinner_new[0] = zinner_new[0] + 999 
                zinner  = np.vstack([zinner, zinner_new])

                zouter_new = np.mean(merge_zo[:,0:], axis = 0) 
                zouter_new[0] = zouter_new[0] + 999 
                zouter  = np.vstack([zouter, zouter_new])
                ##########################################################################
                ri_merge = pos_inner[merge_i]
                ro_merge = pos_outer[merge_o]

                ri_new = np.mean(ri_merge[:,0:], axis = 0)
                ro_new = np.mean(ro_merge[:,0:], axis = 0)
                
                ri_new[0] = ri_new[0] + 999
                ro_new[0] = ro_new[0] + 999

                pos_inner  = np.vstack([pos_inner, ri_new])
                pos_outer  = np.vstack([pos_outer, ro_new])

                new_track = [[0,[cfi_new[0] - 999]], [1,[cfo_new[0] - 999]], 999]
                #print "new track ", new_track
          
                if track1 in trackcoll_temp:
                    trackcoll_temp.remove(track1)
                if track2 in trackcoll_temp:
                    trackcoll_temp.remove(track2)

                trackcoll_temp.append(new_track)           
                #break
                
        trackcoll = trackcoll_temp
        
        #print "trackcoll_temp"
        #for temp in trackcoll_temp:
        #    print temp 
    #print "lennnnnnnnnnnnn", len(pos_outer)
    
   
    #print cfs_i

    #print cfs_o

    #print pos_inner[:,0]
    #print pos_outer[:,0]

    #print "#################################"    
    
    c = 299792458.0*1000.0/1.0e9
    
    #print "trackcoll"
    #for temp in trackcoll:
    #       print temp     
    
    track_meantimes = [] 
       
    # mean loop
    for n, track in enumerate(trackcoll): 
        #print "track ", track
                
        if track[2] == 999:
            track[0][1][0] = track[0][1][0] + 999
            track[1][1][0] = track[1][1][0] + 999


        inner_bar = track[0][1]
        outer_bar = track[1][1]
        #print inner_bar
        #print outer_bar
        
        cfi_track = cfs_i[np.nonzero(np.in1d(cfs_i[:,0],np.array(inner_bar)))[0]]
        #print "cfi", cfi_track
        zi_track = zinner[np.nonzero(np.in1d(zinner[:,0],np.array(inner_bar)))[0]] #[0][1]
        ri_track = pos_inner[np.nonzero(np.in1d(pos_inner[:,0],np.array(inner_bar)))[0]] 
    
        cfi_track = np.mean(cfi_track[:,1:], axis = 0)  # make mean of columns (in case of cluster!)
        #print "cfi track mean: ", cfi_track
        zi_track = np.mean(zi_track[:,1:], axis = 0)
        ri_track = np.mean(ri_track[:,0:], axis = 0)
                   
        ri_track =  np.sqrt((ri_track[1] - vertex[0])*(ri_track[1] - vertex[0]) + (ri_track[2] - vertex[1])*(ri_track[2] - vertex[1]) )
        track_corr = np.sqrt(ri_track*ri_track + zi_track*zi_track)/c #"""
        cfi_track = np.mean(cfi_track) - track_corr
        
        cfo_track = cfs_o[np.nonzero(np.in1d(cfs_o[:,0], np.array(outer_bar)))[0]]  
        zo_track = zouter[np.nonzero(np.in1d(zouter[:,0],np.array(outer_bar)))[0]]
        ro_track = pos_outer[np.nonzero(np.in1d(pos_outer[:,0],np.array(outer_bar)))[0]] 
                
        cfo_track = np.mean(cfo_track[:,1:], axis = 0)  # make mean of columns (in case of cluster!)
        zo_track = np.mean(zo_track[:,1:], axis = 0)
        ro_track = np.mean(ro_track[:,0:], axis = 0)
        ro_track =  np.sqrt((ro_track[1] - vertex[0])*(ro_track[1] - vertex[0]) + (ro_track[2] - vertex[1])*(ro_track[2] - vertex[1]) )

        track_ocorr = np.sqrt(ro_track*ro_track + zo_track*zo_track)/c
        cfo_track = np.mean(cfo_track) - track_ocorr
        
        io_mean = 0.5*(cfo_track + cfi_track)
        
        track_meantimes.append(io_mean)


    #print "track_meantimes ", track_meantimes
    event_mean = np.mean(np.array(track_meantimes))
    
    event_max = np.max(track_meantimes)
    event_min = np.min(track_meantimes)
    
    #print "event mean: ", event_mean
    #print "event max: ", event_max   
    #print "event min: ", event_min
        
    event_mtd = 0.0
    for tr_mt in track_meantimes:
        event_mtd = event_mtd + (tr_mt - event_mean)*(tr_mt - event_mean)
        
    #print "event mtd", event_mtd    
    event_mtd = math.sqrt(event_mtd/(1.0*len(trackcoll) - 1.0))
    event_min_max = (event_max - event_min)
    #print "event mtd ", event_mtd   
    #print "event min max ", event_min_max    
    
    #if len(angles) < 3:
    #    angles = np.append(angles,999)
    #
    return event_mtd, angles, anglesY, event_min_max, [upwards, downwards, horizont]

##################################################################################################################################################################
##################################################################################################################################################################
################################################################################################################################################################## 
# FIXME tis is from the distributions of the 2016 data   
def calc_pvals_mtds(event_mtd):    
    # FOR THE COMBINED MTD
    # cosmics: PRELIM!!!! 
    #cosmics
    #Beta: [ 0.98275134  0.78375044  0.25303222]
    #Beta Std Error: [ 0.02289982  0.00657134  0.00429869]
    mean_c = 0.78375044
    sigma_c = 0.25303222

    # pbars: PRELIM !!!! 
    #Beta: [ 0.99472814 -0.12751322  0.252425  ]
    #Beta Std Error: [ 0.02723027  0.00780577  0.00518992]

    mean_p = -0.12751322
    sigma_p = 0.252425 


    #cosmic distribution:    
    #distr_cosmic = stats.norm(loc = mean_c, scale = sigma_c)      
    # pbar distribution:    
    #distr_pbar = stats.norm(loc = mean_p, scale = sigma_p)  
        
    z_score_c = (event_mtd - mean_c) / sigma_c
    z_score_p = (event_mtd - mean_p) / sigma_p
             
    pval_c = stats.norm.sf(abs(z_score_c))
    pval_p = stats.norm.sf(abs(z_score_p)) 
 
    #pval_p = stats.norm.cdf(event_mtd + 50, loc = mean_p, scale = sigma_p)
        
    print "z-score, p-value cosmic: ", z_score_c, pval_c      
    print "z-score, p-value pbar: ", z_score_p, pval_p  
    
    return pval_c, pval_p
##############################################################################################################
def is_it_a_hit(h, cable_corr_i, cable_corr_o): 

    #print cable_corr_i
    conv_factor_i = 0.014177797004048582995951417004048582995951417004048582995 # inner to mm        
    conv_factor_o = 0.014127641127819548872180451127819548872180451127819548872 # outer to mm

    mt_cuts_inner = [105,120] # 2019
    mt_cuts_outer = [105,122] # 2019
                                    # cuts based on meantimes and zpos distributions for whole hodor with cut on bgo E of 20MeV
    z_cuts_inner = [-183,218] # 2019
    z_cuts_outer = [-300,300] # 2019
    
    t_glob_o = 1.39542395
    t_glob_i = 1.84575840
    
    # get timestamps of active bars
    new_active_bars = []
    nI, nO, I, O = h.getActiveBar()     
    numbers = np.arange(32)  
    currI = numbers[I]
    currO = numbers[O] 
    
    np.set_printoptions(precision=2, suppress=True)
    
    """print "currO vorher: ", currO
    print "currI vorher: ", currI"""
    
    ######
    t_i_u = h.i.CFU[currI] + cable_corr_i[currI][:,0] - t_glob_i 
    t_i_d = h.i.CFD[currI] + cable_corr_i[currI][:,1]    
    t_i_ud =  np.transpose([currI,t_i_u,t_i_d])
    
   # print t_i_ud
   
    #####
    t_o_u = h.o.CFU[currO] + cable_corr_o[currO][:,0]  - t_glob_o 
    t_o_d = h.o.CFD[currO] + cable_corr_o[currO][:,1]    
    t_o_ud =  np.transpose([currO,t_o_u,t_o_d])
    
    """print "#### OUTER "
    print "vorher", t_o_ud
    
    print "zpos ", (t_o_ud[:,1] - t_o_ud[:,2])/conv_factor_o
    
    print "meant",  (t_o_ud[:,1] + t_o_ud[:,2])*0.5
    
    print "#### INNER "
    print "vorher", t_i_ud
    
    print "zpos ", (t_i_ud[:,1] - t_i_ud[:,2])/conv_factor_i
    
    print "meant",  (t_i_ud[:,1] + t_i_ud[:,2])*0.5"""
        
    t_i_ud_new = t_i_ud[np.logical_and(
                    np.logical_and( (t_i_ud[:,1] - t_i_ud[:,2])/conv_factor_i > z_cuts_inner[0], (t_i_ud[:,1] - t_i_ud[:,2])/conv_factor_i < z_cuts_inner[1]),
                    np.logical_and( (t_i_ud[:,1] + t_i_ud[:,2])*0.5 > mt_cuts_inner[0], (t_i_ud[:,1] + t_i_ud[:,2])*0.5< mt_cuts_inner[1])
                    )]
    
    t_o_ud_new = t_o_ud[np.logical_and(
                    np.logical_and( (t_o_ud[:,1] - t_o_ud[:,2])/conv_factor_o > z_cuts_outer[0], (t_o_ud[:,1] - t_o_ud[:,2])/conv_factor_o < z_cuts_outer[1]),
                    np.logical_and( (t_o_ud[:,1] + t_o_ud[:,2])*0.5  > mt_cuts_outer[0], (t_o_ud[:,1] + t_o_ud[:,2])*0.5  < mt_cuts_outer[1])
                    )]      
                    
    ###########
    currI = t_i_ud_new[:,0].astype(int)
    currO = t_o_ud_new[:,0].astype(int)
    
    #print "nachher ", t_o_ud_new#
        
    z_o_pos = (- t_o_ud_new[:,1] + t_o_ud_new[:,2])/conv_factor_o         # # FIXME this is not a fixme. the minus is so that it is confrom with the fibre orietation
    z_o_pos = np.transpose([currO,z_o_pos])
    z_i_pos = (- t_i_ud_new[:,1] + t_i_ud_new[:,2])/conv_factor_i
    z_i_pos = np.transpose([currI,z_i_pos])
      
    #print "****** zpos_i", z_o_pos
    #print "meant ", (t_o_ud[:,1] + t_o_ud[:,2])*0.5
    
    return t_i_ud_new, t_o_ud_new, z_i_pos, z_o_pos, currI, currO      
##################################################################################################################################################################
##################################################################################################################################################################
def is_it_an_orphan(currI, currO):

    currI_c = np.array(currI)
    currO_c = np.array(currO)
    
    list_to_del = []

    for n,c in enumerate(currI):
        array_clustercheck = []
        #print " cluster: ", c 
        #c = c[0]
                
        if c == 31:            
            c0 = 1
            c2 = 29
        elif c == 0:
            c0 = 2
            c2 = 30
        else:
            c0 = c - 2
            c2 = c + 2
                                       
        array_clustercheck = np.arange(c0,c2+1)
                
        #check if this intersects with the other layer
        good_hits = list(set(array_clustercheck).intersection(currO))
                       
        if len(good_hits) == 0:
            #print "ja, inner!!!!!!!!!!!", c, n
            list_to_del.append(n)
            #print currI_c           
    
    currI_c = np.delete(currI, list_to_del)   
    
    list_to_del = []
           
    for m,c in enumerate(currO):
        array_clustercheck = []
        #print "1 cluster: ", c 
        #c = c[0]
                    
        if c == 31:            
            c0 = 1
            c2 = 29
        elif c == 0:
            c0 = 2
            c2 = 30
        else:
            c0 = c - 2
            c2 = c + 2
                                          
        array_clustercheck = np.arange(c0,c2+1)
                    
        #check if this intersects with the other layer
        good_hits = list(set(array_clustercheck).intersection(currI))
                    
            
        if len(good_hits) == 0:
            #print "ja, outer!!!!!!!!!!!", c
            list_to_del.append(m)
            #currO_c = np.delete(currO_c, m)    


    currO_c = np.delete(currO, list_to_del)   
    #print "oprph", currI_c
    #print currO_c
    
    return currI_c, currO_c
##################################################################################################################################################################


############################################################################################################## 
def rotate(b, ang):
    x,y = b
    return x*np.cos(ang)-y*np.sin(ang), x*np.sin(ang)+y*np.cos(ang)    
############################################################################################################## 
##############################################################################################################  
def getCornerMap():
    # left down corner of upper panel
    InnerCoordinatesLeftDown = [(103-5 ,-40-0.375/2), (103-5,-20), (103-5, 0+0.375), (103-5, 20+0.375*2)]
    OuterCoordinatesLeftDown = [(175-5,-70-0.375/2), (175-5,-35), (175-5,0+0.375), (175-5,35+0.375*2)]

    InnerCoordinatesLeftUp = [(103 ,-40-0.375/2), (103,-20), (103, 0+0.375), (103, 20+0.375*2)]
    OuterCoordinatesLeftUp = [(175,-70-0.375/2), (175,-35), (175,0+0.375), (175,35+0.375*2)]

    InnerCoordinatesRightDown = [(103-5 ,-40-0.375/2 + 20), (103-5,-20+20), (103-5, 0+0.375+20), (103-5, 20+0.375*2+20)]
    OuterCoordinatesRightDown = [(175-5,-70-0.375/2+35), (175-5,-35+35), (175-5,0+0.375+35), (175-5,35+0.375*2+35)]

    InnerCoordinatesRightUp = [(103 ,-40-0.375/2 + 20), (103,-20 + 20), (103, 0+0.375 + 20), (103, 20+0.375*2 + 20)]
    OuterCoordinatesRightUp = [(175,-70-0.375/2+35), (175,-35+35), (175,0+0.375+35), (175,35+0.375*2+35)]

    angles = [x*2*np.pi/8 for x in range(8)]
    InnerCoordinates = []
    OuterCoordinates = []

    for i,ang in enumerate(angles):
        for j,s,k,h in zip(InnerCoordinatesLeftDown, InnerCoordinatesLeftUp,InnerCoordinatesRightDown, InnerCoordinatesRightUp):
            InnerCoordinates.append([rotate(j, ang),rotate(s, ang),rotate(k, ang),rotate(h, ang)])
        for j,s,k,h in zip(OuterCoordinatesLeftDown, OuterCoordinatesLeftUp,OuterCoordinatesRightDown, OuterCoordinatesRightUp):
            OuterCoordinates.append([rotate(j, ang),rotate(s, ang),rotate(k, ang),rotate(h, ang)])

    return np.array(InnerCoordinates), np.array(OuterCoordinates)    
##############################################################################################################      
##############################################################################################################
def Get_Fibre_XYBar_Proj(): # calculate the projection of central bar points (x,y) of the bars onto the fibre layer
            mapI, mapO = getPositionMap()
            #fig, ax = plt.subplots()            
            rI = 84.5 # mm, radius inner fibre detector (center of fibre layer)
            rO = 147.0 # mm, radius outer fibre detector (center of fibre layer) 
            
            dist_arr_I = []
            dist_arr_O = []
            
            for i in mapI:
                 # distacne between point and circle:
                d_pc = abs(math.sqrt(i[0]**2 + i[1]**2) - rI)
                dist_arr_I.append(d_pc)
                #ax.scatter(i[0], i[1])
              
            for i in mapO:
                # distacne between point and circle:
                d_pc = abs(math.sqrt(i[0]**2 + i[1]**2) - rO)
                dist_arr_O.append(d_pc)
                #ax.scatter(i[0], i[1])
                
            dist_arr_I = np.array(dist_arr_I).reshape(32,1)
            dist_arr_O = np.array(dist_arr_O).reshape(32,1)    
            
            unit_vecs_I = - mapI / np.linalg.norm(mapI, axis = 1).reshape(32,1) 
            unit_vecs_O = - mapO / np.linalg.norm(mapO, axis = 1).reshape(32,1)
            
            fibre_points_I = mapI + unit_vecs_I*dist_arr_I
            fibre_points_O = mapO + unit_vecs_O*dist_arr_O    
            
            """for i in fibre_points_I:
                ax.scatter(i[0], i[1])
            for i in fibre_points_O:
                ax.scatter(i[0], i[1])            
            plt.show()"""
            
            return fibre_points_I, fibre_points_O

##############################################################################################################      
##############################################################################################################       
def prepare_PosMap(I,O, nI, nO):   
    posmaptestI, posmaptestO = getPositionMapTracking()
    posmaptestI, posmaptestO = posmaptestI[I], posmaptestO[O]
    I_hitpos = [i for i, x in enumerate(I) if x] # lists with numbers of bars with hits
    O_hitpos = [i for i, x in enumerate(O) if x]    

    #3 spalten mit x, y, barnr
    inneruntert = 21 #21
    outeruntert = 36  #36

    pmapI = np.zeros((inneruntert*nI,3))
    pmapO = np.zeros((outeruntert*nO,3))
            
    k=0
    for i in range(0,nI):
        for j in range(0,inneruntert):
            pmapI[k] = [posmaptestI[i][j][0], posmaptestI[i][j][1], I_hitpos[i]]
            k=k+1

    k=0
    for i in range(0,nO):
        for j in range(0,outeruntert):
            pmapO[k] = [posmaptestO[i][j][0], posmaptestO[i][j][1], O_hitpos[i]]
            k=k+1

    pmapI = np.array(pmapI)
    pmapO = np.array(pmapO) 
    
    return pmapI, pmapO 
##############################################################################################################      
##############################################################################################################


##############################################################################################################      
##############################################################################################################        
def prepare_PosMap_nottracking(I,O, nI, nO):   
    posmaptestI, posmaptestO = getPositionMap()
    posmaptestI, posmaptestO = posmaptestI[I], posmaptestO[O]
    I_hitpos = [i for i, x in enumerate(I) if x] # lists with numbers of bars with hits
    O_hitpos = [i for i, x in enumerate(O) if x]    

    pmapI = np.zeros((nI,3))
    pmapO = np.zeros((nO,3))
            
    k=0
    for i in range(0,nI):
            pmapI[k] = [posmaptestI[i][0], posmaptestI[i][1], int(I_hitpos[i])]
            k=k+1

    k=0
    for i in range(0,nO):
            pmapO[k] = [posmaptestO[i][0], posmaptestO[i][1], int(O_hitpos[i])]
            k=k+1

    pmapI = np.array(pmapI)
    pmapO = np.array(pmapO) 
    
    return pmapI, pmapO 
##############################################################################################################      
############################################################################################################## 
          
def prepare_CornerMap(I,O, nI, nO):  
    I_hitpos = [i for i, x in enumerate(I) if x] # lists with numbers of bars with hits
    O_hitpos = [i for i, x in enumerate(O) if x]      
        
    cornermapI, cornermapO = getCornerMap()
    I_hitpos = np.array(I_hitpos)
    if len(I_hitpos) > 0:
        cornermapI = cornermapI[I_hitpos]
    if len(O_hitpos) > 0:
        cornermapO = cornermapO[O_hitpos]
 
    cmapI = []
    for i in range(0,len(I_hitpos)):
        cmapI.append([cornermapI[i], I_hitpos[i]])   #TODO is append slow?
        
    cmapO = []
    for i in range(0,len(O_hitpos)):
        cmapO.append([cornermapO[i], O_hitpos[i]]) 

    cmapI = np.array(cmapI)
    cmapO = np.array(cmapO)
    if np.sum(O)== 0:
        cmap = cmapI               
    else:
        cmap =  np.vstack((cmapI,cmapO)) 

    return cmap         
##############################################################################################################      
##############################################################################################################  
def getPositionMapTracking():

    idI = np.arange(0,21.0,1.0)
    idO = np.arange(0,36.0,1.0)
    #print idI
    #print idO

    innerC = []
    outerC = []
    InnerCoordinates = []
    OuterCoordinates = []
    currentcoord = []
    InnerCoordinatesStart = [(103-2.5 ,-40-0.375/2), (103-2.5,-20), (103-2.5, 0+0.375), (103-2.5, 20+0.375*2)]
    OuterCoordinatesStart = [(175-2.5,-70-0.375/2), (175-2.5,-35), (175-2.5,0+0.375), (175-2.5,35+0.375*2)]     
    angles = [x*2*np.pi/8 for x in range(8)]

    for i,ang in enumerate(angles):
        for ind in range(0,len(InnerCoordinatesStart)):
            InnerCoordinates = []
            for k in idI:
                currentcoord = [InnerCoordinatesStart[ind][0], InnerCoordinatesStart[ind][1]+k] 
                InnerCoordinates.append(rotate(currentcoord, ang))
    
            innerC.append(InnerCoordinates)
        

    for i,ang in enumerate(angles):
        for ind in range(0,len(OuterCoordinatesStart)):
            OuterCoordinates = []
            for k in idO:
                currentcoord = [OuterCoordinatesStart[ind][0], OuterCoordinatesStart[ind][1]+k]
                OuterCoordinates.append(rotate(currentcoord, ang))

            outerC.append(OuterCoordinates)
            
    return np.array(innerC), np.array(outerC)
    
##############################################################################################################
##############################################################################################################  
def getPositionMap():
    InnerCoordinatesStart = [(103-2.5 ,-40-0.375/2+10), (103-2.5,-20+10), (103-2.5, 0+0.375+10), (103-2.5, 20+0.375*2+10)]
    OuterCoordinatesStart = [(175-2.5,-70-0.375/2+17.5), (175-2.5,-35+17.5), (175-2.5,0+0.375+17.5), (175-2.5,35+0.375*2+17.5)]
    angles = [x*2*np.pi/8 for x in range(8)]
    InnerCoordinates = []
    OuterCoordinates = []

    for i,ang in enumerate(angles):
        for j in InnerCoordinatesStart:
            #print j 
            InnerCoordinates.append(rotate(j, ang))
        for j in OuterCoordinatesStart:
            OuterCoordinates.append(rotate(j, ang))

    return np.array(InnerCoordinates), np.array(OuterCoordinates)
##############################################################################################################
