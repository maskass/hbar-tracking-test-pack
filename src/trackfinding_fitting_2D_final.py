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

from ROOT import TF3, Fit, Minuit2
import ROOT


##############################################################################################################
def get_hit_position2D(posmap, barnr, ch): # barnr is a list with the hitbars
                                                
    if len (barnr) > 1 :
        new_ch = []
        ch = np.array(ch)
        ch = ch.transpose()

        ch = list(ch)
        for i in barnr:
            for j in ch:

                if i == j[0]:
                    new_ch.append(0.5*(j[1]+j[2]))

        x = [posmap[i][0] for i in barnr]
        #print x
        y = [posmap[i][1] for i in barnr]  
        #print y
        
        x = np.dot(new_ch,x) / np.sum(new_ch)
        y = np.dot(new_ch,y) / np.sum(new_ch)   
        #print "x", x
        #print "y", y             
        #x = sum(x)/len(x)
        #y = sum(y)/len(y)
        #z = sum(z)/len(z)
    
    else:
        #print "else!"
        x = posmap[barnr][0][0]
        y = posmap[barnr][0][1]

    return x,y
##############################################################################################################
####################################################################################
####################################################################################
def Calc_xy_error_fit(barnr, layer):
    #print "trololo"
    #InnerCoordinatesStart = [(103-2.5 ,-40-0.375/2+10), (103-2.5,-20+10), (103-2.5, 0+0.375+10), (103-2.5, 20+0.375*2+10)]
    #OuterCoordinatesStart = [(175-2.5,-70-0.375/2+17.5), (175-2.5,-35+17.5), (175-2.5,0+0.375+17.5), (175-2.5,35+0.375*2+17.5)]
    #angles = [x*2*np.pi/8 for x in range(8)]
    #print angles
    
    tot_ang = []
    rep = [1.0]*4
    for i in range(8):
        for j in rep:
            angle = i*2*np.pi/8
            tot_ang. append(angle)


    if layer == 1:
        b_ = 35.0 # outer
    else:
        b_ = 20.0
        
    sigmau2 = b_*b_/12.0
    sigmav2 = 5.0*5.0/12.0
    
    V1 = np.array([[sigmau2, 0.0],[0.0, sigmav2]])
    for i in tot_ang:
        alpha = i
        rota = np.array([[math.cos(alpha),-math.sin(alpha)],[math.sin(alpha), math.cos(alpha)]])
        rota_t = np.array([[math.cos(alpha),math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])
        mat = np.dot(rota,np.dot(V1,rota_t))
        #print math.sqrt(mat[0][0]), math.sqrt(mat[1][1])

    
    return math.sqrt(mat[0][0]), math.sqrt(mat[1][1])

############################################################################################################## 
##############################################################################################################
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
##############################################################################################################
def dist_p_line2d(x,y,p):    
    dist = abs(p[0]*x + p[1]*y + p[2])/(math.sqrt(p[0]*p[0] + p[1]*p[1]))
    return dist*dist

##############################################################################################################
def line_2d(x,p):
    y = -p[0]/p[1]*x - p[2]/p[1]    
    return y     
##############################################################################################################   
def do_le_2Dfitting(hitcoll, I, O, vertex, vertex_err, innerA, outerA):
    
    posmapI, posmapO = getPositionMap()
    cornermapI, cornermapO = getCornerMap()
    cmapI = np.array(cornermapI)
    cmapO = np.array(cornermapO)   
        
    trackscoll = []
    chisqcoll = []
    numbertr = len(hitcoll)

    rltrack = [] # helper list to distinguish between tracks with cluster and normal ones    



    final_lines_params = []
    final_lines_points = []
    
    bgox = vertex[0]
    bgoy = vertex[1]
    
    vert_err_x = vertex_err[0]
    vert_err_y = vertex_err[1]
    
    # now loop over all resulting tracks, saved in 'cluster'
    for ctr in hitcoll:
        #print "--- AFTER CLUSTER CHECK: ", ctr
        trackpoints = []
        trackpoints.append([bgox, bgoy, vert_err_x, vert_err_y])
        #print " "
        for hit in ctr:
            if type(hit) is list:
                #print "hit        ", len(hit)
                #barnr = []
                if type(hit[1]) is list: # cluster
                    barnr = hit[1]
                    #print "1", barnr
                else:                    # no cluster
                    barnr = hit[1:]
                    #print "2", barnr
                #print "len: ", len(barnr), "barnr: ", barnr
                
                if int(hit[0]) == 0:
                    #get x,y,z pos of cluster-> mean val ->add to trackpos
                    x,y = get_hit_position2D(posmapI, barnr, innerA)
                    xerr, yerr = Calc_xy_error_fit(barnr,hit[0])
                                   
                if int(hit[0]) == 1:
                    #get x,y,z pos of cluster-> mean val ->add to trackpos
                    x,y = get_hit_position2D(posmapO, barnr, outerA)
                    xerr, yerr = Calc_xy_error_fit(barnr, hit[0])
                    #print "errrrrrr", xerr, yerr
                    
                trackpoints.append([x,y,xerr,yerr])
            
        trackscoll.append(trackpoints) 
        
    number_of_tracks = len(trackscoll)   
    
    #print "trackscoll: "
    #for i in trackscoll:
    #    print i
    
    #save_results = []
             
    for n, currtrack in enumerate(trackscoll):
           
        currtrack =  np.array(currtrack)   
        samplesize = len(trackscoll[n])
  
        points = currtrack[:,:2]
        #print "points: "
        #print points
        errors = currtrack[:,2:]
        #print " errors", errors            

        looppoints = [points[:,[0,1]]]
        looperrors = [errors[:,[0,1]]]  

        for n,(p,e) in enumerate(zip(looppoints, looperrors)):
            
                #print "pppppppppppppppp ", n, p
            
                pars = np.array([1.0,1.0,1.0])    # start values
                step = [0.01,0.01,0.01]
                minimizer = ROOT.Minuit2.Minuit2Minimizer(0)

                errorfct = chi2fct_track2d()

                errorfct.SetVals(samplesize, p, e)
                minimizer.SetFunction(errorfct)

                minimizer.SetMaxFunctionCalls(10000000)
                minimizer.SetMaxIterations(1000000)
                minimizer.SetTolerance(0.001)

                minimizer.SetVariable(0,"par0", pars[0], step[0])
                minimizer.SetVariable(1,"par1", pars[1], step[1])
                minimizer.SetVariable(2,"par2", pars[2], step[2])

               
                retval_min = minimizer.Minimize()
                
                minval = minimizer.MinValue()

                #minimizer.PrintResults()
                stat = minimizer.Status() # status 0 means that the errors of the fit are ok. status 1 means, you cant trust the errorbars
                #print "ooooooooooooooooooooooooooooo ", minval, minimizer.Status() #"""

                ret_params = minimizer.X()
                ret_errors = minimizer.Errors()
                final_params = []
            
           
                for i in range(0,pars.size):
                    final_params.append(ret_params[i])
                #final_errors = np.zeros(3)
                final_lines_params.append(final_params)   
                final_lines_points.append(p)

    """angle = [999,999,999]
    if len(hitcoll) == 2:
        #angle = angle_2_lines(final_lines_params, trackscoll)
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        angle = angle_btw_lines(final_lines_params, trackscoll)
        print "------------------------------------------------------------------------------------------"
        print "ANGLE: ", angle, "degree" 
        print "------------------------------------------------------------------------------------------"

    if len(hitcoll) == 3:
        angle = angle_3_lines(final_lines_params, trackscoll)
        print "------------------------------------------------------------------------------------------"
        print "ANGLES: ", angle, "degree" 
        print "SUM: ", sum(angle)
        print "------------------------------------------------------------------------------------------" """

    return final_lines_params, final_lines_points, trackscoll #, angle
##############################################################################################################
class chi2fct_track2d( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        #print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.samplesize = 0
        self.points = np.zeros(10)
        self.errors = np.zeros(10)

    def NDim( self ):
        #print 'PYTHON NDim called'
        return 3

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
        
        for i in range(0,self.samplesize):
            chisq = chisq + dist_p_line2d(self.points[i][0],self.points[i][1], pars)/errsum

        #print chisq

        return chisq             

##############################################################################################################
##############################################################################################################
##############################################################################################################
def set_list_intersection(set_list):
  if not set_list:
    return set()
  result = set_list[0]
  for s in set_list[1:]:
    result &= s
  return result
##############################################################################################################    
def find_consecutive_nr(numbers):    
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
            nwclusters.remove(c)
        if 0 in c:
            tempc.append(c)
            nwclusters.remove(c)
    
    #print "clusters temp::::::::::: ", tempc                
    if len(tempc) >= 1: # if both clusters with 31 and 0 exist, merge them
        newcluster = [item for sublist in tempc for item in sublist]
        nwclusters.append(newcluster)
     
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
            nwclusters.remove(c)
        if 0 in c:
            tempc.append(c)
            nwclusters.remove(c)
    
    #print "clusters temp::::::::::: ", tempc                
    if len(tempc) >= 1: # if both clusters with 31 and 0 exist, merge them
        newcluster = [item for sublist in tempc for item in sublist]
        nwclusters.append(newcluster)

    return yesno, nwclusters    # bool: found cluster: yes or no?     

##############################################################################################################  
def track_finding_circle(pmapI, cmap, currI, currO, circle, plot_trackfinding, bgo_hit_bounds, innerCFs, outerCFs, innerzPos, outerzPos, BGOXPos, BGOYPos):  # with density of lines added
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
        yesnoI, cI = find_consecutive_nr(currI)  # TODO this is only used for sort tracks function below... check if really needed
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
        FWHM_it = 0.8 # ns
        FWHM_ip = 59 # mm
        FWHM_ip = 30 # mm
        #print "FWHM inner z-pos: ", FWHM_ip
        #print "FWHM inner mtd: ", FWHM_it
        cl_i_merge  = []   
        added_bars = []
        for cluI in clusI:
            cl_i_cfs = innerCFs[np.nonzero(np.in1d(innerCFs[:,0],np.array(cluI)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluI) == 1:
                cl_i_merge.append([0, [int(cluI[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_i_cfs):
                    for m,j in enumerate(cl_i_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = innerzPos[innerzPos[:,0]==bar1][0][1] - innerzPos[innerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle                            
                            meantime1 = 0.5*((innerCFs[innerCFs[:,0]==bar1])[0][1] + (innerCFs[innerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((innerCFs[innerCFs[:,0]==bar2])[0][1] + (innerCFs[innerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2
                            #print meantime1
                            #print meantime2
                            #print meantimediff
                            all_combs.append([0,[bar1,bar2], zposdiff, meantimediff])
                            #if [0, [bar1]] not in cl_i_merge: # and bar1 not in added_bars:
                            #        cl_i_merge.append([0, [bar1]])
                            #if [0, [bar2]] not in cl_i_merge: # and bar2 not in added_bars:
                            #        cl_i_merge.append([0, [bar2]])    

                
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                    
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_ip and abs(ac[3]) < FWHM_it:
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        #print "yas", ac
                        #print pot_new_cl
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
                                          
                #print "pot all : "
                for i in all_pot_cls:
                    ccc = [0,i]
                    if not ccc in cl_i_merge:
                        cl_i_merge.append(ccc)
                    #print i

        #print "inner after neighbour z-pos check: "
        #print cl_i_merge
       
        

        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"            
        FWHM_op = 73 # mm   
        FWHM_op = 50 # mm        
        FWHM_ot = 0.9 #ns
        #print "FWHM outer z-pos: ", FWHM_op    
        #print "FWHM outer mtd: ", FWHM_ot          
        cl_o_merge = []
        added_bars = []
        for cluO in clusO:
            cl_o_cfs = outerCFs[np.nonzero(np.in1d(outerCFs[:,0],np.array(cluO)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluO) == 1:
                cl_o_merge.append([1, [int(cluO[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_o_cfs):
                    for m,j in enumerate(cl_o_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = outerzPos[outerzPos[:,0]==bar1][0][1] - outerzPos[outerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle
                            #print outerCFs
                            meantime1 = 0.5*((outerCFs[outerCFs[:,0]==bar1])[0][1] + (outerCFs[outerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((outerCFs[outerCFs[:,0]==bar2])[0][1] + (outerCFs[outerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2                            
                            all_combs.append([1,[bar1,bar2], zposdiff, meantimediff])
                            #if [1, [bar1]] not in cl_o_merge: # and bar1 not in added_bars:
                            #        cl_o_merge.append([1, [bar1]])
                            #if [1, [bar2]] not in cl_o_merge: #  and bar2 not in added_bars:
                            #        cl_o_merge.append([1, [bar2]])    
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_op and abs(ac[3]) < FWHM_ot:
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
                                                                                                          
        #print "outer after neighbour z-pos check: "
        #print cl_o_merge

        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
       
        # if merge_cl_i still contains neighbouring bars, check again: ... while? recursive?
        # (cluster with more than two bars hit)


        # if a cluster is recognized as two hits, they are not yet added.... so add them now:        
        cl_i_list = [item for sublist in [c[1] for c in cl_i_merge] for item in sublist] # list of all bars in clusters (flatten cluster list)
        all_i_list = [item for sublist in  clusI for item in sublist] # list of ALL bars
        
        
        cl_o_list = [item for sublist in [c[1] for c in cl_o_merge] for item in sublist]
        all_o_list = [item for sublist in clusO for item in sublist] # list of all bars (flatten cluster list)
        
        #print "ALL inner: ", cl_i_list
        #print "ALL outer: ", cl_o_list
        
        for b in all_i_list:
            if b not in cl_i_list:
                cl_i_merge.append([0,[b]])
                
                
        for b in all_o_list:
            if b not in cl_o_list:
                cl_o_merge.append([1,[b]])
        
        #print "cl_i_merge ", cl_i_merge      
        #print "cl_o_merge ", cl_o_merge                             
        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
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
        len_bgo_list = 20.0
        
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
                                # TODO if the track contains a cluster with a bar e.g. 6, then dont add [1,6] again!!
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
                #print h
        #print "##################################"    
         
        if len(hitcollection2) == 0:
            return 999, 999, 999     
            
        final_tracks = new_select_tracks_for_2D(hitcollection2, currI, currO) 
        
        if final_tracks == 999:
            return 999, 999, 999
        
        return final_tracks, line_points_collection, bgo_points_coll
    
    
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
    
    final_tracks = new_select_tracks_for_2D(hitcollection2,currI, currO)
    
    #plot_all_bars_2D(circle,cornermapI,cornermapO,union_cl_I, union_cl_O, line_points_collection)
    
    
    return final_tracks, line_points_collection #"""

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
# calculate the error on x and y coordinate depending on bar number and layer
def Calc_xy_error(barnr, layer):

    tot_ang = []
    rep = [1.0]*4
    for i in range(8):
        for j in rep:
            angle = i*2*np.pi/8
            tot_ang. append(angle)
    
    if layer == 1:
        b_ = 35.0 # outer
    else:
        b_ = 20.0
        
    sigmau2 = b_*b_/12.0
    sigmav2 = 5.0*5.0/12.0
    
    V1 = np.array([[sigmau2, 0.0],[0.0, sigmav2]])
    for i in tot_ang:
        alpha = i
        rota = np.array([[math.cos(alpha),-math.sin(alpha)],[math.sin(alpha), math.cos(alpha)]])
        rota_t = np.array([[math.cos(alpha),math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])
        mat = np.dot(rota,np.dot(V1,rota_t))
        #print math.sqrt(mat[0][0]), math.sqrt(mat[1][1])
   
    return math.sqrt(mat[0][0]), math.sqrt(mat[1][1])       
##################################################################################################################################
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
########################################################################################################################################################
########################################################################################################################################################
# flatten a list
def flatten(l): return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]

########################################################################################################################################################
def track_finding_polygons_tpx(pmapI, cmap, currI, currO, polygon_hitlist, plot_trackfinding, bgo_hit_bounds, innerCFs, outerCFs, innerzPos, outerzPos, verbose_): #, BGOXPos, BGOYPos):  # with density of lines added

        nI = len(currI)
        nO = len(currO)
       
        cmapI = cmap[:nI,:]
        cmapO = cmap[nI:,:] 
        
        # use this for no unification of rects
        yesnoI, cI = find_consecutive_nr(currI)  # TODO this is only used for sort tracks function below... check if really needed
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
        FWHM_it = 0.8 # ns
        FWHM_ip = 59 # mm
        FWHM_ip = 30 # mm
        #print "FWHM inner z-pos: ", FWHM_ip
        #print "FWHM inner mtd: ", FWHM_it
        cl_i_merge  = []   
        added_bars = []
        for cluI in clusI:
            cl_i_cfs = innerCFs[np.nonzero(np.in1d(innerCFs[:,0],np.array(cluI)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluI) == 1:
                cl_i_merge.append([0, [int(cluI[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_i_cfs):
                    for m,j in enumerate(cl_i_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = innerzPos[innerzPos[:,0]==bar1][0][1] - innerzPos[innerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle                            
                            meantime1 = 0.5*((innerCFs[innerCFs[:,0]==bar1])[0][1] + (innerCFs[innerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((innerCFs[innerCFs[:,0]==bar2])[0][1] + (innerCFs[innerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2
                            #print meantime1
                            #print meantime2
                            #print meantimediff
                            all_combs.append([0,[bar1,bar2], zposdiff, meantimediff])
                            #if [0, [bar1]] not in cl_i_merge: # and bar1 not in added_bars:
                            #        cl_i_merge.append([0, [bar1]])
                            #if [0, [bar2]] not in cl_i_merge: # and bar2 not in added_bars:
                            #        cl_i_merge.append([0, [bar2]])    

                
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                    
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_ip and abs(ac[3]) < FWHM_it:
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        #print "yas", ac
                        #print pot_new_cl
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
                                          
                #print "pot all : "
                for i in all_pot_cls:
                    ccc = [0,i]
                    if not ccc in cl_i_merge:
                        cl_i_merge.append(ccc)
                    #print i

        #print "inner after neighbour z-pos check: "
        #print cl_i_merge
       
        

        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"            
        FWHM_op = 73 # mm   
        FWHM_op = 50 # mm        
        FWHM_ot = 0.9 #ns
        #print "FWHM outer z-pos: ", FWHM_op    
        #print "FWHM outer mtd: ", FWHM_ot          
        cl_o_merge = []
        added_bars = []
        for cluO in clusO:
            cl_o_cfs = outerCFs[np.nonzero(np.in1d(outerCFs[:,0],np.array(cluO)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluO) == 1:
                cl_o_merge.append([1, [int(cluO[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_o_cfs):
                    for m,j in enumerate(cl_o_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = outerzPos[outerzPos[:,0]==bar1][0][1] - outerzPos[outerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle
                            #print outerCFs
                            meantime1 = 0.5*((outerCFs[outerCFs[:,0]==bar1])[0][1] + (outerCFs[outerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((outerCFs[outerCFs[:,0]==bar2])[0][1] + (outerCFs[outerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2                            
                            all_combs.append([1,[bar1,bar2], zposdiff, meantimediff])
                            #if [1, [bar1]] not in cl_o_merge: # and bar1 not in added_bars:
                            #        cl_o_merge.append([1, [bar1]])
                            #if [1, [bar2]] not in cl_o_merge: #  and bar2 not in added_bars:
                            #        cl_o_merge.append([1, [bar2]])    
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_op and abs(ac[3]) < FWHM_ot:
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
                                                                                                          
        #print "outer after neighbour z-pos check: "
        #print cl_o_merge

        
       
        # if merge_cl_i still contains neighbouring bars, check again: ... while? recursive?
        # (cluster with more than two bars hit)


        # if a cluster is recognized as two hits, they are not yet added.... so add them now:        
        cl_i_list = [item for sublist in [c[1] for c in cl_i_merge] for item in sublist] # list of all bars in clusters (flatten cluster list)
        all_i_list = [item for sublist in  clusI for item in sublist] # list of ALL bars
        
        
        cl_o_list = [item for sublist in [c[1] for c in cl_o_merge] for item in sublist]
        all_o_list = [item for sublist in clusO for item in sublist] # list of all bars (flatten cluster list)
        
        #print "ALL inner: ", cl_i_list
        #print "ALL outer: ", cl_o_list
        
        for b in all_i_list:
            if b not in cl_i_list:
                cl_i_merge.append([0,[b]])
                
                
        for b in all_o_list:
            if b not in cl_o_list:
                cl_o_merge.append([1,[b]])
        
        #print "cl_i_merge ", cl_i_merge      
        #print "cl_o_merge ", cl_o_merge                             
        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
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
        len_bgo_list = 30
        
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
                                # TODO if the track contains a cluster with a bar e.g. 6, then dont add [1,6] again!!
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
                #print h
        #print "##################################"    
         
        if len(hitcollection2) == 0:
            return 999, 999, 999     
            
        final_tracks = new_select_tracks_for_2D(hitcollection2, currI, currO, verbose_) 
        
        if final_tracks == 999:
            return 999, 999, 999
        
        return final_tracks, line_points_collection, bgo_points_coll
    
    
############################################################################################################################################################################ 
############################################################################################################################################################################ 
def track_finding_polygons(pmapI, cmap, currI, currO, polygon_hitlist, plot_trackfinding, bgo_hit_bounds, innerCFs, outerCFs, innerzPos, outerzPos, verbose_): #, BGOXPos, BGOYPos):  # with density of lines added
        ### First: Find Clusters ###
        #print "####################"
        #print "inner bars hit: ", cmap[:nI,1] # corner map, all rects with corner points plus nr
        #print "outer bars hit: ", cmap[nI:,1]     
        #plt.show()        
        #print "pmapI ", pmapI
        
        print "1"
        
        #print "Outer CFs "
        #print outerCFs

        nI = len(currI)
        nO = len(currO)
       
        cmapI = cmap[:nI,:]
        cmapO = cmap[nI:,:] 
        
        # use this for no unification of rects
        yesnoI, cI = find_consecutive_nr(currI)  # TODO this is only used for sort tracks function below... check if really needed
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
        FWHM_it = 0.8 # ns
        FWHM_ip = 59 # mm
        FWHM_ip = 30 # mm
        #print "FWHM inner z-pos: ", FWHM_ip
        #print "FWHM inner mtd: ", FWHM_it
        cl_i_merge  = []   
        added_bars = []
        for cluI in clusI:
            cl_i_cfs = innerCFs[np.nonzero(np.in1d(innerCFs[:,0],np.array(cluI)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluI) == 1:
                cl_i_merge.append([0, [int(cluI[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_i_cfs):
                    for m,j in enumerate(cl_i_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = innerzPos[innerzPos[:,0]==bar1][0][1] - innerzPos[innerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle                            
                            meantime1 = 0.5*((innerCFs[innerCFs[:,0]==bar1])[0][1] + (innerCFs[innerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((innerCFs[innerCFs[:,0]==bar2])[0][1] + (innerCFs[innerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2
                            #print meantime1
                            #print meantime2
                            #print meantimediff
                            all_combs.append([0,[bar1,bar2], zposdiff, meantimediff])
                            #if [0, [bar1]] not in cl_i_merge: # and bar1 not in added_bars:
                            #        cl_i_merge.append([0, [bar1]])
                            #if [0, [bar2]] not in cl_i_merge: # and bar2 not in added_bars:
                            #        cl_i_merge.append([0, [bar2]])    

                
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                    
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_ip and abs(ac[3]) < FWHM_it:
                        bar1 = ac[1][0]
                        bar2 = ac[1][1]
                        #print "yas", ac
                        #print pot_new_cl
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
                                          
                #print "pot all : "
                for i in all_pot_cls:
                    ccc = [0,i]
                    if not ccc in cl_i_merge:
                        cl_i_merge.append(ccc)
                    #print i

        #print "inner after neighbour z-pos check: "
        #print cl_i_merge
       
        

        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"            
        FWHM_op = 73 # mm   
        FWHM_op = 50 # mm        
        FWHM_ot = 0.9 #ns
        #print "FWHM outer z-pos: ", FWHM_op    
        #print "FWHM outer mtd: ", FWHM_ot          
        cl_o_merge = []
        added_bars = []
        for cluO in clusO:
            cl_o_cfs = outerCFs[np.nonzero(np.in1d(outerCFs[:,0],np.array(cluO)))[0]] # cfs of clusters
            #nr_of_comb = len(cluI) - 1 # check all neighbouring bars 
            #cl_i_m_d = np.zeros((nr_of_comb,4))  #array to save [bar1, bar2, mtd, zposdiff] of neighbours      
            if len(cluO) == 1:
                cl_o_merge.append([1, [int(cluO[0])]])
            else:      
                all_combs = []          
                for n,i in enumerate(cl_o_cfs):
                    for m,j in enumerate(cl_o_cfs):
                        if n < m: # and abs(m -n) in [1,31] :
                            bar1 = int(i[0]) # bar number 1
                            bar2 = int(j[0]) # bar number 2
                            zposdiff = outerzPos[outerzPos[:,0]==bar1][0][1] - outerzPos[outerzPos[:,0]==bar2][0][1] # diff of zpos... will be around zero if hit frorm same particle
                            #print outerCFs
                            meantime1 = 0.5*((outerCFs[outerCFs[:,0]==bar1])[0][1] + (outerCFs[outerCFs[:,0]==bar1])[0][2])
                            meantime2 = 0.5*((outerCFs[outerCFs[:,0]==bar2])[0][1] + (outerCFs[outerCFs[:,0]==bar2])[0][2])
                            meantimediff = meantime1 - meantime2                            
                            all_combs.append([1,[bar1,bar2], zposdiff, meantimediff])
                            #if [1, [bar1]] not in cl_o_merge: # and bar1 not in added_bars:
                            #        cl_o_merge.append([1, [bar1]])
                            #if [1, [bar2]] not in cl_o_merge: #  and bar2 not in added_bars:
                            #        cl_o_merge.append([1, [bar2]])    
                 
                #print "all combinations: "
                #for a in all_combs:
                #    print a
                
                all_pot_cls = []
                pot_new_cl = [] 
                for ac in all_combs:   
                    #pot_new_cl = [] 
                    if abs(ac[2]) < FWHM_op and abs(ac[3]) < FWHM_ot:
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
                                                                                                          
        #print "outer after neighbour z-pos check: "
        #print cl_o_merge

        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
       
        # if merge_cl_i still contains neighbouring bars, check again: ... while? recursive?
        # (cluster with more than two bars hit)


        # if a cluster is recognized as two hits, they are not yet added.... so add them now:        
        cl_i_list = [item for sublist in [c[1] for c in cl_i_merge] for item in sublist] # list of all bars in clusters (flatten cluster list)
        all_i_list = [item for sublist in  clusI for item in sublist] # list of ALL bars
        
        
        cl_o_list = [item for sublist in [c[1] for c in cl_o_merge] for item in sublist]
        all_o_list = [item for sublist in clusO for item in sublist] # list of all bars (flatten cluster list)
        
        #print "ALL inner: ", cl_i_list
        #print "ALL outer: ", cl_o_list
        
        for b in all_i_list:
            if b not in cl_i_list:
                cl_i_merge.append([0,[b]])
                
                
        for b in all_o_list:
            if b not in cl_o_list:
                cl_o_merge.append([1,[b]])
        
        #print "cl_i_merge ", cl_i_merge      
        #print "cl_o_merge ", cl_o_merge                             
        #print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
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
        len_bgo_list = 0.7*len(polygon_hitlist)
        
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
                    for i in range(0,int(len_bgo_list)):  # draw XX points from the circle / polygons aka all pmt pixel with a hit. make it dependent of the size of the hit shape
                    
                        #p1 = draw_pt_from_circle(circle.radius,[BGOXPos,BGOYPos])
                        p1 = draw_point_from_bgohit(polygon_hitlist, bgo_hit_bounds)
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
                                # TODO if the track contains a cluster with a bar e.g. 6, then dont add [1,6] again!!
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
                #print h
        #print "##################################"    
         
        if len(hitcollection2) == 0:
            return 999, 999, 999     
            
        final_tracks = new_select_tracks_for_2D(hitcollection2, currI, currO, verbose_) 
        
        if final_tracks == 999:
            return 999, 999, 999
        
        return final_tracks, line_points_collection, bgo_points_coll
    
    
############################################################################################################## 
##############################################################################################################  
##############################################################################################################   
def new_select_tracks_for_2D(hitcollection, currI, currO, verbose_): # has the density of lines added!

    if verbose_:
        print "select tracks ############################################################################"


    # sort according to probability        
    hitcollection = list(reversed(sorted(hitcollection, key=itemgetter(-1))))

    if verbose_:
        print "------------------------------------------------------------------------------------------"  
                             
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
                            
                        """elif trh in currI and len(currI) != 0 and len(currO) == 0:
                            print "in 2 "
                            finalhitcoll.append(tr) # but here we append with line count
                            currI = currI[currI != trh]
                            
                        elif trho in currO and len(currO) != 0 and len(currI) == 0:
                            print "in 3 "
                            finalhitcoll.append(tr) # but here we append with line count
                            currO = currO[currO != trho]
                        
                        elif len(currI) == 0 and len(currO) == 0: 
                            print "in 4 "
                            break"""
                        
                        
                        
                        

     
    #print "currI", currI
    #print "currO", currI 

    
    
    #copy finalhitcoll
    """finalhitcoll_copy = list(finalhitcoll)
      
    ultimatehitcoll = []
    ultimatehitcoll.append(finalhitcoll[0]) # append track with highest probabiltiy
    del finalhitcoll_copy[0] # remove that element from the copy
    
                
    for t in finalhitcoll_copy:  
        yes = False   # check all other tracks for overlaps  
        if t in ultimatehitcoll:
            continue
             
        if not t in ultimatehitcoll:
            track = t[:-1]
            yes = False
            for hit in track:
                #print "hit ", hit 
                for utr in ultimatehitcoll:
                    #print hit, len(hit[1])
                    if len(hit[1])> 1:   # sort out tracks that have clusters and are already added
                        #print "in ", [hit[1][0]], [hit[1][1]]
                        #print ultimatehitcoll
                        if [hit[0],[hit[1][0]]] in utr or [hit[0],[hit[1][1]]] in utr:
                            #print "yes!"
                            yes = True
                            break
                                                   
                    if hit in utr:
                        yes = True
                        break

            if yes == False:
                ultimatehitcoll.append(t)        
    
      
    # remove dublicates:  # works only with ints, I think..
    #
    
    # remove dublicates the itertool way:
    #finalhitcoll = [s for s in finalhitcoll if len(s[0]) in [2] ]
    #finalhitcoll.sort()
    
 
    #finalhitcoll = list(finalhitcoll for finalhitcoll,_ in itertools.groupby(finalhitcoll))
    
    #finalhitcoll = [s for s in finalhitcoll if len(s[0]) in [2] ]
    #finalhitcoll.sort()
    #finalhitcoll = list(finalhitcoll for finalhitcoll,_ in itertools.groupby(finalhitcoll))"""
    

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






#########################################################################################################################################################################
#########################################################################################################################################################################
#                                                                           FUNCTION CEMETERY
#                                                   ......................some might come back ..................
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################

def track_selection_fitting_2D(central_hit, central_hit_err, point_collections):

    print " "

    print "-------------------------------- 3D TRACK SELECTION START --------------------------------"
    final_lines_params = [] # store parameter of fitted lines 
    for n, currcoll in enumerate(point_collections):
                print " "
                points = []
                errors = []
                for i in currcoll:
                    print i
                    points.append([i[0][0], i[1][0]]) # [z, x/y]
                    errors.append([i[0][1], i[1][1]])
       
                    
                points.append(central_hit)
                errors.append(central_hit_err)

                samplesize = len(points)

                pars = np.array([1.0,1.0])    # start values
                step = [0.01,0.01]
                minimizer = ROOT.Minuit2.Minuit2Minimizer(0)

                errorfct = chi2fct2d()

                errorfct.SetVals(samplesize, points, errors)
                minimizer.SetFunction(errorfct)

                minimizer.SetMaxFunctionCalls(10000000)
                minimizer.SetMaxIterations(1000000)
                minimizer.SetTolerance(0.0001)

                minimizer.SetVariable(0,"par0", pars[0], step[0])
                minimizer.SetVariable(1,"par1", pars[1], step[1])
                #minimizer.SetVariable(2,"par2", pars[2], step[2])

               
                retval_min = minimizer.Minimize()
                
                minval = minimizer.MinValue()

                minimizer.PrintResults()
                stat = minimizer.Status() # status 0 means that the errors of the fit are ok. status 1 means, you cant trust the errorbars
                print "ooooooooooooooooooooooooooooo ", minval, minimizer.Status()#"""

                ret_params = minimizer.X()
                ret_errors = minimizer.Errors()
                final_params = []
            
           
                for i in range(0,pars.size):
                    final_params.append(ret_params[i])
                #final_errors = np.zeros(3)
                final_lines_params.append(final_params)   
                #final_lines_points.append(p)
                
    return final_lines_params                

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
##############################################################################################################   

######################################################################################################################################################################### 
######################################################################################################################################################################### 
def create_all_polygons(fibres_no_hits, z_fibres, bar_coords, curr_bars, width_x, width_y):
    #print "haha"
    
    #print len(z_fibres), len(fibres_no_hits)
    #print len(curr_bars), len(bar_coords)
    
    poli_side_list = []
    poli_top_list = []
    helper_list_combinations = []
    
    for iz, fi in zip(z_fibres[fibres_no_hits], fibres_no_hits):
        for bar, point in zip(curr_bars, bar_coords):
            #print bar, fi
            if  [bar, fi] in helper_list_combinations:
                continue
                
            #print bar, fi
            xerr = width_x[bar]*0.5
            yerr = width_y[bar]*0.5
            poli_side = Polygon([(iz + 2.0, point[1]-yerr), (iz + 2.0, point[1]+yerr), (iz - 2.0, point[1]+yerr), (iz - 2.0, point[1]-yerr)])
            poli_top = Polygon([(iz + 2.0, point[0]-xerr), (iz + 2.0, point[0]+xerr), (iz - 2.0, point[0]+xerr), (iz - 2.0, point[0]-xerr)])
            poli_side_list.append(poli_side)
            poli_top_list.append(poli_top)
            helper_list_combinations.append([bar, fi])

    #print helper_list_combinations

    return poli_side_list, poli_top_list             

######################################################################################################################################################################### 
#########################################################################################################################################################################  
def rotate_line_through(polygons_no_hit_i, polygons_no_hit_o, plotyes, BGO_pol, pol_inner, pol_outer, poli_p, polo_p, max_, min_, z_bgo_max, z_bgo_min):
    
    #c_points = points_on_circ(300.0, 200)
    
    # pol_inner/outer:  array with elements of structure: [polygon, bar, fibre]
    # polo/polo_p: [z pos, x/y (depending whether side or top), bar, fibre]
    
    
    #print "FIBRES WITH HITS, INNER:", pol_inner
    #print "FIBRES WITH HITS, OUTER:", pol_outer
    
    #print " "
    
    #print "FIBRES WITH HITS, INNER 2:", poli_p
    #print "FIBRES WITH HITS, OUTER 2:", polo_p
    
    line_collection = [] 
    point_coll = []

    nr_bgo_pts = 100

    for i in range(0, nr_bgo_pts):
        #print "bgo point loop..."
       
        p1 = draw_point_from_bgohit([BGO_pol], [[100,-100], [100, -100]])  
        #ax.plot(p1[0],p1[1])
        #print p1
        
        
        
        for polo in polo_p: # these are points
            #print polo[2], p1
            linecoll_elem = []
            #print "polo loop..."
            
            line = [p1, polo[0]] # point and coords
            path = shapely.geometry.LineString(line)
                       
            linecoll_elem.append([1,polo[1], polo[2]]) # [layer, bar, fibre]
            
            
            for poli in pol_inner:
            #    print "here!"
                intersection = path.intersection(poli[0])
                
                """nohit_intersec_inner = False
                for nohit_i in polygons_no_hit_i:
                    intersec = path.intersection(nohit_i)
                    if intersec.type == 'LineString':
                        nohit_intersec_inner = True
                        break
                        
                        
                        
                nohit_intersec_outer  = False
                for nohit_o in polygons_no_hit_o:
                    intersec = path.intersection(nohit_o)
                    if intersec.type == 'LineString':
                        nohit_intersec_outer = True
                        break
                
                if nohit_intersec_outer == True or nohit_intersec_inner == True:
                    continue"""
                    
                #print intersection.type
                if intersection.type == 'LineString':  
                    #plt.plot([p1[0], c[0]], [p1[1], c[1]], marker = '_')
                    #print "int line", intersection
                    #print poli[1], poli[2]
                    linecoll_elem.append([0,poli[1], poli[2]])
                    
            if len(linecoll_elem) > 1 :
                line_collection.append(linecoll_elem)
                if plotyes == True:
                    point_coll.append([p1, polo[0]])
            
        
        
        """for c in c_points:
            #print c
            #exit(0)
            #ax.scatter(c[0], c[1], marker = '.' , c = 'black')
            line = [p1, c]
            path = shapely.geometry.LineString(line)
            
            linecoll_elem = []
            
            for poli in pol_inner:
            #    print "here!"
                intersection = path.intersection(poli[0])
                #print intersection.type
                if intersection.type == 'LineString' and len(linecoll_elem) == 0:   # not two of the inner
                    #plt.plot([p1[0], c[0]], [p1[1], c[1]], marker = '_')
                    #print "int line", intersection
                    #print poli[1], poli[2]
                    linecoll_elem.append([0,poli[1], poli[2]])
                    
            if len(linecoll_elem) == 0:
                continue
   
            for polo in pol_outer:
            #    print "here!"
                intersection = path.intersection(polo[0])
                #print intersection.type
                if intersection.type == 'LineString' and len(linecoll_elem) == 1:  
                    #plt.plot([p1[0], c[0]], [p1[1], c[1]], marker = '_')
                    #print "int line", intersection
                    #print polo[1], polo[2]
                    linecoll_elem.append([1,polo[1], polo[2]])
                 
            #linecoll_elem = np.array(linecoll_elem)
            if len(linecoll_elem) == 2:
                line_collection.append(linecoll_elem)
                if plotyes == True:
                    point_coll.append([p1, c]) # """
      
    # remove duplicates
    #line_collection = [s for s in line_collection if len(s[0]) in [3] ]
    
    
    print "end loops..."
    
    
    line_collection.sort()
    line_collection = list(line_collection for line_collection,_ in itertools.groupby(line_collection))
   

    #line_collection = np.array(line_collection)

    #print line_collection
    #print line_collection.shape
    
    
    #for l in line_collection:
    #    print "straight line", l 
        
        
    return line_collection, point_coll


    #######################################################################################################################################################
    #######################################################################################################################################################
    
#########################################################################################################################################################################
def line_param(t,p): # parameter form of a straight line 
    x = p[0] + p[3]*t
    y = p[1] + p[4]*t
    z = p[2] + p[5]*t
    return x,y,z    
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

def get_xy_of_line_angles(z,p):
    
    t = (z - p[1])/math.cos(p[2])*math.sin(p[3])
    
    y =        math.sin(p[2])*math.sin(p[3])*t    
    z = p[1] + math.cos(p[3])*t
        
    return y,z
#########################################################################################################################################################################
def get_xy_of_line_new(z,p): # not in use atm   
    #x = p[0] + p[1]*z
    #y = p[1] + p[3]*z    
    y = (z - p[2])/p[3]
    x = p[0] + p[1]*y
        
    return x,y
######################################################################################################################################################################### 
def get_xy_of_line(z,p): # not in use atm
    t = (z - p[2]) / p[5]    
    x = p[0] + p[3]*t
    y = p[1] + p[4]*t
        
    return x,y
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


def angle_3Dtracks(v1, v2):
# v1 is your firsr vector
# v2 is your second vector
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

    angle2 =  2 * np.pi - angle

    return min(angle,angle2)

#########################################################################################################################################################################
#########################################################################################################################################################################

def  dist_lines_volume(lines_fit_):
    print "-------------------------- distance to points in volume around center --------------------------" 
    # we are checking in a volume around the BGO, the distance to the lines of the points in the volume

    """fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection='3d')    
    ax_3d.set_xlim(-300, 300)
    ax_3d.set_ylim(-300, 300)   
    ax_3d.set_zlim(-300, 300) """
    
    # create points of volume:
    
    xbound_min = -100
    xbound_max =  100
    xwidth = 5

    ybound_min = -100
    ybound_max =  100
    ywidth = 5
    
    zbound_min = -200
    zbound_max =  200
    zwidth = 10       
    
    x_p = np.arange(xbound_min, xbound_max, xwidth)
    y_p = np.arange(ybound_min, ybound_max, ywidth)
    z_p = np.arange(zbound_min, zbound_max, zwidth)    
    xyz_grid = np.vstack((ndmesh(x_p,y_p,z_p))).reshape(3,-1).T       
    #ax_3d.scatter(xyz_grid[:,0], xyz_grid[:,1], xyz_grid[:,2])    
    #plt.show()

    nr_l = len(lines_fit_)
    
    dist_mat = np.zeros(len(xyz_grid))
    
    for m, p in enumerate(xyz_grid):
    
        distsum = 0
        for n, params in enumerate(lines_fit_):    
            dist2 = dist_p_line(p[0],p[1],p[2],params)  
            #print dist2  
            distsum = distsum + dist2
  
        dist_mat[m] = distsum/(1.0*nr_l)
        #print dist_mat[m]
    
    index_min = np.argmin(dist_mat)
    
    dist_mini = dist_mat[index_min]
    
    min_point = xyz_grid[index_min]
    
    
    #grid_at_min_xy = xyz_grid[xyz_grid[:,2]==min_point[2]].reshape(200/10,200/10,3)
    #grid_at_min_xz = xyz_grid[xyz_grid[:,1]==min_point[1]].reshape(200/10,200/10,3)
    #grid_at_min_yz = xyz_grid[xyz_grid[:,0]==min_point[0]].reshape(200/10,200/10,3)
    
    
    dist_mat = dist_mat.reshape(2*xbound_max/xwidth,2*ybound_max/ywidth,2*ybound_max/ywidth)
    
    indexmin2 = list(np.unravel_index(dist_mat.argmin(), dist_mat.shape))
    
    print "-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    
    print "MIN POINT (red)", min_point
    
    #print indexmin2[0], indexmin2[1], indexmin2[2]

    print "-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"

    #xyz_grid = xyz_grid.reshape(20,20,20,3)
    
    
    # TODO
    # maybe try this ndimage labelling function (also in vertex_reconstruction class for bgo) to find clusters
    # find min, find max and then set a cut to get clusters of lower distances... then draw and see what happens
    
    
    
    
    
    #print dist_mat.shape
    
    #print dist_mat#[:,1].shape
    
    
    #print dist_mat[:,1]
    """from matplotlib.cm import gnuplot
    fig2, ax2 = plt.subplots(1,3, figsize = (10,3))
    
    extentxy = [xbound_min,xbound_max,ybound_min,ybound_max]
    extentyz = [zbound_min,zbound_max,ybound_min,ybound_max]
    extentxz = [zbound_min,zbound_max, xbound_min,xbound_max]
    
    ax2[0].scatter(min_point[1], min_point[0])
    ax2[1].scatter(min_point[2], min_point[0])
    ax2[2].scatter(min_point[2], min_point[1])
    
    ax2[1].set_title("z vs x")
    ax2[2].set_title("z vs y")
    ax2[0].set_title("y vs x")
    
    ax2[0].imshow(dist_mat[:,:,indexmin2[2]], interpolation='none', extent = extentxy, aspect='auto', cmap = gnuplot, origin = "lower")
    ax2[2].imshow(dist_mat[indexmin2[0],:,:], interpolation='none', extent = extentyz, aspect='auto', cmap = gnuplot, origin = "lower")
    ax2[1].imshow(dist_mat[:,indexmin2[1],:], interpolation='none', extent = extentxz, aspect='auto', cmap = gnuplot, origin = "lower") # """
    
    #plt.imshow(dist_mat[:,1],interpolation='none')
        
    return min_point

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



