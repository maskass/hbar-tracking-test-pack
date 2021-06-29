#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import math
import matplotlib.backends.backend_pdf
from matplotlib.lines import Line2D             
from matplotlib.cm import  gray_r, gray, bone
from matplotlib import lines as mpl_lines 
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
import itertools as it 
from numpy import random

from itertools import groupby
from operator import itemgetter
import itertools
import matplotlib.pyplot as plt
import re

import scipy.stats as stats
#import cv2

import sys
sys.path.append("..")

from src.trackplot_final import *
from src.trackfinding_new_workon_sims_plots_for_thesis import *
from src.hodoscope_trackhelper import *
from src.trackfitting2D_sims import *
#from src.trackfinding_fitting_3D_final import *
from src.vertex_reconstruction_2D import *
from src.hodoscope_tracking import * # same as just hodoscope with a few extra functions for tracking
from src.bgo import * # newest calibration
from src.scalar import *

from numpy import ndarray

cosmic_count = 0
event_count_tot = 0

##############################################################################################################     
##############################################################################################################
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
            c2 = 0
        else:
            c2 = c + 1
                
        if c == 0:
            c0 = 31
        else:
            c0 = c - 1
                                       
        array_clustercheck = [c0,c,c2]
                
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
            c2 = 0
        else:
            c2 = c + 1
                    
        if c == 0:
            c0 = 31
        else:
            c0 = c - 1
                                          
        array_clustercheck = [c0,c,c2]
                    
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
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
if __name__ == "__main__":
    import ROOT
    from rootpy.tree import Tree, TreeChain
    import logging
    # Most verbose log level
    logging.basicConfig(level=logging.DEBUG)

    #from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    from sys import argv

    if len(argv) < 2:
        print "usage: ", argv[0], "file1.root  "
        exit(-1)


    # these are just read in so the same python classes as for the measured data can be used:
    # (below dummy hodoscope and bgo are created)
    bgoCalibData = readBGOCalibrationData("../calibrationData/gainBGOcalibration.dat")
    correctionData = readTimingCorrectionData("../calibrationData/hodoscopeTimingCorrectionTable.dat")  
    #corrtiming = np.loadtxt("./calibrationData/cable_length_corrections.dat")
    cable_corr_inner = np.loadtxt('../calibrationData/cable_corrections_inner_up_down_2_15mev.dat')
    cable_corr_outer = np.loadtxt('../calibrationData/cable_corrections_outer_up_down_2_15mev.dat')
    thresi = np.loadtxt("../calibrationData/inner_charge_cuts_up_down_pbar_ang.dat")
    threso = np.loadtxt("../calibrationData/outer_charge_cuts_up_down_pbar_ang.dat")
    #threso = np.loadtxt("outer_amp_cuts_up_down_cosmic.dat")
    fibre_inner_cuts = np.loadtxt('../calibrationData/inner_fibre_tot_cuts_2019.dat')
    fibre_outer_cuts = np.loadtxt('../calibrationData/outer_fibre_tot_cuts_2019.dat')   
    
    
    midas_run_nr = 999
    
    cand_list = [] # list of candidates, used for filling the root tree

    ####################################################################################################################################################################################
    
    print argv[1]
    
    print " "
    
    # the tree of the read in file:
    T = TreeChain("HbarEventTree",argv[1])  
      
    # dummy tree:  
    Tdummy = TreeChain("HbarEventTree", "./rootfiles/output-run-6351_for_dummy_event.root")
    for n, dummy in enumerate(Tdummy):
        if n > 0:
            continue
        dummy_bgo   = BGO(dummy, bgoCalibData)    
        dummy_hodor = hodoscope(dummy, thresi, threso, fibre_inner_cuts, fibre_outer_cuts) # """
    ####################################################################################################################################################################################
    ####################################################################################################################################################################################

    plotyes = False  # creates event plots
    plotshow = False        # show the event plots during analysis
    
    writeTxt = True
   
    simplecuts = False

    BGOcut = 0.0
    
    fibrebool = True
    
    BGOon = False

    adcTimestamp = 0.0
    mixtime = []

    numbers = np.arange(32)
    
    cuspRunNumber = 0
    startmixAdc = None
    mtd = 999
    #########################################################################################################################################################################
    #########################################################################################################################################################################
    if writeTxt == True:
        
        file_mtds = open("tracking_results_" + argv[1][40:-5] + ".dat", 'w')
        #file_mtds = open('tests.dat', 'w')
        #file_dists = open('dists.dat', 'w')
        
        #file_mtds = open('test_pbars_0_0_smaller_errors_eps150.dat', 'w')
        #file_dists = open('sims_pbars_0_0_smaller_error_fibres_dists.dat', 'w')  
        
        """file_amps_o = open('ftf_boc.dat', 'w')    
        file_amps_i = open('ftf_bic.dat', 'w')    
        #file_ofc = open('ftf_ofc.dat', 'w')    
        #file_ifc = open('ftf_ifc.dat', 'w') # """  
        
        
        """file_oc = open('ftfbic_oc.dat', 'w')    
        file_ic = open('ftfbic_ic.dat', 'w')    
        file_ofc = open('ftfbic_ofc.dat', 'w')    
        file_ifc = open('ftfbic_ifc.dat', 'w')"""
        """file_oc = open('ftfpbert_oc.dat', 'w')    
        file_ic = open('ftfpbert_ic.dat', 'w')    
        file_ofc = open('ftfpbert_ofc.dat', 'w')    
        file_ifc = open('ftfpbert_ifc.dat', 'w')        
        #file_mtds = open('test.dat', 'w')     """
    #########################################################################################################################################################################
    
    tot_events = 0
    #########################################################################################################################################################################
    #########################################################################################################################################################################    
    #########################################################################################################################################################################
    #########################################################################################################################################################################
    active = np.zeros(32)

    # mapping of hodoscope channels of the simulations to the real numbering
    # NOT needed anymore, they already come out of the simulations in the right order  
    #sim_real_map = [11,10,9,8,7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12]

    #file_amps_i = open('./sims_amplitudes/chips_bar_amps_inner.dat', 'w')
    #file_amps_o = open('./sims_amplitudes/chips_bar_amps_outer.dat', 'w')
    
    for evtnumber, event in enumerate(T): 
    
       
        #if BGOon == True:
        bgo_E = event['BGOEnergy']*1.0
        BGOYPos = event['BGOXPos']*1000.0
        BGOXPos = event['BGOYPos']*1000.0
        BGOZPos = event['BGOZPos']*-1000.0
        #else:
        #    bgo_E = -999
        #    BGOXPos = 999
        #    BGOYPos = 999
        #    BGOZPos = 999
            
        ####################################
        ####################################
        # first we read in the data from the current event:
        ####################################
        ####################################
        
        zerrI = 59# measured z error of bars, inner layer (FWHM), in mm
      
        zerrO = 73# measured z error of bars, outer layer (FWHM), in mm
        
        inner_z_sigma = zerrI/2.355   
        outer_z_sigma = zerrO/2.355
        
        InnerCharge = np.array(event['ChargeInnerLayer'])#[sim_real_map]
        InnerPosY = np.array(event['YPositionInnerLayer'])*1000.0#[sim_real_map]*1000.0 
        InnerPosX = np.array(event['XPositionInnerLayer'])*1000.0#[sim_real_map]*1000.0
        InnerPosZ = np.array(event['ZPositionInnerLayer'])*-1000.0#[sim_real_map]*-1000.0 
        
        for n in range(0, len(InnerPosZ)):
            InnerPosZ[n] = InnerPosZ[n] + np.random.normal(0.0, inner_z_sigma, 1 )
               
        OuterCharge = np.array(event['ChargeOuterLayer'])#[sim_real_map]
        OuterPosY = np.array(event['YPositionOuterLayer'])*1000.0#[sim_real_map]*1000.0
        OuterPosX = np.array(event['XPositionOuterLayer'])*1000.0#[sim_real_map]*1000.0
        OuterPosZ = np.array(event['ZPositionOuterLayer'])*-1000.0#[sim_real_map]*-1000.0   
        
        for n in range(0, len(OuterPosZ)):
            OuterPosZ[n] = OuterPosZ[n] + np.random.normal(0.0, outer_z_sigma, 1 )

        InnerFibreCharge = np.array(event['ChargeInnerFibre'])
        InnerFibrePosY = np.array(event['YPositionInnerFibre'])*1000.0
        InnerFibrePosX = np.array(event['XPositionInnerFibre'])*1000.0
        InnerFibrePosZ = np.array(event['ZPositionInnerFibre'])*-1000.0  
        
        #print "-----------------", len(InnerFibreCharge)
        
        OuterFibreCharge = np.array(event['ChargeOuterFibre'])
        OuterFibrePosY = np.array(event['YPositionOuterFibre'])*1000.0
        OuterFibrePosX = np.array(event['XPositionOuterFibre'])*1000.0
        OuterFibrePosZ = np.array(event['ZPositionOuterFibre'])*-1000.0  
        
        
        # writing out the deposited energy in the bars and fibres
              
        """for fi_ in InnerFibreCharge:
            file_ifc.write("%s " % fi_)
        file_ifc.write("\n")
        file_ifc.flush()  
        
        for fo_ in OuterFibreCharge:
            file_oc.write("%s " % fo_)
        file_ofc.write("\n")
        file_ofc.flush()  # """  


        """for fi_ in InnerCharge:
            file_ic.write("%s " % fi_)
        file_ic.write("\n")
        file_ic.flush()  
        
        for fo_ in OuterCharge:
            file_oc.write("%s " % fo_)
        file_oc.write("\n")
        file_oc.flush()  # """          
     
        #continue
        
        sims_fibre_inner_hits = np.array([InnerFibrePosX, InnerFibrePosY, InnerFibrePosZ])
        sims_fibre_outer_hits = np.array([OuterFibrePosX, OuterFibrePosY, OuterFibrePosZ])
        sims_BGO_hit = np.array([BGOXPos, BGOYPos, BGOZPos])
  
        for n,i in enumerate(InnerCharge):
            if i < 0.7:
                InnerCharge[n] = 0.0
                
        for n,i in enumerate(OuterCharge):
            if i < 0.7:
                OuterCharge[n] = 0.0  #"""

        print "-----------------", len(InnerFibreCharge)        
                   
        I = InnerCharge[:]>0.0
        nI = np.sum(I)
        O = OuterCharge[:]>0.0
        nO = np.sum(O)
        
        
        
        if nI < 1:
            continue

        if nO < 1:
            continue

        if nI > 16:
            continue

        if nO > 16:
            continue    
            
        #print I
        #print O
        #print nI, nO
        #print np.logical_and(I,O)
        #if np.sum(np.logical_and(I,O)) == 0:# == True
        #    continue           
        
        #if bgo_E < 0.5:
        #    continue
        
        for n,i in enumerate(InnerFibreCharge):
            if i < 0.53:
                InnerFibreCharge[n] = 0.0
               
        for n,i in enumerate(OuterFibreCharge):
            if i < 0.53:
                OuterFibreCharge[n] = 0.0  #"""


        fI = InnerFibreCharge > 0
        fO = OuterFibreCharge > 0
        
        
        print "-----------------", len(InnerFibreCharge)
        
       
        currI = numbers[I]
        currO = numbers[O]
        

        print "##########################################################################################"
        print "########################              NEW EVENT!             #############################"
        print "##########################################################################################"
        print "CUSP, Midas, Event:", cuspRunNumber, midas_run_nr, evtnumber
        event_count_tot = event_count_tot + 1


        print "BGO ENERGY (MeV): ", bgo_E
        print "I:", nI, currI
        print "O:", nO, currO
                
        print "currI : ", currI
        print "currO : ", currO
                

        #currI_uncut = currI
        #currO_uncut = currO
        #####################################################################################################################
        #####################################################################################################################
               
        currI_o, currO_o = is_it_an_orphan(currI, currO) # orphans excluded

        print "currI for tracking (after orphan removing and cuts):", currI_o
        print "currO for tracking (after orphan removing and cuts):", currO_o
        
        number_of_orphans_I = len(currI) - len(currI_o)
        number_of_orphans_O = len(currO) - len(currO_o)
        
        print "number of orphan hits inner: ", number_of_orphans_I
        print "number of orphan hits outer: ", number_of_orphans_O
        
  
        if nI > 0 or nO > 0:
                          
            cornerplotI, cornerplotO = getCornerMap()
            pmapI, pmapO = prepare_PosMap(I,O,nI, nO)
            
            cmap = prepare_CornerMap(I,O, nI, nO)
            #####################################################################################################################
            #####################################################################################################################         
            ynI, clusI = find_consecutive_nr(currI)
            ynO, clusO = find_consecutive_nr(currO)
                        
            clusI = [s for s in clusI if len(s) > 1]
            clusO = [s for s in clusO if len(s) > 1]
            
            #if len(clusI) < 1:
            #    continue
            
            print "Inner cluster: ", clusI
            print "Outer cluster: ", clusO
            #####################################################################################################################
            #####################################################################################################################  """       

        #exit(0)
        print " "
        #####################################################################################################################
        #####################################################################################################################
        
        nr_of_cluster = 999
        cluster_dist = 999
        polygon_hitlist = []
        vertex = [0.0,0.0]
        vertex2 = [0.0,0.0]
        
        print "VEEEEEEEEEERTEX", vertex2
        

        #####################################################################################################################
        #####################################################################################################################

        print "###########################         FIND TRACKS:        ##################################"    

        ######################################################################################
        ######################################################################################  
        plot_trackfinding = [True,True] # points, lines   # save rotation lines and their rotation points (in bgo) for plotting  
        trackcollection = 999
        linepoints = 999
        bgopoints = 999
        
        #vertex_BGO = [BGOXPos, BGOYPos, BGOZPos ]

        if BGOon != True:
            tpx_length = 28.0 # mm      
            tpx_hit_bounds = [[-0.5*tpx_length,0.5*tpx_length],[-0.5*tpx_length,0.5*tpx_length]]
            ext = [(0.5*tpx_length, 0.5*tpx_length),  (-0.5*tpx_length, 0.5*tpx_length), (-0.5*tpx_length, -0.5*tpx_length), (0.5*tpx_length, -0.5*tpx_length)]            
            tpx_polygon = Polygon(ext)  
            
            bounds_for_2D_vertex = [-0.5*tpx_length, 0.5*tpx_length] 
               
            vertex2 = [0.0,0.0]
            vertex_err = [tpx_length,tpx_length]
        
        
            if nI > 0 or nO > 0:
                trackcollection, linepoints, bgopoints, union_cl_I, union_cl_O = track_finding_polygons_tpx(pmapI, cmap, currI, currO, tpx_polygon,
                                                         plot_trackfinding, tpx_hit_bounds, InnerPosZ, OuterPosZ) #, xm, ym)"""
                                                         
                                                         
                if trackcollection != 999:
                    print "number of 2D lines:", len(trackcollection)
                else: 
                    print "number of 2D lines:", 0     
               
        else:
            if nI > 0 or nO > 0:
                hit_circle = Circle(xy=(BGOXPos, BGOYPos), radius=20.0,  color = 'white', lw=2, alpha = 0.5)
                
                trackcollection, linepoints, bgopoints, union_cl_I, union_cl_O =  track_finding_circle(pmapI, cmap, currI, currO, hit_circle,BGOXPos, BGOYPos,
                                                         plot_trackfinding, InnerPosZ, OuterPosZ) #, xm, ym)
                

        
        

        nr_of_fit_lines_ = 999
        mean_d_pt_ = [999,999,999]
        av_d_pt_ = [999,999,999]
        n_cluster_ = 999
        biggest_cluster_ = [999,999,999,999,999]
        vertex_fitted_ = [999,999,999]
        vf_covar_ = [[999,999,999],[999,999,999],[999,999,999]] 
        cov_l_cl_ = [[999,999,999],[999,999,999],[999,999,999]] 
        tot_cov_ =[[999,999,999],[999,999,999],[999,999,999]] 
        hess_bool_ = 999
        lines_in_l_cl_ = 999
        
        #print "OOOOOOOOOOI", nI, nO
        if fibrebool == True and nI > 0 and nO > 0:
            #innerzPos = -innerzPos
            #outerzPos = -outerzPos
            
            #print "OI!"
                           

           
            #print InnerFibreCharge
            circle = Circle(xy=(0, 0), radius=45.0,  color = 'white', lw=2, alpha = 0.5)
            ax, mps__, nr_of_fit_lines_, mean_d_pt_, av_d_pt_, n_cluster_, biggest_cluster_, dists_, vertex_fitted_, vf_covar_, cov_l_cl_, tot_cov_, av_dist_cl_, hess_bool_, lines_in_l_cl_, mean_clus_xyz_, total_covar_mean_,  covar_largest_cl_mn_ = track_finding_4layers(polygon_hitlist, dummy_hodor,
                                                                                         nO, nI, I, O,
                                                                                         InnerFibreCharge, OuterFibreCharge,
                                                                                         zerrI, zerrO, cuspRunNumber, evtnumber,
                                                                                         circle, union_cl_I, union_cl_O,
                                                                                         sims_fibre_inner_hits, sims_fibre_outer_hits, sims_BGO_hit,
                                                                                         InnerPosX,InnerPosY,InnerPosZ,
                                                                                         OuterPosX,OuterPosY,OuterPosZ)
                                                                                         
         
         
            numbers_fi = np.arange(63) # number array of the inner fibre channels
            numbers_fo = np.arange(100) # number array of the outer fibre channels 
            
            
            #fO = OuterFibreCharge > 0 # active fibres outer layer (boolian array)
            #fI = InnerFibreCharge > 0 # active finbres inner layer (boolian array)
            
            curr_fI = numbers_fi[fI] # channels numbers of hits inner fibres
            curr_fO = numbers_fo[fO] # channel numbers of hits outer fibres
            
            
            """for ix, iy, iz in zip(InnerPosX[I],InnerPosY[I],InnerPosZn[I]):
                        ax[2].plot(ix, iy, 'o', markerfacecolor = 'black',  markersize=3)
                        ax[1].plot(iz, ix, 'o', markerfacecolor = 'black',  markersize=3)    
                        ax[0].plot(iz, iy, 'o', markerfacecolor = 'black',  markersize=3)             

            for ox, oy, oz in zip(OuterPosX[O],OuterPosY[O],OuterPosZn[O]):
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
                ax[0].plot(sims_if[2], sims_if[1], 'o', markerfacecolor = 'black',  markersize=3)
                
                plt.show()   #"""
                
                #av_all_dists = np.mean(dists_)
                # write out the distances between fround 3D tracks
            """if type(mps__).__name__ != 'int':
            
                print "------------------------------------------------------------------------------------------------------------------------------------------", mps__.shape, dists_.shape
                print mps__
                print dists_
         
         
                for di_, dis_ in zip(mps__, dists_):
                    file_dists.write("%s " % evtnumber)
                    file_dists.write("%s " % bgo_E)
                    file_dists.write("%s " % nI)
                    file_dists.write("%s " % nO)
                    #file_dists.write("%s " % biggest_cluster_[0])
                    #file_dists.write("%s " % av_dist_cl_)
                    #file_dists.write("%s " % av_all_dists)
                    file_dists.write("%s " % di_[0]) 
                    file_dists.write("%s " % di_[1]) 
                    file_dists.write("%s " % di_[2]) 
                    file_dists.write("%s " % dis_) 
                    
                    file_dists.write("\n")
                file_dists.flush()  # """                                                                        
            
            print "BIGGEST CLUSTER", biggest_cluster_
            print "BIGGEST CLUSTER", biggest_cluster_[0]
            print "BIGGEST CLUSTER", biggest_cluster_[1]


        print "Final Hitcollection:      ################################################################" 
        
        if trackcollection != 999:
            nr_of_tracks = len(trackcollection)
        
            if nr_of_tracks > 0:
                for i in trackcollection:
                    print i
                       
                print " "
        else:
            nr_of_tracks = 0
        print "##############################        Fitting 2D:       ##################################"
        print "Fitting", nr_of_tracks, "tracks "
        line_params = []
        #angle1 = 999
        angle = 999
        vertex = [999,999]
               
        vertex_err = [50,50]
        bounds_for_2D_vertex = [-0.5*50, 0.5*50] 
        
        #angle = do_le_anglecalc(trackcollection, I, O, xm, ym, chargeI, chargeO, bgo_E)
        if nr_of_tracks > 1:
            line_params, fitlines_points, trackscoll = do_le_2Dfitting_sims(trackcollection, I, O, sims_BGO_hit, vertex_err, InnerCharge, OuterCharge) 
        
        print "#######################        VERTEX 2D through fitting:       ##########################"
        
        
        angle1 = 999
        angle2 = 999
        angle3 = 999
        angleY = 999
        mean_angleY = 999 
        orientatio = [999,999,999]
        mtd3 = 999
        mtd_min_max = 999
        fitlines_points = 999
        line_paras = 999
        vertex_points = 999
        vertex_2d_lines = [999,999]
                    
        #vertex_points = 999
        """if nr_of_tracks > 1:
            vertex_points, vertex_2d_lines = determine_vertex(line_params, bounds_for_2D_vertex) # vertex is the arithm mean of vertex points

        if math.isnan(vertex[0]) or math.isnan(vertex[1]):
            vertex_2d_lines = [999,999]
        
        #mtds2, angle2 = calc_mtds_of_tracks_corr(trackcollection, innerCFs, outerCFs, innerzPos, outerzPos, vertex2) #, rinner, router)
        
        
        
        print "Vertex: ", vertex_2d_lines"""
        
    
        
        
        if nr_of_tracks == 1:
            angleY, orientatio = calc_angle_for_1track(trackcollection, vertex2)
            
            
        #orientatio = None
        """if nr_of_tracks > 1:
            mtd3, angle_all, angleY, mtd_min_max, orientatio = calc_mtds_of_tracks_corr_std(trackcollection, innerCFs, outerCFs, innerzPos, outerzPos, vertex2) #, rinner, router)
            print "all", angle_all    
            angle1 = np.amax(angle_all)
            print "ang1", angle1
            
            angle_all = np.setdiff1d(angle_all,angle1)
            print angle_all
            angle2 = 0
            angle3 = 0        
            
            if nr_of_tracks  == 2:
                angle = angle1
            
            if nr_of_tracks > 2:
                angle2 = np.amax(angle_all)
                angle_all = np.setdiff1d(angle_all,angle2)
                if len(angle_all) > 1:
                    angle3 = np.amax(angle_all)
                else:
                    angle3 = angle2
                    angle2 = angle1
                    #angle3 = angle1
                #angle_all = np.setdiff1d(angle_all,angle3)
                    
        mean_angleY = np.mean(np.array(angleY))    """           
        
        if simplecuts == True:
            if nr_of_tracks < 2:
                continue
                
            if bgo_E < BGOcut:
                continue
                              
            if nr_of_tracks == 2 and angle1 > 160:
                continue
                    
            if mtd > 0.4:
                continue 
       
        print "------------------------------------------------------------------------------------------"
        print "ANGLE: ", angle1, angle2, angle3, "degree" 
        print "Y ANGLE: ", angleY, "degree, mean: ", mean_angleY
        
        print "------------------------------------------------------------------------------------------"
        
        if orientatio:
            print "Orientation: ", orientatio

        #pval_c, pval_p = calc_pvals_mtds(mtd_combi)
            
        """if plotyes == True: # and timestampx < 0.0219999 and timestampx > 0.021000:

                bgopoints = np.array(bgopoints)
                plot_all_bars_2D(circle, union_cl_I, union_cl_O, linepoints, trackcollection, line_params, vertex_BGO) 
                
                plt.show()  """  
        
        if writeTxt == True:
        
            #for nb in noise_bars_O:
        
                    #file_mtds.write("%s " % nb)
                    #file_mtds.write("%s " % h.o.CFU[nb])
                    #file_mtds.write("%s " % h.o.CFD[nb])
                    #file_mtds.write("%s " % h.o.AmpU[nb])
                    #file_mtds.write("%s " % h.o.AmpD[nb])
                    file_mtds.write("%s " % cuspRunNumber)
                    file_mtds.write("%s " % midas_run_nr)
                    file_mtds.write("%s " % evtnumber)
                    file_mtds.write("%s " % bgo_E)    
                    file_mtds.write("%s " % nI)    
                    file_mtds.write("%s " % nO)                    
                    file_mtds.write("%s " % nr_of_tracks)  
                    file_mtds.write("%s " % angle1)
                    file_mtds.write("%s " % angle2)
                    file_mtds.write("%s " % angle3)
                    file_mtds.write("%s " % mean_angleY)
                    file_mtds.write("%s " % orientatio[0])
                    file_mtds.write("%s " % orientatio[1])
                    file_mtds.write("%s " % orientatio[2])               
                    file_mtds.write("%s " % tot_events) # 14
          
                    file_mtds.write("%s " % mtd3)
                    file_mtds.write("%s " % mtd_min_max)
                    
                    file_mtds.write("%s " % vertex_2d_lines[0])
                    file_mtds.write("%s " % vertex_2d_lines[1])

                    file_mtds.write("%s " % vertex2[0])
                    file_mtds.write("%s " % vertex2[1])
                    
                    
                    file_mtds.write("%s " % nr_of_cluster) # 21
                    file_mtds.write("%s " % cluster_dist)
                    file_mtds.write("%s " % 999)
                    file_mtds.write("%s " % 999)
                    file_mtds.write("%s " % number_of_orphans_I)
                    file_mtds.write("%s " % number_of_orphans_O)
                    
                    file_mtds.write("%s " % nr_of_fit_lines_) # 27
                    
                    file_mtds.write("%s " % mean_d_pt_[0]) # 28
                    file_mtds.write("%s " % mean_d_pt_[1])
                    file_mtds.write("%s " % mean_d_pt_[2]) 
                    
                    file_mtds.write("%s " % av_d_pt_[0]) # 31
                    file_mtds.write("%s " % av_d_pt_[1])
                    file_mtds.write("%s " % av_d_pt_[2]) 
                    
                    file_mtds.write("%s " % n_cluster_) 
                    file_mtds.write("%s " % biggest_cluster_[2])  # 35
                    file_mtds.write("%s " % biggest_cluster_[3])
                    file_mtds.write("%s " % biggest_cluster_[4]) 
                    
                    file_mtds.write("%s " % biggest_cluster_[0]) # 38
                    file_mtds.write("%s " % biggest_cluster_[1]) 

                    file_mtds.write("%s " % BGOXPos) # 40
                    file_mtds.write("%s " % BGOYPos)  
                    file_mtds.write("%s " % BGOZPos)                   
                                                     
                    file_mtds.write("%s " % vertex_fitted_[0]) # 43
                    file_mtds.write("%s " % vertex_fitted_[1])  
                    file_mtds.write("%s " % vertex_fitted_[2]) 

                    file_mtds.write("%s " % vf_covar_[0][0]) # 46
                    file_mtds.write("%s " % vf_covar_[1][1])  
                    file_mtds.write("%s " % vf_covar_[2][2])  
                    
                    file_mtds.write("%s " % cov_l_cl_[0][0]) # 49
                    file_mtds.write("%s " % cov_l_cl_[1][1])  
                    file_mtds.write("%s " % cov_l_cl_[2][2]) 
                    
                    file_mtds.write("%s " % tot_cov_[0][0]) # 52
                    file_mtds.write("%s " % tot_cov_[1][1])  
                    file_mtds.write("%s " % tot_cov_[2][2]) # 54
                    file_mtds.write("%s " % hess_bool_)  
                    file_mtds.write("%s " % lines_in_l_cl_   )     #56  
                    
                    file_mtds.write("%s " % mean_clus_xyz_[0]) # 57
                    file_mtds.write("%s " % mean_clus_xyz_[1])  
                    file_mtds.write("%s " % mean_clus_xyz_[2]) 
                    
                    file_mtds.write("%s " % total_covar_mean_[0][0]) # 60
                    file_mtds.write("%s " % total_covar_mean_[1][1])  
                    file_mtds.write("%s " % total_covar_mean_[2][2]) # 62
                    
                    file_mtds.write("%s " % covar_largest_cl_mn_[0][0]) 
                    file_mtds.write("%s " % covar_largest_cl_mn_[1][1])  
                    file_mtds.write("%s " % covar_largest_cl_mn_[2][2]) 
                                                                               
                                                             
                    if writeTxt == True:    
                            file_mtds.write("%s " % 999) 
                            file_mtds.write("\n")
                            file_mtds.flush()   
            
                
        """cand_list.append([cuspRunNumber, evtnumber, bgo_E, nI_all, nO_all, nr_of_tracks, timestampx, 0.0, mtd3, mtd_min_max, angle1, angle2, angle3,
        vertex[0], vertex[1], vertex2[0], vertex2[1], nr_of_cluster, cluster_dist, mean_angleY, orientatio[0], orientatio[1], orientatio[2], number_of_noise_hits_I,
        number_of_noise_hits_O, number_of_orphans_I, number_of_orphans_O])"""
              
        print "------------------------------------------------------------------------------------------"    
        print "Total events in run:", tot_events #, timestampx
        print "------------------------------------------------------------------------------------------"
            
    ########################################################################################################################################################################
    ############################################################################################################################################################################ 
    ############################################################################################################################################################################
    ############################################################################################################################################################################ 
