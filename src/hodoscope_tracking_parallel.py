#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.cm import coolwarm,bone,gray
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import colorConverter
from itertools import *
from operator import *

#TDC0_inner = [58,50,42,34,8,16,24,32,1,9,17,25,60,52,44,36,59,51,43,35,61,53,45,37]
#TDC0_outer = [98,90,82,74,66,58,22,30,46,7,15,23,31,39,50,97,89,81,73,65,57,49,99,91,83,75,67,59,51]

#TDC1_inner = [57,49,41,33,56,48,40,6,14,22,30,7,15,23,31,62,54,46,38,63,55,39,47,28,20,12,4,29,21,13,5,27,19,11,3,26,18,10,2]
#TDC1_outer = [38,47,8,16,24,32,40,48,1,9,17,25,33,41,2,10,19,27,11,3,42,34,26,18,44,36,28,20,12,4,43,35,100,92,84,76,68,60,52,53,54,61,62,69,70,71,77,78,79,85,86,87,93,94,95,5,13,21,29,37,45,6,14,63,55,96,88,80,72,64,56] 

###############################################################################################################################
class fibre_layer:
    def __init__(self, lf_list, pos, threshs): # pos can be "Inner" or "Outer"
    
    
        LE_UPPER_CUT = -60.0
        LE_LOWER_CUT = -100.0
        LE_CUT_bool = False
        
        TOT_CUT = 0        
        
       
        if pos == "Inner":
            self.nOfFibres = 63
            self.diameter = 167.0 # mm, inner diameter
            self.chWidth = 4.0 
            self.borderUp = 127.0 # mm, outer border of most outer fibre bundle, seen from center of hodoscope
            self.borderDown = 125.0 # mm 4x63 = 252, 127+125 = 252
            
            self.LE = np.zeros(64*128)
            self.ToT = np.zeros(64*128)
            for b in range(0, 64*128):      
                self.LE[b] = lf_list.TdcLEInner.GetValue(b)
                
                self.ToT[b] = lf_list.TdcWidthInner.GetValue(b)
                              
            # rows are the stamps and columns the channels
            self.LE = self.LE.reshape(128,64)
            self.ToT = self.ToT.reshape(128,64)
                       
            self.LE_first = (self.LE)[0,:] 
            self.ToT_first = (self.ToT)[0,:]   
               
            
            """if np.sum(self.ToT_first[38]) != 0:
                
                print "tot", self.LE
                print "tot shape", (self.LE).shape
              

                print "shape2", np.max(self.ToT, axis = 0)
                
                
                print "lol ?", self.LE[np.argmax(self.ToT, axis = 0)][0]
                
                #print "len1", len(self.ToT[np.argmax(self.ToT[:,:]),:])
                #print "len2", len(np.argmax(self.ToT, axis = 0))
                
                print " "
                
                #print np.argmax(self.ToT[:,:])
                
                exit(0)"""
            
            """if self.LE_first[1] != 0:
                print self.LE_first[1], self.ToT_first[1]       """

            # I added azero in the beginning of this array, cos we start at channel 1!
            self.which_TDC = np.array([0,       0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1])
            #self.TDC0_inner = [ 1, 8, 9, 16, 17, 24, 25, 32, 34, 35, 36, 37, 42, 43, 44, 45, 50, 51, 52, 53, 58, 59, 60, 61] 
            #self.TDC1_inner = [ 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 30, 31, 33, 38, 39, 40, 41, 46, 47, 48, 49, 54, 55, 56, 57, 62, 63] # 
            
            trigger = np.array([lf_list.TdcTriggerTime.GetValue(0),lf_list.TdcTriggerTime.GetValue(1)])
            
            trigger_sub_i = self.which_TDC
            
            
            #######################################################################################
            self.ToT_largest = self.ToT[np.argmax(self.ToT, axis = 0)][0]
            self.LE_largestToT = self.LE[np.argmax(self.ToT, axis = 0)][0]  
            
            self.ToT_largest = (np.delete(self.ToT_largest, 0))
            self.ToT_largest = self.ToT_largest[:63]
            
            les_il = (self.LE_largestToT - trigger_sub_i)  # substract trigger time
            les_il[self.LE_largestToT==0] = 0   # set those entries zero where there was no hit
           
            les_il = np.delete(les_il, 0)
            
            les_il = les_il[:63]
            self.LE_largestToT = les_il
            #######################################################################################
                
            trigger_sub_i[self.which_TDC  == 0] = trigger[0]  # set trigger times in the array
            trigger_sub_i[self.which_TDC  == 1] = trigger[1]
       
            les_i = (self.LE_first - trigger_sub_i)  # substract trigger time
            les_i[self.LE_first==0] = 0   # set those entries zero where there was no hit

            les_i = np.delete(les_i, 0)
            self.ToT_first = np.delete(self.ToT_first, 0)
                        
            if LE_CUT_bool == True:
                le_mask_i = np.logical_and(les_i > LE_UPPER_CUT, les_i < LE_LOWER_CUT)
                les_i[le_mask_i] = -999
                self.ToT_first[le_mask_i] = -999     
                                     
            les_i = les_i[:63]
            
            self.ToT_first = self.ToT_first[:63]
        
            self.LEsub = les_i
            
            """if self.LE_first[1] != 0:
                #print  "trig", trigger_sub_i[1]
                #print "sub", self.LEsub[0], self.ToT_first[0]   """    
                
            #print "leeeeeeeeeeeeeeeeeeen tot i", len(self.ToT_first), len(threshs) 
                      
            self.active = np.logical_and(self.LEsub != -999, np.greater(self.ToT_first, threshs))

           
        elif pos == "Outer":
            self.nOfFibres = 100
            self.diameter = 293.0 # mm, inner diameter
            self.chWidth = 4.0 
            self.borderUp = 201.75 # mm, outer border of most outer fibre bundle, seen from center of hodoscope
            self.borderDown = 198.25 # mm 4x100 = 400, 201.75+198.25 = 400
            
            self.LE = np.zeros(128*128)
            self.ToT = np.zeros(128*128)
            for b in range(0, 128*128):      
                self.LE[b] = lf_list.TdcLEOuter.GetValue(b)
                #print lf_list.TdcLEOuter.GetValue(b)
                self.ToT[b] = lf_list.TdcWidthOuter.GetValue(b)

            self.LE = self.LE.reshape(128,128) 
            self.ToT = self.ToT.reshape(128,128)
            self.LE_first = (self.LE)[0,:]  
            self.ToT_first = (self.ToT)[0,:]  
 
            # I added azero in the beginning of this array, cos we start at channel 1!and I added 27 0 (zeros) in the end to match the length of 128
            self.which_TDC = np.array([0,     1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1,    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

            
            trigger = np.array([lf_list.TdcTriggerTime.GetValue(0),lf_list.TdcTriggerTime.GetValue(1)])
            
            trigger_sub_o = self.which_TDC
            
            #######################################################################################
            self.ToT_largest = self.ToT[np.argmax(self.ToT, axis = 0)][0]
            self.LE_largestToT = self.LE[np.argmax(self.ToT, axis = 0)][0]  
            
            self.ToT_largest = (np.delete(self.ToT_largest, 0))
            self.ToT_largest = self.ToT_largest[:100]
            
            les_ol = (self.LE_largestToT - trigger_sub_o)  # substract trigger time
            les_ol[self.LE_largestToT==0] = 0   # set those entries zero where there was no hit
           
            les_ol = np.delete(les_ol, 0)
            
            les_ol = les_ol[:100]
            self.LE_largestToT = les_ol
            #######################################################################################
                
            trigger_sub_o[self.which_TDC  == 0] = trigger[0]  # set trigger times in the array
            trigger_sub_o[self.which_TDC  == 1] = trigger[1]
       
            les_o = (self.LE_first - trigger_sub_o)  # substract trigger time
            les_o[self.LE_first==0] = 0   # set those entries zero where there was no hit
            
            self.LEsub = les_o

            les_o = np.delete(les_o, 0)     
            self.ToT_first = (np.delete(self.ToT_first, 0))
            
            if LE_CUT_bool == True:
                le_mask_o = np.logical_and(les_o > LE_UPPER_CUT, les_o < LE_LOWER_CUT)
                les_o[le_mask_o] = -999
                self.ToT_first[le_mask_o] = -999  
                                
            les_o = les_o[:100]
            
            #print "leeeeeeeeeeeeeeeeeeeeen", len(les_o)
            self.ToT_first = self.ToT_first[:100]
                         
            self.LEsub = les_o
            
            #print "leeeeeeeeeeeeeeeeeeen tot o", len(self.ToT_first), len(threshs) 
                                  
            self.active = np.logical_and(self.LEsub != -999, np.greater(self.ToT_first, threshs))
            #self.active = self.ToT_first > 0
            #print "active", self.active
            #self.LEsub[self.active] = 
            
            #self.TDC0_outer = [ 7, 15, 22, 23, 30, 31, 39, 46, 49, 50, 51, 57, 58, 59, 65, 66, 67, 73, 74, 75, 81, 82, 83, 89, 90, 91, 97, 98, 99]
            #self.TDC1_outer =[1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,24,25,26,27,28,29,32,33,34,35,36,37,38,40,41,42,43 ,44,45,47,48,52 ,53,54,55,56,60,61,62, 63,64,68,69,70,71,72,76,77,78,79,80,84,85,86,87,88,92,93,94,95,96, 100]
            


###############################################################################################################################
    def getPositionMap(self):
        PosZ = []
        center_first_fibre = self.borderUp - self.chWidth*0.5
        for f in range(self.nOfFibres):
            PosZ.append(center_first_fibre - f*self.chWidth) # this gives the center of the fibre bundle as z pos
          
        return np.array(PosZ)
###############################################################################################################################        
###############################################################################################################################        
###############################################################################################################################
class layer:
    def __init__(self, lf_list, pos, thres_layer):
        assert type(pos) is str, "position has to be a string! but got {0}".format(type(pos))

        self.nOfBars = 32
        self.chargeU = np.zeros(self.nOfBars)
        self.chargeD = np.zeros(self.nOfBars)
        self.AmpU    = np.zeros(self.nOfBars)
        self.AmpD    = np.zeros(self.nOfBars)
        self.LEU     = np.zeros(self.nOfBars)
        self.LED     = np.zeros(self.nOfBars)
        self.CFU     = np.zeros(self.nOfBars)
        self.CFD     = np.zeros(self.nOfBars)
        
        self.active = np.zeros(self.nOfBars)
        
        for b in range(0, self.nOfBars):
        
            if pos == 'Inner':
                    self.chargeU[b] = lf_list.ChargeInnerLayerUpstream.GetValue(b)
                    self.chargeD[b] = lf_list.ChargeInnerLayerDownstream.GetValue(b)
                    self.AmpU[b] = lf_list.AmplitudeInnerLayerUpstream.GetValue(b)
                    self.AmpD[b] = lf_list.AmplitudeInnerLayerDownstream.GetValue(b)
                    self.LEU[b] = lf_list.LEtimeStampsInnerLayerUpstream.GetValue(b)
                    self.LED[b] = lf_list.LEtimeStampsInnerLayerDownstream.GetValue(b)
                    self.CFU[b] = lf_list.CFtimeStampsInnerLayerUpstream.GetValue(b)
                    self.CFD[b] = lf_list.CFtimeStampsInnerLayerDownstream.GetValue(b)  
                    
                                     
            elif pos == 'Outer':            
                    self.chargeU[b] = lf_list.ChargeOuterLayerUpstream.GetValue(b)
                    self.chargeD[b] = lf_list.ChargeOuterLayerDownstream.GetValue(b)
                    self.AmpU[b] = lf_list.AmplitudeOuterLayerUpstream.GetValue(b)
                    self.AmpD[b] = lf_list.AmplitudeOuterLayerDownstream.GetValue(b)
                    self.LEU[b] = lf_list.LEtimeStampsOuterLayerUpstream.GetValue(b)
                    self.LED[b] = lf_list.LEtimeStampsOuterLayerDownstream.GetValue(b)
                    self.CFU[b] = lf_list.CFtimeStampsOuterLayerUpstream.GetValue(b)
                    self.CFD[b] = lf_list.CFtimeStampsOuterLayerDownstream.GetValue(b)           
        

        self.activeLE = np.logical_and(self.LEU > 0, self.LED > 0)
        self.activeCF = np.logical_and(self.CFU > 0, self.CFD > 0)
        #self.active  = self.activeCF 
        #self.active  = self.activeLE
        self.active  = self.activeAmplitude(pos, thres_layer)
        #self.trigger = self.AmpTrigger(trigger_thres)
###############################################################################################################################
    def activeAmplitude(self, c, thresh):
        active = np.zeros(self.nOfBars,dtype=np.bool)


        for i in range(self.nOfBars):
            #print i
            #print thresh[i]
            res = True
            
            if (self.chargeU[i]<=thresh[i][0] or self.chargeD[i]<=thresh[i][1]):
            #if (self.AmpU[i]<=0.2 or self.AmpD[i]<=0.2):
                        res = False    
                
            active[i]=res    
        


        ######################

        return active

####################################################################################################################################

    def findNbOfActiveBarsWithClusters(self):
        hits = []
        if len(self.active) == 0:
            return [len(hits),hits]

        cluster= []
        
        for i in range(len(self.active)-1):
            cluster.append(self.active[i])
            if (self.active[i+1]-self.active[i]) == 1:
                pass
            else:
                hits.append(cluster)
                cluster = []

        # treat last element
        cluster.append(self.active[-1])

        if cluster[-1]==31 and self.active[0]==0:
            if len(hits)==0:
                hits.append(cluster)
            else:
                hits[0][len(hits[0]):] = cluster
        else:
            hits.append(cluster)
        
        return [len(hits),hits]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
class hodoscope:
    def __init__(self, lf_list, thresi, threso, threshi_f, thresho_f):
        self.i = layer(lf_list, "Inner", thresi)
        self.o = layer(lf_list, "Outer", threso)
        self.fi = fibre_layer(lf_list, "Inner", threshi_f)
        self.fo = fibre_layer(lf_list, "Outer", thresho_f)

        self.BG = []


####################################################################################################################################

    def getActiveBar(self):
        return np.sum(self.i.active), np.sum(self.o.active), self.i.active, self.o.active

####################################################################################################################################        
    def setActiveBars(self, currI, currO):
        #print "-------------------------------------------- ", self.i.active 
        self.i.active = currI
        self.o.active = currO
####################################################################################################################################

    def getActiveCluster(self):
        return 0,0,0,0
        numI, clusterI = self.i.findNbOfActiveBarsWithClusters()
        numO, clusterO = self.o.findNbOfActiveBarsWithClusters()
        return numI, numO, clusterI, clusterO

####################################################################################################################################

    def getActiveClusterCount(self):
        numI, dump = self.i.findNbOfActiveBarsWithClusters()
        
        numO, dump = self.o.findNbOfActiveBarsWithClusters()
        return [numI, numO]

####################################################################################################################################

    def getCharge(self):
        retValInner = [
            np.arange(len(self.i.chargeU))[self.i.active],
            self.i.chargeU[self.i.active], 
            self.i.chargeD[self.i.active] ]
        retValOuter = [                 
            np.arange(len(self.o.chargeU))[self.o.active],
            self.o.chargeU[self.o.active], 
            self.o.chargeD[self.o.active]  ]
        return np.array(retValInner), np.array(retValOuter) ;

####################################################################################################################################

    def getAmplitude(self, corrTable = []):
        if len(corrTable) == 0:
            corrTable = np.zeros(len(self.i.CFU)*5).reshape(32,5)
            corrTable[:,3] = np.ones(32)
            corrTable[:,4] = np.ones(32)

        retValInner = [
            np.arange(len(self.i.AmpU))[self.i.active],
            self.i.AmpU[self.i.active], 
            self.i.AmpD[self.i.active]*corrTable[self.i.active,3] ]
        retValOuter = [                 
            np.arange(len(self.o.AmpU))[self.o.active],
            self.o.AmpU[self.o.active], 
            self.o.AmpD[self.o.active]*corrTable[self.o.active,4]  ]
        return np.array(retValInner), np.array(retValOuter) ;

####################################################################################################################################

    def getTotalAmplitude(self, corrTable = []):
        if len(corrTable) == 0:
            corrTable = np.zeros(len(self.i.CFU)*5).reshape(32,5)
            corrTable[:,3] = np.ones(32)
            corrTable[:,4] = np.ones(32)

        retVal = [ 
            self.i.AmpU, 
            self.i.AmpD * corrTable[:,3] , 
            self.o.AmpU, 
            self.o.AmpD * corrTable[:,4]  ]
        return np.array(retVal);

####################################################################################################################################

    def getTimingLE(self):
        retValInner = [
            np.arange(len(self.i.LEU))[self.i.active],
            self.i.LEU[self.i.active], 
            self.i.LED[self.i.active] ]
        retValOuter = [                 
            np.arange(len(self.o.LEU))[self.o.active],
            self.o.LEU[self.o.active], 
            self.o.LED[self.o.active]  ]
        return np.array(retValInner), np.array(retValOuter) ;

####################################################################################################################################

    def getTimingCF(self, corrTable = []):
        if len(corrTable) == 0:
            corrTable = np.zeros(len(self.i.CFU)*7).reshape(32,7)
            corrTable[:,3] = np.ones(32)
            corrTable[:,4] = np.ones(32)

        retValInner = [
            np.arange(len(self.i.CFU))[self.i.active],
            self.i.CFU[self.i.active]                             , 
            self.i.CFD[self.i.active] + corrTable[self.i.active,1]]
        retValOuter = [                 
            np.arange(len(self.o.CFU))[self.o.active],
            self.o.CFU[self.o.active]                             , 
            self.o.CFD[self.o.active] + corrTable[self.o.active,2]]
        return np.array(retValInner), np.array(retValOuter) 
 
#################################################################################################################################### 
        
    def getZPositionHelper(self, timingData, amplitudeData, params):
        dist              = params[0] # mm

        timingVariance    = params[1] # mm^2
        timingScaling     = params[2]
        timingOffset      = params[3] # mm
        amplitudeVariance = params[4] # mm^2
        amplitudeScaling  = params[5]
        amplitudeOffset   = params[6] # mm

        variance       = 1./(1./timingVariance + 1./amplitudeVariance)
  
        timeUpstream   = timingData   [0,:]
        timeDownstream = timingData   [1,:]
        amplUpstream   = amplitudeData[0,:]
        amplDownstream = amplitudeData[1,:]

        tDiff = timeUpstream - timeDownstream
        aComp = amplUpstream/(amplUpstream+amplDownstream)

        Zpos  = (tDiff-timingOffset   )*timingScaling   *dist/timingVariance
        Zpos -= (aComp-amplitudeOffset)*amplitudeScaling*dist/amplitudeVariance
        Zpos *= variance
 
        return Zpos

####################################################################################################################################

    def getZPosition(self, corrTable = []):
        paramsInner = [99.5,   73.7**2, 0.84,  0.  ,  61.1**2, 13.5, 0.5]
        paramsOuter = [172.5, 126.**2,  0.45, -0.02, 116.**2,   9.2, 0.5]

        cfI,  cfO  = self.getTimingCF (corrTable)
        ampI, ampO = self.getAmplitude(corrTable)
        retValInner = [
            np.arange(len(self.i.LEU))[self.i.active],
            self.getZPositionHelper(cfI[1:,:], ampI[1:,:], paramsInner)]
        retValOuter = [
            np.arange(len(self.o.LEU))[self.o.active],
            self.getZPositionHelper(cfO[1:,:], ampO[1:,:], paramsOuter)]

        return np.array(retValInner), np.array(retValOuter)

####################################################################################################################################

    def getDistO(self, p):
        if len(p) != 2:
            return 0
        PosMap = self.getPositionMap()
        d = PosMap[1][int(p[0])]-PosMap[1][int(p[1])]
        d = np.sqrt(np.sum(d*d))
        return d

####################################################################################################################################

    def getDistI(self, p):
        if len(p) != 2:
            return 0
        PosMap = self.getPositionMap()
        d = PosMap[0][int(p[0])]-PosMap[0][int(p[1])]
        d = np.sqrt(np.sum(d*d))
        return d

####################################################################################################################################
    def GetInnerAngle(self):
        if np.count_nonzero(self.i.active) > 2:
            if len(self.GetClusters())==2:
                hits = []
                print self.GetClusters()
                for i in range(0,2):   
                    hits.append(np.mean(self.GetClusters()[i]))   
                #print "clusters angle: ", abs(hits[0] - hits[1])
                if abs(hits[0] - hits[1]) in range(12,20):
                    return False

            else:
                numar = np.zeros(32)
                for i in range(0,32):
                    numar[i] = i
                AI = numar[self.i.active]     
                if abs(AI[0] - AI[1]) in range(12,20):
                    return False
        return True
####################################################################################################################################       
    def GetClusters(self):
        numar = np.zeros(32)
        for i in range(0,32):
            numar[i] = i
        AI = numar[self.i.active]
        c_count = 0
        cluster_list = []
        for k, g in groupby(enumerate(AI), lambda (i,x):i-x):
            c_count = c_count + 1
            cluster_list.append(map(itemgetter(1), g))
        return cluster_list
                  
####################################################################################################################################
    def rotate(self, b, ang):
        x,y = b
        return x*np.cos(ang)-y*np.sin(ang), x*np.sin(ang)+y*np.cos(ang)
        
####################################################################################################################################

    def getPositionMap(self):
        InnerCoordinatesStart = [(103-2.5 ,-40-0.375/2+10), (103-2.5,-20+10), (103-2.5, 0+0.375+10), (103-2.5, 20+0.375*2+10)]
        OuterCoordinatesStart = [(175-2.5,-70-0.375/2+17.5), (175-2.5,-35+17.5), (175-2.5,0+0.375+17.5), (175-2.5,35+0.375*2+17.5)]
        angles = [x*2*np.pi/8 for x in range(8)]
        InnerCoordinates = []
        OuterCoordinates = []

        for i,ang in enumerate(angles):
            for j in InnerCoordinatesStart:
                #print j 
                InnerCoordinates.append(self.rotate(j, ang))
            for j in OuterCoordinatesStart:
                OuterCoordinates.append(self.rotate(j, ang))

        return np.array(InnerCoordinates), np.array(OuterCoordinates)
        
####################################################################################################################################

    def getPositionMapTracking(self):

        idI = np.arange(0,21,2)
        idO = np.arange(0,36,5)

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
                    InnerCoordinates.append(self.rotate(currentcoord, ang))
        
                innerC.append(InnerCoordinates)
            

        for i,ang in enumerate(angles):
            for ind in range(0,len(OuterCoordinatesStart)):
                OuterCoordinates = []
                for k in idO:
                    currentcoord = [OuterCoordinatesStart[ind][0], OuterCoordinatesStart[ind][1]+k]
                    OuterCoordinates.append(self.rotate(currentcoord, ang))

                outerC.append(OuterCoordinates)
                
        return np.array(innerC), np.array(outerC)
        
####################################################################################################################################

    def getCornerMap(self):
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
                InnerCoordinates.append([self.rotate(j, ang),self.rotate(s, ang),self.rotate(k, ang),self.rotate(h, ang)])
            for j,s,k,h in zip(OuterCoordinatesLeftDown, OuterCoordinatesLeftUp,OuterCoordinatesRightDown, OuterCoordinatesRightUp):
                OuterCoordinates.append([self.rotate(j, ang),self.rotate(s, ang),self.rotate(k, ang),self.rotate(h, ang)])

        return np.array(InnerCoordinates), np.array(OuterCoordinates)


####################################################################################################################################

    def drawHodoscope(self, colourMap=coolwarm):
    
   
        # left down corner of upper panel
        InnerCoordinatesStart = [(103 ,-40-0.375/2), (103,-20), (103, 0+0.375), (103, 20+0.375*2)]
        OuterCoordinatesStart = [(175,-70-0.375/2), (175,-35), (175,0+0.375), (175,35+0.375*2)]
        angles = [x*2*np.pi/8 for x in range(8)]
        InnerCoordinates = []
        OuterCoordinates = []
        maxData = 0
        for i,ang in enumerate(angles):
            for j in InnerCoordinatesStart:
                InnerCoordinates.append(self.rotate(j, ang))
            for j in OuterCoordinatesStart:
                OuterCoordinates.append(self.rotate(j, ang))

        # create primitive boxes, inner layer
        patches = []
        values  = []
        bg = []
        for i,val in enumerate(self.i.active):
            if len(bg)<64:
                bg.append(Rectangle(xy=InnerCoordinates[i], width=20, height=5, angle=angles[i//4]*180/np.pi+90, alpha=0.03, color='g'))
            if val == True:
                patches.append(Rectangle(xy=InnerCoordinates[i], width=20, height=5, angle=angles[i//4]*180/np.pi+90))
                values.append(self.i.AmpU[i]+self.i.AmpD[i])
                
        innerCollection = PatchCollection(patches, cmap=colourMap, zorder = 100)
        innerCollection.set_array(np.array(values))
        if len(values)>0:
            maxdata = max(values)
        else: 
            maxdata = 0 

        # create primitive boxes, outer layer
        patches = []
        values  = []
        for i,val in enumerate(self.o.active):
            if len(bg)<64:
                bg.append(Rectangle(xy=OuterCoordinates[i], width=35, height=5, angle=angles[i//4]*180/np.pi+90, alpha=0.03, facecolor='r'))

            if val == True:
                patches.append(Rectangle(xy=OuterCoordinates[i], width=35, height=5, angle=angles[i//4]*180/np.pi+90))
                values.append(self.o.AmpU[i]+self.o.AmpD[i])
                #values.append(1)
        outerCollection = PatchCollection(patches, cmap=colourMap, zorder= 100)
        outerCollection.set_array(np.array(values))

        values.append(maxData)
        maxData = max(values)

        BG = PatchCollection(bg, cmap=bone, alpha = 0.1, zorder= 100)
        BG.set_array(np.ones(64))
        
        innerCollection.set_clim(0,maxData)
        outerCollection.set_clim(0,maxData)
        return innerCollection, outerCollection, BG

####################################################################################################################################

    def drawHodoscope3D(self, ax):
        innerCollection, outerCollection, BG = self.drawHodoscope()
        innerZlength = 200 #mm
        outerZLength = 350 #mm
        BG.set_alpha(0.5)
        ax.add_collection3d(poly, zs=zs, zdir='y')

####################################################################################################################################

    def Get_HitBoxes(self):
        # left down corner of upper panel
        InnerCoordinatesStart = [(103 ,-40-0.375/2), (103,-20), (103, 0+0.375), (103, 20+0.375*2)]
        OuterCoordinatesStart = [(175,-70-0.375/2), (175,-35), (175,0+0.375), (175,35+0.375*2)]
        angles = [x*2*np.pi/8 for x in range(8)]
        InnerCoordinates = []
        OuterCoordinates = []
        maxData = 0
        for i,ang in enumerate(angles):
            for j in InnerCoordinatesStart:
                InnerCoordinates.append(self.rotate(j, ang))
            for j in OuterCoordinatesStart:
                OuterCoordinates.append(self.rotate(j, ang))

        # create primitive boxes, inner layer
        patches = []
        values  = []
        bg1 = []
        bg2 = []
        for i,val in enumerate(self.i.active):
            if len(bg1)<64:
                bg1.append(Rectangle(xy=InnerCoordinates[i], width=20, height=5, angle=angles[i//4]*180/np.pi+90, alpha=0.03, color='black'))
            if val == True:
                patches.append(Rectangle(xy=InnerCoordinates[i], width=20, height=5, angle=angles[i//4]*180/np.pi+90))
                values.append(1)
                
        #innerCollection = PatchCollection(patches, cmap=colourMap)
        #innerCollection.set_array(np.array(values))
        if len(values)>0:
            maxdata = max(values)
        else: 
            maxdata = 0 

        # create primitive boxes, outer layer
        patches = []
        values  = []
        for i,val in enumerate(self.o.active):
            if len(bg2)<64:
                bg2.append(Rectangle(xy=OuterCoordinates[i], width=35, height=5, angle=angles[i//4]*180/np.pi+90, alpha=0.03, facecolor='black'))

            if val == True:
                patches.append(Rectangle(xy=OuterCoordinates[i], width=35, height=5, angle=angles[i//4]*180/np.pi+90))
                values.append(1)
                #values.append(1)
        #outerCollection = PatchCollection(patches, cmap=colourMap)
        #outerCollection.set_array(np.array(values))

        values.append(maxData)
        maxData = max(values)

        return bg1, bg2
       
####################################################################################################################################

def readTimingCorrectionData(filename):
    return np.loadtxt(filename)

if __name__ == "__main__":
    from rootpy.tree import TreeChain
    import matplotlib.pyplot as plt
    import logging
    # Most verbose log level
    logging.basicConfig(level=logging.DEBUG)

    from sys import argv

    if len(argv) < 2:
        print "usage: ", argv[0], "rootfile1.root rootfile2.root ... rootfileN.root  "
        exit(-1)

    T = TreeChain("HbarEventTree", argv[1:])
    for event in T:  
        h = hodoscope(event)
        I, O = h.getTimingCF()
        #nI, nO, I, O=h.getActiveCluster()
        inner,outer = h.drawHodoscope()
        #print I
        #print O
        #print nI, nO
        #print h.getPositionMap()
        fig, ax = plt.subplots()
        ax.set_xlim(-200,200)
        ax.set_ylim(-200,200)
        ax.add_collection(inner)
        ax.add_collection(outer)

        plt.show()

   
        """which_TDC = []
        for i in range(1,101):
            if i in self.TDC0_outer:
                which_TDC.append(0)
            elif i in self.TDC1_outer:
                which_TDC.append(1)

        print which_TDC
        print len(which_TDC)
        exit(0)"""

        """last_les = []
        for i in self.LE:
            last = i[i!=0]
            if len(last) > 0:
                last_les.append(last[-1])
            else:
                last_les.append(0)

        last_tots = []
        for j in self.ToT:
            last = j[j!=0]
            if len(last) > 0:
                last_tots.append(last[-1])
            else:
                 last_tots.append(0)

        self.LE_last = np.array(last_les)
        self.ToT_last = np.array(last_tots) """
        
        
        
        
                   
        """last_les = []
        for i in self.LE:
            last = i[i!=0]
            if len(last) > 0:
                last_les.append(last[-1])
            else:
                last_les.append(0)

        last_tots = []
        for j in self.ToT:
            last = j[j!=0]
            if len(last) > 0:
                last_tots.append(last[-1])
            else:
                last_tots.append(0
)
        self.LE_last = np.array(last_les)
        self.ToT_last = np.array(last_tots)"""
