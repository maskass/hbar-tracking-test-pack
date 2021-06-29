#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.cm import coolwarm, gray_r



class BGOPortPainter:
    numOfChannels = 64
    pkgSize = 52.0 # mm
    def __init__(self, port,colourMap):
        self.preparePlot()

        if port == 'A':
            self.coordinates += [-self.pkgSize/2,+self.pkgSize/2]
        if port == 'B':
            self.coordinates += [-self.pkgSize/2,-self.pkgSize/2]
        if port == 'C':       
            self.coordinates += [+self.pkgSize/2,+self.pkgSize/2]
        if port == 'D':
            self.coordinates += [+self.pkgSize/2,-self.pkgSize/2]

        #print port
        self.patches = [self.dispatchRect(i) for i in range(self.numOfChannels)]
        self.collection = PatchCollection(self.patches, cmap=colourMap)

    def get_rects_for_tracking(hitarr):
        patches = self.patches[hitarr]
        return rects

    def plot(self, data, maxData):
        #print maxData
        self.collection.set_clim(0,maxData)
        self.collection.set_array(data.reshape(self.numOfChannels,))
        return self.collection

    def preparePlot(self):
        # define primitives
        self.corner = Rectangle(xy=[0,0], width=6.26, height=6.26, angle=0)
        self.updown = Rectangle(xy=[0,0], width=6.08, height=6.26, angle=0)
        self.side   = Rectangle(xy=[0,0], width=6.26, height=6.08, angle=0)
        self.middle = Rectangle(xy=[0,0], width=6.26, height=6.08, angle=0)

        # create coordinate table
        self.coordinates = np.zeros(self.numOfChannels*2).reshape(self.numOfChannels,2)
        
        # fancy indices
        self.corners = np.zeros(64, dtype=bool)
        self.corners[0]  = True
        self.corners[7]  = True
        self.corners[56] = True
        self.corners[63] = True
        self.corners = self.corners.reshape(8,8)
        self.sides = np.ones((8,8), dtype=bool)
        self.sides = self.sides - self.corners
        self.topbot= np.copy(self.sides)
        self.topbot[1:7,0:8] = False
        self.sides[0:8,1:7] = False
        self.middle = np.zeros((8,8), dtype=bool)
        self.middle[1:7,1:7] = True
        self.corners = self.corners.reshape(64,1)
        self.topbot = self.topbot.reshape(64,1)
        self.sides = self.sides.reshape(64,1)
        self.middle = self.middle.reshape(64,1)

        # coordinates for placement
        coordinateHelp = self.coordinates.reshape(8,8,2)
        # corner
        coordinateHelp[0,0] = [-3*6.08 - 6.26, -3*6.08 - 6.26]
        coordinateHelp[7,0] = [+3*6.08       , -3*6.08 - 6.26]
        coordinateHelp[0,7] = [-3*6.08 - 6.26, +3*6.08       ]
        coordinateHelp[7,7] = [+3*6.08       , +3*6.08       ]
        # sides
        coordinateHelp[0,1:7] = [[-3*6.08 - 6.26, i*6.08]for i in range(-3,3)]
        coordinateHelp[7,1:7] = [[+3*6.08       , i*6.08]for i in range(-3,3)]
        coordinateHelp[1:7,0] = [[i*6.08, -3*6.08 - 6.26]for i in range(-3,3)]
        coordinateHelp[1:7,7] = [[i*6.08, +3*6.08       ]for i in range(-3,3)]
        # middle
        coordinateHelp[1:7,1:7] = [ [[j*6.08, i*6.08]for i in range(-3,3)] for j in range(-3,3)]

    def dispatchRect(self,i):
        rect = None
        if self.corners[i] == True:
            rect = Rectangle(xy=self.coordinates[i], width=6.26, height=6.26, angle=0)
        if self.sides[i] == True:
            rect = Rectangle(xy=self.coordinates[i], width=6.08, height=6.26, angle=0)
        if self.topbot[i] == True:
            rect = Rectangle(xy=self.coordinates[i], width=6.26, height=6.08, angle=0)
        if self.middle[i] == True:
            rect = Rectangle(xy=self.coordinates[i], width=6.08, height=6.08, angle=0)
        return rect
        
######################################################################    

class BGOPort:
    numOfChannels = 64
    pkgSize = 52.0 # mm
    def __init__(self, event, port, calibData, pedestal):

        # comment: averaging is already done in the function that reads the calibration data, see below
        Qce = 0.500095/256.0
        signalscale_glob = 2.48972/1000.0
        
        #signalscale_glob = 0.00318569
        #Qce = 2.00012/256.0 
        
        
        self.signalscale = signalscale_glob

        pedestal = 10.0 # minus ten in nagatas code
        #pedestal = 0.0
    
        self.rawdata = np.array(event['BGODataPort{0}'.format(port)]) # rawdata of port A/B/C/D
        
        #print self.rawdata
        
        self.rawdata -= pedestal #*np.ones(64)
        self.rawdata /= calibData
       
        datatmp = np.copy(self.rawdata.reshape(8,8))

        if port == 'A':
            self.data = datatmp[::-1]
        if port == 'B':
            self.data = datatmp[::-1].T[::-1]
        if port == 'C':
            self.data = datatmp.T
        if port == 'D':
            self.data = datatmp.T[::-1].T
 
        #self.data = self.rawdata/self.signalscale
        plttmp = []

        #self.data *= 0.00290034  # convert to MeV
        #self.data *= 0.003
        self.data *= self.signalscale # convert to MeV
        self.data -= Qce
        

        for i in range(8):
            for j in range(8):
                plttmp.append(self.data[i,j])
        self.plotdata = np.array(plttmp)

    def getRawCharge(self):
        return np.sum(self.rawdata)

    def getCharge(self):
        return np.sum(self.data)

    def getMaxCharge(self):
        return np.max(self.data)

####################################################################
class BGOPainter:
    ports = ['A','B','C','D']
    def __init__(self,colourMap):
        self.portPainter = []
        for i,j in zip(self.ports, range(len(self.ports)) ):
            self.portPainter.append(BGOPortPainter(i,colourMap))

    def plot(self, data):
        maxData = data.getMaxCharge()
        return [i.plot(data.portData[idx].plotdata, maxData) for idx,i in enumerate(self.portPainter)]

class BGO:
    ports = ['A','B','C','D']
    pkgSize = 52.0 # mm
    def __init__(self, event, calibrationData):
        self.portData = []           
        for i,j in zip(self.ports, range(len(self.ports)) ):
            self.portData.append(BGOPort(event, i, calibrationData[j*64:j*64+64], 10) )
            
        #self.portData[0].preparePlot()
        self.coord   = np.zeros(256*2).reshape(16,16,2)
        self.dataMap = np.zeros(256  ).reshape(16,16,)
        self.prepareCords()
        self.dataMap[0: 8,8:16] = self.portData[0].data
        self.dataMap[0: 8,0: 8] = self.portData[1].data
        self.dataMap[8:16,8:16] = self.portData[2].data
        self.dataMap[8:16,0: 8] = self.portData[3].data

        
    def getMaxCharge(self):
        return np.max([i.getMaxCharge() for i in self.portData])

    def getRawCharge(self):
        retSum = 0
        for i in self.portData:
            retSum += i.getRawCharge()
        return retSum

    def getCharge(self):
        retSum = 0
        for i in self.portData:
            retSum += i.getCharge()
        return retSum

    def setData(self, dataMap):
        self.portData[0].data     = dataMap[0: 8,8:16]
        self.portData[0].plotdata = dataMap[0: 8,8:16]
        self.portData[1].data     = dataMap[0: 8,0: 8]
        self.portData[1].plotdata = dataMap[0: 8,0: 8]
        self.portData[2].data     = dataMap[8:16,8:16]
        self.portData[2].plotdata = dataMap[8:16,8:16]
        self.portData[3].data     = dataMap[8:16,0: 8]
        self.portData[3].plotdata = dataMap[8:16,0: 8]

    def prepareCords(self):
        # coordinates
        coordinates = np.zeros(64*2).reshape(64,2)
        coordinateHelp = coordinates.reshape(8,8,2)
        # corner
        coordinateHelp[0,0] = [-3*6.08 - 6.26, -3*6.08 - 6.26]
        coordinateHelp[7,0] = [+3*6.08       , -3*6.08 - 6.26]
        coordinateHelp[0,7] = [-3*6.08 - 6.26, +3*6.08       ]
        coordinateHelp[7,7] = [+3*6.08       , +3*6.08       ]
        # sides
        coordinateHelp[0,1:7] = [[-3*6.08 - 6.26, i*6.08]for i in range(-3,3)]
        coordinateHelp[7,1:7] = [[+3*6.08       , i*6.08]for i in range(-3,3)]
        coordinateHelp[1:7,0] = [[i*6.08, -3*6.08 - 6.26]for i in range(-3,3)]
        coordinateHelp[1:7,7] = [[i*6.08, +3*6.08       ]for i in range(-3,3)]
        # middle
        coordinateHelp[1:7,1:7] = [ [[j*6.08, i*6.08]for i in range(-3,3)] for j in range(-3,3)]
        coordA = coordinateHelp + [-self.pkgSize/2,+self.pkgSize/2]
        coordB = coordinateHelp + [-self.pkgSize/2,-self.pkgSize/2]
        coordC = coordinateHelp + [+self.pkgSize/2,+self.pkgSize/2]
        coordD = coordinateHelp + [+self.pkgSize/2,-self.pkgSize/2]

        self.coord[0: 8,8:16] = coordA
        self.coord[0: 8,0: 8] = coordB
        self.coord[8:16,8:16] = coordC
        self.coord[8:16,0: 8] = coordD


def readBGOCalibrationData(fileName):
    retVal = np.zeros(256)
    f = open(fileName,"r")
    for line in f:
        line = line.replace(","," ")
        vals = line.split()
        retVal[int(vals[0])] = float(vals[1])

    retVal = retVal / ( np.sum(retVal)/len(retVal) )
    return retVal

####################################################################
def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im.get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

if __name__ == "__main__":
    from rootpy.tree import TreeChain
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import logging
    from scalar import *
    # Most verbose log level
    logging.basicConfig(level=logging.DEBUG)

    from sys import argv
    from hodoscope import *

    if len(argv) < 2:
        print "usage: ", argv[0], "rootfile1.root rootfile2.root ... rootfileN.root  "
        exit(-1)

    bgoCalibData = readBGOCalibrationData("../calibrationData/gainBGOcalibration.dat")
    correctionData = readTimingCorrectionData("../calibrationData/hodoscopeTimingCorrectionTable.dat")  

    SclT = TreeChain("ScalarDataTree", argv[1:])

    startmixtime = 0.0
    endmixtime = 0.0

    for event in SclT:
        s = scalar(event)
        if s.mixingStart>0.0:
              startmixtime = s.midasTime
        if s.mixingStop>0.0:
              endmixtime = s.midasTime

    print "start and end timestamp of mixing: ", startmixtime, endmixtime, " Cusp run number: ", s.cuspRunNumber
   
    print 'start=', startmixtime
    print 'end=', endmixtime
    count = 0
    T = TreeChain("HbarEventTree", argv[1:])
    for event in T:  
        ts = event["midasTimeStamp"][0]
        #if ts < 1417923065:
        #    continue
        b = BGO(event, bgoCalibData)
        h = hodoscope(event)
        nI, nO, I, O=h.getActiveBar()
        print nI, nO
        
        if nI<2:
            continue
        #if nI>10 or nO>10:
        #    continue
        #if b.getCharge() < 30:
        #if b.getCharge() < 200:
        #    continue
        if b.getCharge() < 20:
            continue
        #print nI, nO
        if ts < startmixtime or ts > endmixtime:
            continue
        inner,outer,BG = h.drawHodoscope()
        datag= b.coord.reshape(256,2)
        data= b.dataMap.reshape(256,)
        bgoPainter = BGOPainter(coolwarm)
        bgoHandle=b
        a,b,c,d = bgoPainter.plot(b)
        fig, ax = plt.subplots(frameon=False)
        fig.patch.set_alpha(0.0)
        ax.set_xlim(-200,200)
        ax.set_ylim(-200,200)
        ax.add_collection(a    )
        ax.add_collection(b    )
        ax.add_collection(c    )
        ax.add_collection(d    )
        ax.add_collection(BG   )
        ax.add_collection(inner)
        ax.add_collection(outer)
        #hodoColourBar = fig.colorbar(inner)
        #hodoColourBar.set_label("Hodoscope energy deposit")
        #BGOColourBar = fig.colorbar(a)
        #BGOColourBar.set_label("BGO energy deposit")
        ax.set_aspect(1)
        ax.set_ylabel("vertical position [mm]")
        ax.set_xlabel("horizontal position [mm]")
        ax.patch.set_facecolor='w'
        #print ax.get_images()
        plt.tight_layout()
        ni,no,nd,nf=h.getActiveBar()
        
        ax.text(-70, -145, unicode('BGO Edep: {0:.2f}MeV'.format(bgoHandle.getCharge()), 'latin-1'))
        ax.text(-70, +145, unicode('Time: {0}s'.format(ts-startmixtime), 'latin-1'))

        print ts
        count = count+1
        print "number of event: ", count
        plt.savefig("./event_plots/event_{0}.pdf".format(count), facecolor=fig.get_facecolor())
        plt.savefig("./event_plots/event_{0}.svg".format(count), facecolor=fig.get_facecolor())
        #plt.show()
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax = fig.gca(projection='3d')
        #select = data>-1
        #ax.scatter(datag[select,1],datag[select,0],data[select])
        #plt.show()
        #print b.getRawCharge()

