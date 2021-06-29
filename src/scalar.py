#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

class scalar:
    def __init__(self, event):
        self.signalTrigger = event['ScalarData'][ 0]
        self.notInUse      = event['ScalarData'][ 1]
        self.hodoscopeRate = event['ScalarData'][ 2]
        self.ADPulse       = event['ScalarData'][ 3]
        self.PIPSRate      = event['ScalarData'][ 4]
        self.mixingStart   = event['ScalarData'][ 5]
        self.mixingStop    = event['ScalarData'][ 6]
        self.dsScinti      = event['ScalarData'][ 7]
        self.cuspDAQClock  = event['ScalarData'][ 8]
        self.bresciaRate   = event['ScalarData'][ 9]
        self.normScinti1   = event['ScalarData'][10]
        self.normScinti2   = event['ScalarData'][11]
        self.normScinti3   = event['ScalarData'][12]
        self.normScinti4   = event['ScalarData'][13]
        self.normalisation = event['ScalarData'][14]

        self.midasTime     = event['midasTimeStamp'][0]
        self.cuspRunNumber = event['CUSPRunNumber' ][0]
	self.timeStamp     = event['Timestamp'][0]

def plotPulse(ax, data, idx, colour):
    for i in data[data[:,idx]>0,0]:
        ax.axvline(i,c=colour,lw=1,alpha=0.5)

def annotateCuspRunNumbers(ax,number,idx,data,dataIDX):
    for i,j in zip(number,idx):
        ax.annotate(s="CUSPRunNumber={0}".format(number),xy=(data[idx,0][0],data[idx,dataIDX]+1000), xycoords='data')#,xytext=(0, 0), textcoords='offset points', arrowprops=dict(arrowstyle="->"))


def findMixingCycles(data):
    startIDX = data[:,2]>0 
    stopIDX  = data[:,3]>0 

    startTime = data[startIDX,0]
    mixing = []
    for start in startTime:
        idx = data[:,0]>start
        if np.sum(stopIDX & idx) == 0:
            continue
        stop = data[(stopIDX & idx),0][0]
        #print start, stop, stop-start
        if stop-start > 100 and stop-start<=251:
            mixing.append([data[data[:,0]==start,1][0], start, stop])

    return np.array(mixing)

def correctDuplicateMixings(mixing):
    if len(mixing) == 0:
        return 0
    #print mixing
    retval = []
    mix, idx = np.unique(mixing[:,0],return_index=True)
    for i,j in zip(mix,idx):
        select = mixing[:,0] == i
        retval.append(mixing[select,:])
        #print mixing[select,2]-mixing[select,1]
        #print np.sum(select), i
        
    return np.array(retval)

if __name__ == "__main__":
    import ROOT
    ROOT.PyConfig.StartGuiThread = False
    from rootpy.tree import TreeChain
    import matplotlib.pyplot as plt
    
    fig,ax = plt.subplots(nrows=3,ncols=1, sharex=True,frameon=False)
    fig.subplots_adjust(hspace=0.1,wspace=0)
    
    import logging
    # Most verbose log level
    logging.basicConfig(level=logging.DEBUG)

    from sys import argv

    if len(argv) < 2:
        print "usage: ", argv[0], "rootfile1.root rootfile2.root ... rootfileN.root  "
        exit(-1)

    T = TreeChain("ScalarDataTree", argv[1:])
    
    data = []
    for event in T:
        scl = scalar(event)
        data.append([scl.midasTime, scl.cuspRunNumber, scl.mixingStart, scl.mixingStop, scl.ADPulse, scl.bresciaRate, scl.hodoscopeRate, scl.signalTrigger])
        del(scl)
    del(T)
    
    data = np.array(data)#,dtype='float64')
    ts, idx = np.unique(data[:,0],return_index=True)
    sumdata = []
    for i,j in zip(ts,idx):
        print i, ts[-1]
        select = (data[:,0]==i)
        #data[select,0] = data[select,0]+np.linspace(i,i+1,np.sum(select))
        tmp = [i,data[j,1]]
        tmp.extend(np.sum(data[select,2:],axis=0).tolist()[:])
        sumdata.append(tmp)
    
    
    sumdata = np.array(sumdata)
    
    CuspNumbers, cuspIDX = np.unique(data[:,1],return_index=True)

    mixing = findMixingCycles(sumdata)

    #mixing = correctDuplicateMixings(mixing)
    with open("mixing.dat","w") as f:
        f.write("0 {0} -1\n".format(len(ts)))
        for i in mixing:
            f.write("{0} {1} {2}\n".format(i[0], i[1], i[2]))
    print mixing
    #exit(0)

    ax[0].plot(sumdata[:,0],sumdata[:,5])
    ax[0].set_title('Brescia detector rate')
    ax[0].set_ylabel("Counts [Hz]")
    #annotateCuspRunNumbers(ax[0],CuspNumbers, cuspIDX, sumdata, 5)
    print CuspNumbers
    plotPulse(ax[0], sumdata, 4, 'g')
    plotPulse(ax[0], sumdata, 2, 'r')
    plotPulse(ax[0], sumdata, 3, 'c')

    ax[1].plot(sumdata[:,0],sumdata[:,6])
    ax[1].set_title('Hodoscope only rate')
    ax[1].set_ylabel("Counts [Hz]")
    plotPulse(ax[1], sumdata, 4, 'g')
    plotPulse(ax[1], sumdata, 2, 'r')
    plotPulse(ax[1], sumdata, 3, 'c')

    ax[2].plot(sumdata[:,0],sumdata[:,7])
    ax[2].set_title('BGO+Hodoscope trigger rate')
    ax[2].set_ylabel("Counts [Hz]")
    ax[2].set_xlabel("Time [s]")
    plotPulse(ax[2], sumdata, 4, 'g')
    plotPulse(ax[2], sumdata, 2, 'r')
    plotPulse(ax[2], sumdata, 3, 'c')

    fig.subplots_adjust(hspace=0.2,wspace=0)

    #plt.show()
    plt.savefig("test.pdf")
