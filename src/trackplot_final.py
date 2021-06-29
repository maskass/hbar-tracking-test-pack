#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
#import Image
import math
from matplotlib.lines import Line2D             
from matplotlib import lines as mpl_lines 
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle, Wedge


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as M3
from matplotlib import cm

from hodoscope_trackhelper import *
#import shapely
#from shapely.geometry import Polygon

from shapely.ops import cascaded_union
from descartes.patch import PolygonPatch
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
import scipy.stats
from src.bgo import *

from pyPdf import PdfFileWriter, PdfFileReader

def append_pdf(input,output):
    [output.addPage(input.getPage(page_num)) for page_num in range(input.numPages)]

##############################################################################################################
##############################################################################################################

def plot_event_2D_2(evtnumber, cuspRunNumber, h, bgo, linepoints, trackcollection, bgohitlist, bgopoints, coords_f,
line_params, fitlines_points, mtd, timebla, plot_trackfinding, vertex_points, vertex, vertex2):

        inner,outer,BG = h.drawHodoscope()
        
        datag= bgo.coord.reshape(256,2)
        data= bgo.dataMap.reshape(256,)
        bgoPainter = BGOPainter(coolwarm)
        bgoHandle=bgo
        #print bgoPainter.middle
        a,b,c,d = bgoPainter.plot(bgo)
        
        
        fig, ax = plt.subplots( frameon=False, figsize = (7,7))
        
        fig.patch.set_alpha(0.0)
        ax.set_xlim(-200,200)
        ax.set_ylim(-200,200)
       
        ax.add_collection(BG   )
        ax.add_collection(inner)
        ax.add_collection(outer)
        ax.add_collection(a    )
        ax.add_collection(b    )
        ax.add_collection(c    )
        ax.add_collection(d    )

        ax.set_aspect(1)    
        ax.set_ylabel("vertical position (mm)", fontsize = 15)
        ax.set_xlabel("horizontal position (mm)", fontsize = 15)
        ax.patch.set_facecolor='w'

      ################################################################################################
        # plot lines from the rotation
        if plot_trackfinding[1] == True:
            plot_tracklines2D(linepoints,ax[0])
        ################################################################################################
        
        ################################################################################################
        # plot white circle to indicate BGO 
        r = 50.0
        circ = Circle([0,0], r)
        p = PatchCollection([circ], alpha=1.0)
        p.set_edgecolor("w")
        p.set_facecolor("None")
        ax.add_collection(p   )
        ################################################################################################
        plt.tight_layout()
        #ni,no,nd,nf=h.getActiveBar()
        #print startmixtime - endmixtime       
        
        ax.text(-70, -135, unicode('BGO Edep: {0:.2f} MeV'.format(bgo.getCharge()), 'latin-1'), fontsize = 14)
        if mtd != 999:
            ax.text(-70, -150, unicode('Time diff: {0:.3f} ns'.format(mtd), 'latin-1'), fontsize = 14)
        ax.text(-70, +145, unicode('Time: {0:.5f} s'.format(timebla), 'latin-1'), fontsize = 14)
        ax.text(-70, +130, unicode('Event #: {0}'.format(evtnumber), 'latin-1'), fontsize = 14)
        ax.text(-70, +115, unicode('Cusp #: {0}'.format(cuspRunNumber), 'latin-1'), fontsize = 14)
        #ax.set_title('Event:')
        
        plt.tick_params(axis='both', which='major', labelsize=15)
        
        ################################################################################################
        # plot fitted tracks 

        
        positionmap = getPositionMap()
        pos_outer = positionmap[1]
        add = np.arange(32)
        pos_outer = np.column_stack((add, pos_outer))
        
        if trackcollection != None:
            for pars,tracks in zip(line_params, trackcollection):
                    x1 = 0.0
                    x2 = tracks[1][1][0] 
                    
                    if x2 is list:
                        x2 = x2[0]
                       
                    x2 = pos_outer[pos_outer[:,0] == x2][0][1]    
                                      
                               
                    y1 = line_2d(x1,pars)
                    y2 = line_2d(x2,pars)               
                    p0 = [x1,y1]
                    p1 = [x2,y2]
                    slope = slope_from_points(p0, p1)
                    intercept = p0[1] - slope*p0[0] 
                    
                                    
                    x = np.zeros(2)
                    x[0] = x1
                    x[1] = x2   
                    data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
                    line = mpl_lines.Line2D(x, data_y, color='red', lw = 1)
                    ax.add_line(line) 
                
                
                
                
                
        """for pars,points in zip(line_params, fitlines_points):
                x1 = vertex[0]
                #x2 = 
        
                x1 = 10
                x2 = -10
                
                y1 = line_2d(x1,pars)
                y2 = line_2d(x2,pars)               
                p0 = [x1,y1]
                p1 = [x2,y2]
                slope = slope_from_points(p0, p1)
                intercept = p0[1] - slope*p0[0] 
                   
                #print "points ", points 
                 
                x = np.zeros(2)
                x[0] = points[0][0]
                x[1] = points[2][0]
                if x[1] < 0:
                    x[1] = x[1] - 50
                else:
                    x[1] = x[1] + 50
                        

  
                data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
                line = mpl_lines.Line2D(x, data_y, color='red', lw = 1)
                ax[0].add_line(line)"""                     
        ################################################################################################     """  
        # plot intersection of lines         
        #for vert in vertex_points:
        #    plt.plot(vert[0], vert[1], 'o', color = 'red') 
        """if vertex != None:   
            vertex = np.array(vertex)                   
            ax[0].plot(vertex[0], vertex[1], 'o', color = 'green')         
          
        if vertex2 != None:
            vertex2 = np.array(vertex2)                               
            ax[0].plot(vertex2[0], vertex2[1], 'o', color = 'blue')       
        ################################################################################################     """
        ################################################################################################
        
        

        """theta = (90 - 70.0)*math.pi/180.0

        for i in range(-1000,1000):

            x = vertex2[0] + i*math.cos(theta)
            y = vertex2[1] + i*math.sin(theta)
            
            
            ax[0].scatter(x,y, s = 0.5)
        
        theta = (90 - 110.0)*math.pi/180.0

        for i in range(-1000,1000):

            x = vertex2[0] + i*math.cos(theta)
            y = vertex2[1] + i*math.sin(theta)
            
            
            ax[0].scatter(x,y, s = 0.5)"""
        
        
        
        
        
        
        
        #plt.savefig("./event_onlySMI_{0}_{1}.png".format(cuspRunNumber, evtnumber), facecolor=fig.get_facecolor()) 
        #{0}'.format(cuspRunNumber)  
        #plt.savefig("bgo_hit_map.pdf", facecolor=fig2.get_facecolor())
        #plt.show()
        #pdf.savefig( fig )
        #pdf.close()
        #cand_list.append([cuspRunNumber,evtnumber,BGOcharge, nI, nO, (startmixtime-ts),mtd])

        #print ts
        #count = count+1
            
        #print "total events: ", tot_events

        plt.savefig("./event_plots_2D/event2D_{1}_{0}.pdf".format(evtnumber,cuspRunNumber), facecolor=fig.get_facecolor())

        #return fig
        plt.show()


##############################################################################################################
def plot_event_hodoscope(ax, h):

    inner,outer,BG = h.drawHodoscope()

    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)
   
    ax.add_collection(BG   )
    ax.add_collection(inner)
    ax.add_collection(outer)
    
    """tpx_length = 28.0 # mm      
    tpx_hit_bounds = [[-0.5*tpx_length,0.5*tpx_length],[-0.5*tpx_length,0.5*tpx_length]]
    ext = [(0.5*tpx_length, 0.5*tpx_length),  (-0.5*tpx_length, 0.5*tpx_length), (-0.5*tpx_length, -0.5*tpx_length), (0.5*tpx_length, -0.5*tpx_length)]            
    tpx_polygon = Polygon(ext)     
    mpl_poly = PolygonPatch(tpx_polygon, alpha = 0.5, color = 'black')
    ax.add_patch(mpl_poly)"""
    
    # plot black circle to indicate BGO 
    r = 50.0
    circ = Circle([0,0], r)
    p = PatchCollection([circ], alpha=1.0)
    p.set_edgecolor("black")
    p.set_facecolor("None")
    ax.add_collection(p   )

    ax.set_aspect(1)    
    ax.set_ylabel("vertical position (mm)")
    ax.set_xlabel("horizontal position (mm)")
    ax.patch.set_facecolor='w'
    
##############################################################################################################
def plot_event_2D_tpx(evtnumber, cuspRunNumber, h, linepoints, trackcollection, tpx_polygon, 
line_params, fitlines_points, timebla, plot_trackfinding, vertex_points, vertex, nr_of_tracks):

        inner,outer,BG = h.drawHodoscope()
                   
        fig, ax = plt.subplots( frameon=False, figsize = (7,7))
        
        fig.patch.set_alpha(0.0)
        ax.set_xlim(-200,200)
        ax.set_ylim(-200,200)
       
        ax.add_collection(BG   )
        ax.add_collection(inner)
        ax.add_collection(outer)
        
        mpl_poly = PolygonPatch(tpx_polygon, alpha = 0.5, color = 'black')
        ax.add_patch(mpl_poly)

        ax.set_aspect(1)    
        ax.set_ylabel("vertical position (mm)", fontsize = 15)
        ax.set_xlabel("horizontal position (mm)", fontsize = 15)
        ax.patch.set_facecolor='w'
  
        
        ################################################################################################
        # plot lines from the rotation
        if plot_trackfinding[1] == True:
            plot_tracklines2D(linepoints,ax)
        ################################################################################################        
       
        plt.tight_layout()
   
        plt.tick_params(axis='both', which='major', labelsize=15)
        
        ################################################################################################
        # plot fitted tracks         
        positionmap = getPositionMap()
        pos_outer = positionmap[1]
        add = np.arange(32)
        pos_outer = np.column_stack((add, pos_outer))
        
        if trackcollection != None:
            for pars,tracks in zip(line_params, trackcollection):
                    x1 = 0
                    x2 = tracks[1][1][0] 
                    
                    if x2 is list:
                        x2 = x2[0]
                       
                    x2 = pos_outer[pos_outer[:,0] == x2][0][1]    
          
                    y1 = line_2d(x1,pars)
                    y2 = line_2d(x2,pars)               
                    p0 = [x1,y1]
                    p1 = [x2,y2]
                    slope = slope_from_points(p0, p1)
                    intercept = p0[1] - slope*p0[0] 
               
                    x = np.zeros(2)
                    x[0] = x1
                    x[1] = x2   
                    data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
                    line = mpl_lines.Line2D(x, data_y, color='red', lw = 1)
                    ax.add_line(line) 
                                  
        ################################################################################################     """  
        # plot intersection of lines         
        #for vert in vertex_points:
        #    plt.plot(vert[0], vert[1], 'o', color = 'red') 
        #print vertex
        if list(vertex) != [999,999]:   
            vertex = np.array(vertex)                   
            ax.plot(vertex[0], vertex[1], 'o', color = 'green')         
             
        ################################################################################################     """
        ################################################################################################


        plt.savefig("./event_plots/event_{1}_{0}.pdf".format(evtnumber,cuspRunNumber), facecolor=fig.get_facecolor())

        #return fig
        #plt.show()

##############################################################################################################
######################################################################################################################
# old 3 plotting, the new version (3 projections) is in trackfinidng_fitting_3D_final
# this is also cool though, it does a real 3d plot
def plot_event3D(cuspnr, evtnr, h, zposII, zposOO, zerrI, zerrO, bgox, bgoy, bgoz): 


    fig3D = fig = plt.figure(figsize=plt.figaspect(1)) 
    nI, nO, I, O = h.getActiveBar()
 
    #fig1 = plt.figure(figsize=(11,11))
    #ax3D = fig1.add_subplot(111, projection='3d')
    
    #print "boooooooooooooooooooo ", I
    #print "boooooooooooooooooooo ", O          
    # add a white circle for bgo in PMT map:
    r = 50.0
    c = Circle([0,0], r, alpha = 0.7)
    c.set_edgecolor('black')
    c.set_facecolor('black')  
    
    pmapcI, pmapcO = getPositionMap()
    pmapcI, pmapcO = pmapcI[I], pmapcO[O]     
   
    fxI = pmapcI[:,0]
    fxO = pmapcO[:,0]
    fyI = pmapcI[:,1]
    fyO = pmapcO[:,1]
   
    cornermapI, cornermapO = getCornerMap()    
    cmapI = np.array(cornermapI)
    
    cmapO = np.array(cornermapO)           
    cmap3d =  np.vstack((cmapI,cmapO))
    
    cmapI = cmapI[I]
    cmapO = cmapO[O]
    
    list_i = [i for i,x in enumerate(I) if x == True]
    
    list_o = [o for o,x in enumerate(O) if x == True]
    
    #print list_i
    #print list_o
    
    ax3D = fig3D.add_subplot(111, projection='3d')
            

    ax3D.add_patch(c)
    M3.art3d.pathpatch_2d_to_3d(c, z=bgoz, zdir="y")
    
    ax3D.set_xlabel("x (mm)", fontsize=10)
    ax3D.set_ylabel("z (mm)", fontsize=10) 
    ax3D.set_zlabel("y (mm)", fontsize=10)  
        
    ax3Dlim = 300
   
    #plt.yticks(fontsize = 13)
    
    #ax3D.add_patch(circle)
    #M3.art3d.pathpatch_2d_to_3d(circle, z=-10, zdir="y")
    bluelw = '0.4'
    bluelw2 = '0.8'
    bluems = '2'
    bluems2 = '0.9'
    blacklw = '0.3'
     
    ax3D.plot([bgox],[bgoz],[bgoy], linestyle="None", marker=".", c = "green",linewidth=bluelw, markersize = bluems)      
           
    for n, (rec,point) in enumerate(zip(cmapI,pmapcI)):
        zposI = zposII[n]
        
        zerr = zerrI
        ax3D.plot([point[0]], [zposI], [point[1]], linestyle="None", marker=".", c = "blue",linewidth=bluelw, markersize =bluems)
        #xerr, yerr = get_hit_error([rec])
        xerr, yerr = Calc_xy_error(list_i[n],0)
        ax3D.plot([point[0]+xerr, point[0]-xerr], [zposI, zposI], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'blue', markersize =bluems2)
        ax3D.plot([point[0], point[0]], [zposI, zposI], [point[1]+yerr, point[1]-yerr],  marker="_",linewidth=bluelw2, c = 'blue', markersize =bluems2)
        ax3D.plot([point[0], point[0]], [zposI+zerr, zposI-zerr], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'blue',  markersize =bluems2)
        
    for n, (rec,point) in enumerate(zip(cmapO,pmapcO)):
        zposO = zposOO[n]
        #print zposO
        zerr = zerrO
        ax3D.plot([point[0]], [zposO], [point[1]], linestyle="None", marker=".", c = "red",linewidth=bluelw, markersize =bluems)        
        #xerr, yerr = get_hit_error([rec])
        xerr, yerr = Calc_xy_error(list_o[n],1)
        ax3D.plot([point[0]+xerr, point[0]-xerr], [zposO, zposO], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'red', markersize =bluems2)
        ax3D.plot([point[0], point[0]], [zposO, zposO], [point[1]+yerr, point[1]-yerr],  marker="_",linewidth=bluelw2, c = 'red', markersize =bluems2)
        ax3D.plot([point[0], point[0]], [zposO+zerr, zposO-zerr], [point[1], point[1]],  marker="_",linewidth=bluelw2, c = 'red', markersize =bluems2)
    
    #zI, zO = h.getZPosition()             
    #zposI = zI[1]
    #zposO = zO[1]
        
    
    ax3D.set_xlim3d(-ax3Dlim, ax3Dlim)
    ax3D.set_ylim3d(-ax3Dlim, ax3Dlim)
    ax3D.set_zlim3d(-ax3Dlim, ax3Dlim)
        
    #plot all bars in black lines:
    for count,rec in enumerate(cmap3d):# rec[0] = four corner points (x,y), rec[1] number of bar
        for pt in rec: # loop over pt
            for pt2 in rec:

                if count < 32:  
                    zs1 = [135,135]  
                    zs2 = [-135,-135] 
                    zs3 = [135,-135] 
                    if (pt[0] == pt2[0] and pt[1] == pt2[1]):
                        ax3D.plot([pt[0], pt2[0]],zs3,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)
                   
                ax3D.plot([pt[0], pt2[0]],zs1,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)
                ax3D.plot([pt[0], pt2[0]],zs2,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)
                                     # rec[0] = four corner points, rec[1] number of bar  
                if count >=32:
                    dummycount = count - 32
                    zs1 = [235,235]
                    zs2 = [-235,-235] 
                    zs3 = [235,-235] 
                    if (pt[0] == pt2[0] and pt[1] == pt2[1]):
                        ax3D.plot([pt[0], pt2[0]],zs3,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)
                              
                ax3D.plot([pt[0], pt2[0]],zs1,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)
                ax3D.plot([pt[0], pt2[0]],zs2,[pt[1], pt2[1]], alpha = 0.1, c='black', linewidth=blacklw)



    # save axis settings
    axes_ticks = []
    axes_ticks.append(ax3D.xaxis.get_ticklocs())
    axes_ticks.append(ax3D.yaxis.get_ticklocs())
    axes_ticks.append(ax3D.zaxis.get_ticklocs())    


    ax3D.view_init(0, 90)
    ax3D.w_yaxis.line.set_lw(0.)
    ax3D.set_yticks([])
    # this turns off the axes:
    #ax3D._axis3don = False
   
    # this sets the background colour to white
    ax3D.grid(False)
    ax3D.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    fig3D.tight_layout()
    fig3D.savefig("front.pdf", bbox_inches='tight')
    ###########################################################
    ax3D.view_init(0, 0)
    ax3D.w_yaxis.line.set_lw(1.)
    ax3D.get_yaxis().set_ticks(axes_ticks[1])
    ax3D.w_xaxis.line.set_lw(0.)
    ax3D.set_xticks([])
    
    fO = h.fo.active
    fI = h.fi.active
    fi_dia = h.fi.diameter
    fo_dia = h.fo.diameter
    
    z_I_Fibre_hits = h.fi.getPositionMap()[fI]
    z_O_Fibre_hits = h.fo.getPositionMap()[fO]
    
    ToTs_i = (1.0*h.fi.ToT_first/127.0)[fI]
    ToTs_o = (1.0*h.fo.ToT_first/127.0)[fO]
    
    print h.fi.ToT    
    print h.fo.ToT
    
    
    
    print ToTs_i 
    print ToTs_o
    #exit(0)
    
    
    #patches = []
    chw = h.fi.chWidth
    for n, (f, c) in enumerate(zip(z_I_Fibre_hits, ToTs_i)):
        rec = Rectangle(xy = (f - chw*0.5, -fi_dia*0.5-chw - chw*0.5), width = chw, height = fi_dia + 2.0*chw, angle = 0.0 , alpha = c*0.8) # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
        #patches.append(rec)
        rec.set_facecolor('blue')  
        ax3D.add_patch(rec)  
        M3.art3d.pathpatch_2d_to_3d(rec, z=bgoz, zdir="x")

    chw = h.fo.chWidth
    for f, c in zip(z_O_Fibre_hits, ToTs_o):
        rec = Rectangle(xy = (f - chw*0.5, -fo_dia*0.5-chw - chw*0.5), width = chw, height = fo_dia + 2.0*chw, angle = 0.0 , alpha = c*0.8) # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
        #patches.append(rec)
        rec.set_facecolor('red')  
        ax3D.add_patch(rec)  
        M3.art3d.pathpatch_2d_to_3d(rec, z=bgoz, zdir="x")       
    #p = PatchCollection(patches, alpha=0.4)
    #p.set_array(np.array(colors))
    #ax3D.add_collection(p)   
    #ax3D.add_patch()
    #M3.art3d.pathpatch_2d_to_3d(p, z=bgoz, zdir="y")
    
    
    # this turns off the axes:
    #ax3D._axis3don = False
   
    # this sets the background colour to white
    ax3D.grid(False)
    ax3D.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    fig3D.tight_layout()
    fig3D.savefig("side.pdf", bbox_inches='tight') 
    ############################################################
    
    #for p in patches:
    #    p.remove()
    
    chw = h.fi.chWidth
    for f, c in zip(z_I_Fibre_hits, ToTs_i):
        rec = Rectangle(xy =  (-fi_dia*0.5 - chw*0.5, f - chw*0.5), height = chw, width = fi_dia + 2.0*chw, angle = 0.0 , alpha = c*0.8) # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
        #patches.append(rec)
        rec.set_facecolor('blue')  
        ax3D.add_patch(rec)  
        M3.art3d.pathpatch_2d_to_3d(rec, z=bgoz, zdir="z")

    chw = h.fo.chWidth
    for f, c in zip(z_O_Fibre_hits, ToTs_o):
        rec = Rectangle(xy = (-fo_dia*0.5 - chw*0.5, f - chw*0.5), height = chw, width = fo_dia + 2.0*chw, angle = 0.0 , alpha = c*0.8) # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
        #patches.append(rec)
        rec.set_facecolor('red')  
        ax3D.add_patch(rec)  
        M3.art3d.pathpatch_2d_to_3d(rec, z=bgoz, zdir="z")   
    
    ax3D.view_init(90, 90)
    ax3D.w_zaxis.line.set_lw(0.)
    ax3D.w_xaxis.line.set_lw(1.)
    ax3D.get_xaxis().set_ticks(axes_ticks[0])
  
    ax3D.set_zticks([])
    # this turns off the axes:
    #ax3D._axis3don = False
   
    # this sets the background colour to white
    ax3D.grid(False)
    ax3D.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax3D.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    fig3D.tight_layout()
    fig3D.savefig("top.pdf", bbox_inches='tight')   
    
    
    output = PdfFileWriter()

    # Appending two pdf-pages from two different files
    append_pdf(PdfFileReader(open("front.pdf","rb")),output)
    append_pdf(PdfFileReader(open("side.pdf","rb")),output)
    append_pdf(PdfFileReader(open("top.pdf","rb")),output)

    # Writing all the collected pages to a file
    output.write(open("./event_plots/front_side_top_{0}_{1}_n.pdf".format(cuspnr,evtnr),"wb"))
     

  
############################################################################################################## 
##############################################################################################################
def line_2d(x,p):
    y = -p[0]/p[1]*x - p[2]/p[1]    
    return y  
  
##############################################################################################################
def slope_from_points(point1, point2):
    return (point2[1] - point1[1])/(point2[0] - point1[0])
    
##############################################################################################################
##############################################################################################################
def plot_all_bars_2D(ax_pol, circle, union_cl_I, union_cl_O, vertex2) :
    cornermapI, cornermapO = getCornerMap()
    cmapI = np.array(cornermapI)
    cmapO = np.array(cornermapO)  
    #fig_pol, ax_pol = plt.subplots(figsize = (12,10))
    ax_pol.set_xlim(-200,200)
    ax_pol.set_ylim(-200,200)
    
    #ax_pol.set_xlim(-20,50)
    #ax_pol.set_ylim(-30,25)    
    
    for cornerI,cornerO in zip(cmapI,cmapO): 
        x0 = cornerI[2][0]
        y0 = cornerI[2][1]
        x1 = cornerI[0][0]
        y1 = cornerI[0][1]
        x2 = cornerI[1][0]
        y2 = cornerI[1][1]
        x3 = cornerI[3][0]
        y3 = cornerI[3][1]
        ext = [(x0, y0), (x1, y1), (x2, y2), (x3,y3)]
        polygon = Polygon(ext)
        mpl_poly = PolygonPatch(polygon, alpha = 0.1,color = 'black')
        ax_pol.add_patch(mpl_poly)

               
        x0 = cornerO[2][0]
        y0 = cornerO[2][1]
        x1 = cornerO[0][0]
        y1 = cornerO[0][1]
        x2 = cornerO[1][0]
        y2 = cornerO[1][1]
        x3 = cornerO[3][0]
        y3 = cornerO[3][1]
        ext = [(x0, y0), (x1, y1), (x2, y2), (x3,y3)]
        polygon = Polygon(ext)
        mpl_poly = PolygonPatch(polygon, alpha = 0.1, color = 'black')
        ax_pol.add_patch(mpl_poly)
  
    BGO_circle =  Circle(xy=(0.0 ,0.0), radius=50.0,  color = 'black', alpha = 0.5)
    ax_pol.add_patch(BGO_circle)
    ax_pol.set_xlabel("x (mm)")#, fontsize=30)
    ax_pol.set_ylabel("y (mm)")#, fontsize=30)  
    
    ax_pol.set_aspect(1)
    
    inner_fibres = Wedge((0,0), 167.0*0.5, 0, 360, width=4, color = 'gainsboro', alpha = 0.5, lw=0.1, ec = 'grey')
    ax_pol.add_patch(inner_fibres)
    outer_fibres = Wedge((0,0), 293.0*0.5, 0, 360, width=4, color = 'gainsboro', alpha = 0.5, lw=0.1, ec = 'grey')
    ax_pol.add_patch(outer_fibres)
    
    #plt.yticks(fontsize = 20)
    #plt.xticks(fontsize = 20)
    #ax_pol.add_patch(circle)    
    #plot_tracklines2D(linepoints,ax_pol)
    
    """positionmap = getPositionMap()
    pos_outer = positionmap[1]
    add = np.arange(32)
    pos_outer = np.column_stack((add, pos_outer))
    
    if trackcollection != None:
        for pars,tracks in zip(line_params, trackcollection):
                x1 = vertex2[0]
                x2 = tracks[1][1][0] 
                
                if x2 is list:
                    x2 = x2[0]
                   
                x2 = pos_outer[pos_outer[:,0] == x2][0][1]    
                

                
                           
                y1 = line_2d(x1,pars)
                y2 = line_2d(x2,pars)               
                p0 = [x1,y1]
                p1 = [x2,y2]
                slope = slope_from_points(p0, p1)
                intercept = p0[1] - slope*p0[0] 
                
                                
                x = np.zeros(2)
                x[0] = x1
                x[1] = x2   
                data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
                line = mpl_lines.Line2D(x, data_y, color='red', lw = 1)
                ax_pol.add_line(line)""" 
    if union_cl_I[0] !=  999:
        for up in union_cl_I:       
            if up.geom_type == 'MultiPolygon':
                for polygon in up:
                    mpl_poly = PolygonPatch(polygon, alpha = 0.6)
                    ax_pol.add_patch(mpl_poly)

            elif up.geom_type == 'Polygon':
                mpl_polys = PolygonPatch(up, alpha = 0.6)
                ax_pol.add_patch(mpl_polys)
    if union_cl_O[0] !=  999:        
        for up in union_cl_O:       
            if up.geom_type == 'MultiPolygon':
                for polygon in up:
                    mpl_poly = PolygonPatch(polygon, alpha = 0.6)
                    ax_pol.add_patch(mpl_poly)

            elif up.geom_type == 'Polygon':
                mpl_polys = PolygonPatch(up, alpha = 0.6)
                ax_pol.add_patch(mpl_polys)
    

    #plt.savefig("tescht.pdf")
##############################################################################################################
def plot_tracklines2D(linepoints, ax):      
    linepoints =  np.array(linepoints)
    p0 = [0,0]
    p1 = [0,0]
    for points in linepoints:        
        p0[0] = points[0]
        p0[1] = points[1]
        p1[0] = points[2]
        p1[1] = points[3]
        slope = slope_from_points(p0, p1)
        intercept = p0[1] - slope*p0[0] 
        x = ax.get_xlim()
        #x[0] = 
        #y = ax.get_ylim()
        data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
        #line = mpl_lines.Line2D(x, data_y, color='red', lw = 0.1)
        ax.plot(x,data_y, lw = 0.5, alpha = 0.05, color = 'black')
        #ax.add_line(line) 



##############################################################################################################
def plot_event2D(ax,h,B):
    inner,outer,BG = h.drawHodoscope()
    orig = np.copy(B.dataMap)
    inner.set_edgecolor("w")
    inner.set_cmap(gray_r)
    outer.set_edgecolor("w")    
    outer.set_cmap(gray_r)

    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)

    bgoPainter = BGOPainter(gray_r)
    a,b,c,d = bgoPainter.plot(B)
    a.set_edgecolor("w")
    b.set_edgecolor("w")
    c.set_edgecolor("w")
    d.set_edgecolor("w")
    BG.set_facecolor("w")
    BG.set_alpha(0.001)
    ax.add_collection(a    )
    ax.add_collection(b    )
    ax.add_collection(c    )
    ax.add_collection(d    )
    ax.add_collection(BG)
    
    ax.add_collection(inner)
    ax.add_collection(outer)

    axlim = 200
    ax.set_xlim(-axlim,axlim)
    ax.set_ylim(-axlim,axlim)

    r = 50.0
    c = Circle([0,0], r)
    p = PatchCollection([c], alpha=0.5)
    p.set_edgecolor("g")
    ax.add_collection(p   )
  
