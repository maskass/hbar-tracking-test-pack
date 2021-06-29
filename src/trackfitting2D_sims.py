
import ROOT
from ROOT import TF3, Fit #, Minuit2

import matplotlib.pyplot as plt
from src.hodoscope_trackhelper import *
import math

##############################################################################################################
def det(a, b):
    return a[0] * b[1] - a[1] * b[0]
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
        print "x", x
        print "y", y             
        #x = sum(x)/len(x)
        #y = sum(y)/len(y)
        #z = sum(z)/len(z)
    
    else:
        #print "else!"
        x = posmap[barnr][0][0]
        y = posmap[barnr][0][1]

    return x,y
##############################################################################################################
def get_hit_position2D_sims(posmap, barnr, ch): # barnr is a list with the hitbars
                                                
    if len (barnr) > 1 :
        new_ch = []
        ch = np.array(ch)
        ch = ch.transpose()
        
        new_ch = ch[barnr]

        x = [posmap[i][0] for i in barnr]
        y = [posmap[i][1] for i in barnr]  

        x = np.dot(new_ch,x) / np.sum(new_ch)
        y = np.dot(new_ch,y) / np.sum(new_ch)   
           
        #x = sum(x)/len(x)
        #y = sum(y)/len(y)
        #z = sum(z)/len(z)
    
    else:
        #print "else!"
        x = posmap[barnr][0][0]
        y = posmap[barnr][0][1]

    return x,y
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
class chi2fct2d( ROOT.TPyMultiGenFunction ):
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
def dist_p_line2d(x,y,p):    
    dist = abs(p[0]*x + p[1]*y + p[2])/(math.sqrt(p[0]*p[0] + p[1]*p[1]))
    return dist*dist

##############################################################################################################
def line_2d(x,p):
    y = -p[0]/p[1]*x - p[2]/p[1]    
    return y  
##############################################################################################################    
def angle_btw_lines(line_paras, trackscoll):   # only for two lines atm!
    angle = 0
    k1 = -line_paras[0][0]/line_paras[0][1]
    k2 = -line_paras[1][0]/line_paras[1][1]
    
    theta1 = math.degrees(math.atan(k1))
    theta2 = math.degrees(math.atan(k2))

    angle = theta1 - theta2

                #if angle < 160:
                #    continue
    if (theta1 < 0 and theta2 < 0) or (theta1 > 0 and theta2 > 0):
        angle = 180 - abs(angle)
        
    print "TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEST"
    print "ANGLE: ",theta1, theta2, theta1 - theta2, "degree" 
    print "---------------------------------------------------------------"
    #angle = theta1 - theta2
    return [abs(angle),999,999]
     
##############################################################################################################    
def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]
############################################################################################################## 
def angle_2_lines(hitcoll, vertex):

        angle = [999, 999, 999]

        tr = hitcoll[0]
        tr2 = hitcoll[1]

        positionmap = getPositionMap()
    
        pos_inner = positionmap[0]
        pos_outer = positionmap[1]
        add = np.arange(32)

        pos_inner = np.column_stack((add, pos_inner))
        pos_outer = np.column_stack((add, pos_outer))


        merge_o = tr[1][1]
        merge_o2 = tr2[1][1]
        
        ro_merge = pos_outer[merge_o]
        ro_merge2 = pos_outer[merge_o2]
        ro_new = np.mean(ro_merge[:,0:], axis = 0)
        ro_new2 = np.mean(ro_merge2[:,0:], axis = 0)
                     

        vA = [(vertex[0]-ro_new[1]), (vertex[1]-ro_new[2])]
        vB = [(vertex[0]-ro_new2[1]), (vertex[1]-ro_new2[2])]

        dot_prod = dot(vA, vB)

        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
 
        cos_ = dot_prod/magA/magB
                    
        if abs(cos_) > 1:
            return [999, 999, 999]
        else:                    
            angle = math.acos(cos_)
            angle = math.degrees(angle)%360
            return [angle, 999, 999]

############################################################################################################## 
############################################################################################################## 
def angle_3_lines(line_paras, trackscoll):

    angles = []
    track1 = trackscoll[0]
    track2 = trackscoll[1]
    track3 = trackscoll[2]
    
    line_p1 = line_paras[0]
    line_p2 = line_paras[1]
    line_p3 = line_paras[2]
    
    combs = []
    combs.append([track1, track2])
    combs.append([track2, track3])
    combs.append([track1, track3])
    
    line_ps = []
    line_ps.append([line_p1, line_p2])
    line_ps.append([line_p2, line_p3])
    line_ps.append([line_p1, line_p3])
    
    
    for i,j in zip(combs,line_ps):
        track1 = i[0]
        track2 = i[1]
        #print "track1", track1
        #print "track2,", track2

        x1 = track1[1][0]
        x2 = track1[2][0]
        
        y1 = line_2d(x1, j[0])
        y2 = line_2d(x2, j[0])
        
        x3 = track2[1][0]
        x4 = track2[2][0]
        
        y3 = line_2d(x3, j[1])
        y4 = line_2d(x4, j[1])
        # Get nicer vector form
        vA = [(x1-x2), (y1-y2)]
        vB = [(x3-x4), (y3-y4)]
        # Get dot prod
        dot_prod = dot(vA, vB)
        # Get magnitudes
        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
        # Get cosine value
        cos_ = dot_prod/magA/magB
        
        if abs(cos_) > 1:
            angles.append(999)
            continue
            
        
        # Get angle in radians and then convert to degrees
        angle = math.acos(cos_)
        # Basically doing angle <- angle mod 360
        ang_deg = math.degrees(angle)%360

        if ang_deg-180>=0:
        #    # As in if statement
            angles.append(360 - ang_deg)
        else: 
            angles.append(ang_deg)







        """k1 = -j[0][0]/j[0][1]
        k2 = -j[1][0]/j[1][1]
    
        theta1 = math.degrees(math.atan(k1))
        theta2 = math.degrees(math.atan(k2))

        angle = theta1 - theta2

                #if angle < 160:
                #    continue
        if (theta1 < 0 and theta2 < 0) or (theta1 > 0 and theta2 > 0):
            angle = 180 - abs(angle)

        angles.append(angle)"""
    
    return angles
            
############################################################################################################## 
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
        ##print "--- AFTER CLUSTER CHECK: ", ctr
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
                    xerr, yerr = Calc_xy_error(barnr,hit[0])
                                   
                if int(hit[0]) == 1:
                    #get x,y,z pos of cluster-> mean val ->add to trackpos
                    x,y = get_hit_position2D(posmapO, barnr, outerA)
                    xerr, yerr = Calc_xy_error(barnr, hit[0])
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

                errorfct = chi2fct2d()

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
                #print "ooooooooooooooooooooooooooooo ", minval, minimizer.Status()"""

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
           

##################################################################################################################################
def do_le_2Dfitting_sims(hitcoll, I, O, vertex, vertex_err, innerA, outerA):
    ##print " "
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
        ##print "--- AFTER CLUSTER CHECK: ", ctr
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
                    x,y = get_hit_position2D_sims(posmapI, barnr, innerA)
                    xerr, yerr = Calc_xy_error(barnr,hit[0])
                                   
                if int(hit[0]) == 1:
                    #get x,y,z pos of cluster-> mean val ->add to trackpos
                    x,y = get_hit_position2D_sims(posmapO, barnr, outerA)
                    xerr, yerr = Calc_xy_error(barnr, hit[0])
                    #print "errrrrrr", xerr, yerr
                    
                print x, y, xerr, yerr    
                trackpoints.append([x,y,xerr,yerr])
            
        trackscoll.append(trackpoints) 
        
    number_of_tracks = len(trackscoll)   
    
    #print "trackscoll: "
    #for i in trackscoll:
    #    print i
    
    #save_results = []
             
    for n, currtrack in enumerate(trackscoll):
        print currtrack   
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

                errorfct = chi2fct2d()

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
                #print "ooooooooooooooooooooooooooooo ", minval, minimizer.Status()"""

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
           

def Calc_xy_error(barnr, layer):
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
