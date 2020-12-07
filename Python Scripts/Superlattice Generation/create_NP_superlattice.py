#!/usr/bin/env python
# coding: utf-8
import numpy as np
import math
import os
from IPython.display import clear_output
import timeit
import itertools as it
#np.set_printoptions(threshold=np.nan)

#the most general lattice, all other bravais lattices can be derived from this one
class triclinic_lattice:
    def __init__(self,latticeType,currDiameter,paramDict):
        self.npnumber = paramDict['npnumber']
        self.xNum = paramDict['xNum']
        self.yNum = paramDict['yNum']
        self.zNum = paramDict['zNum']
        self.dSTD = paramDict['dSTD']
        self.lSTD = paramDict['lSTD']
        self.elDensity = paramDict['electronDensity']
        self.vDensity = paramDict['vacancyDensity']
        self.gbNumber = paramDict['grainBoundaryNumber']
        self.tbNumber = paramDict['twinBoundaryNumber']
        self.crack = paramDict['crackBool']
        self.rotate = paramDict['rotateBool']
        self.grain_110 = paramDict['110Grain']
        self.gAngles = paramDict['grainAngles']
        self.voidCutoff = paramDict['voidCutoff']
        self.dDiameter = currDiameter
        self.lType = latticeType
        self.crack_length = paramDict['crackLength']
        self.trim = paramDict['trimmingMethod']
        self.boundaryInfo = (0,0)
        self.average_ll = 0.5 #expect ligand length total between np to be 1 nm, so per np = 0.5 nm
        #xiaolei data: diameter = 6 +- 0.4
        self.NP_grain_compression_overlap_thr = 0.05*self.dDiameter #was 0.5*self.dDiameter
        self.NP_twin_compression_overlap_thr = 0.25*self.dDiameter
        self.min_compression_thr_tol = 1e-2
        
        #standard deviation of diameter
        if paramDict['stdMethod'] == "percent":
            self.np_d_std = self.dDiameter*self.dSTD
        elif paramDict['stdMethod'] == "fixed":
            self.np_d_std = self.dSTD
            
        #define upper and lower bounds of size disorder
        self.np_ub = self.dDiameter + 2*self.average_ll
        self.np_lb = self.dDiameter - 2*self.average_ll
        
        #upper and lower bounds of position disorder
        self.npl_ub = self.average_ll/2
        self.npl_lb = -self.average_ll/2
        
        #define where the first nanoparticle will sit
        a = self.dDiameter/2.0 + self.average_ll
        self.baseArray = [a,a,a]
        
        self.np_array = []
        self.np_id_array = []
        self.zNumAdjust = 0
        
    #a method to define all the lattice angles
    def setLatticeAngles(self):
        self.alpha = 99*math.pi/180
        self.beta = 99*math.pi/180
        self.gamma = 99*math.pi/180
      
    #a method to create the lattice vectors from the lattice angles and given spacing
    def setLatticeVectors(self, latticeSpacing):
        (a,b,c) = latticeSpacing
        
        self.a_1 = [0,0,b]
        self.a_2 = [a*math.sin(self.gamma),0,a*math.cos(self.gamma)]
        
        x_3 = c*(math.cos(self.beta)-math.cos(self.alpha)*math.cos(self.gamma))/math.sin(self.gamma)
        z_3 = c*math.cos(self.alpha)
        y_3 = math.sqrt(c**2-x_3**2-z_3**2)
        self.a_3 = [x_3,y_3,z_3]
       
    #a method to set the basis vectors
    def setSpacingLattice(self):
        self.setLatticeAngles()
        
        #in the typical a,b,c representation, z will be in the b direction, x in a, and y in c
        self.sidelength = self.dDiameter + 2*self.average_ll
        latticeSpacing = (self.sidelength, self.sidelength, self.sidelength)
        
        self.setLatticeVectors(latticeSpacing)
        
        #adjust the zNum in order to cut the sample into a cubic shape
        if self.trim == "cubic" and self.lType == 'tric':
            deltaZ = self.a_2[2]*self.xNum
            self.zNumAdjust = int(deltaZ/self.a_1[2])
            self.zNum += abs(self.zNumAdjust)
        
    def setSpacing011Plane(self):
        zVector = np.asarray(self.a_1) + np.asarray(self.a_3)
        zLength = math.sqrt(np.dot(zVector, zVector))
        proj_xzLength = math.sqrt((self.a_3[2] + self.a_1[2])**2 + self.a_3[0]**2)
        
        gammaCorrection = math.asin(abs(self.a_3[0])/(self.a_1[2] + self.a_3[2]))
        alphaCorrection = math.acos(proj_xzLength/zLength)
        
        self.alpha = self.alpha - alphaCorrection
        self.gamma = self.gamma + gammaCorrection
        
        self.setLatticeVectors([self.sidelength, zLength,self.sidelength])
        
        
    #defining the position where the lattice will be generated from
    def initiatePositions(self):
        self.curr_x = self.baseArray[0]
        self.curr_y = self.baseArray[1]
        self.curr_z = self.baseArray[2]

        self.cellx = self.baseArray[0]
        self.celly = self.baseArray[1]
        self.cellz = self.baseArray[2]
        
    #a method to add a nanoparticle to the nanoparticle solid. Called by generateLattice and add_zboundary_randomness
    def appendNP(self, npArray):
        #random size
        self.np_diameter = np.random.normal(self.dDiameter,self.np_d_std)
        while self.np_diameter < self.np_lb or self.np_diameter > self.np_ub:
            self.np_diameter = np.random.normal(self.dDiameter,self.np_d_std)
        
        #to ensure no overlap a priori, radius + jitter cannot exceed desired radius + ligand length
        upperl_bound = (self.np_ub - self.np_diameter)/2.0
        iterations = 0
        
        #random location
        self.x_jit = np.random.normal(0,self.lSTD)
        self.y_jit = np.random.normal(0,self.lSTD)
        self.z_jit = np.random.normal(0,self.lSTD)
        self.total_jit = math.sqrt(self.x_jit**2 + self.y_jit**2 + self.z_jit**2)
        while self.total_jit > upperl_bound:
        #while self.total_jit < -upperl_bound or self.total_jit > upperl_bound:
            self.x_jit = np.random.normal(0,self.lSTD)
            self.y_jit = np.random.normal(0,self.lSTD)
            self.z_jit = np.random.normal(0,self.lSTD)
            self.total_jit = math.sqrt(self.x_jit**2 + self.y_jit**2 + self.z_jit**2)
            iterations += 1
            if iterations > 20:
                self.x_jit = 0.0
                self.y_jit = 0.0
                self.z_jit = 0.0
                break
        
        #append nanoparticle
        npArray.append([self.curr_x + self.x_jit,self.curr_y+self.y_jit,self.curr_z+self.z_jit,self.np_diameter, self.x_jit, self.y_jit, self.z_jit])
        return npArray
        
    #a method to create the nanoparticle solid
    #the id array is used for adding randomness to the left z-edge (needed for sorting the lattice)
    #cell_edge values are used for calculating sample volume
    def generateLattice(self):
        npArray = []
        idArray = []
        base_id = -1
        id = 0
        
        base_x = self.baseArray[0]
        base_z = self.baseArray[2]
        base_z_y = self.baseArray[2]
        
        for y in range(1,self.yNum+1):
            for x in range(1,self.xNum+1):
                base_id += 1
                id = base_id
                for z in range(1,self.zNum + 1):
                    npArray = self.appendNP(npArray)
                    idArray.append(id)
                    id+=self.xNum
                    if y == 1 and x == 1 and z == self.zNum - abs(self.zNumAdjust):
                        #need the side length for Monte Carlo simulations
                        #self.cellz = self.curr_z - self.baseArray[2] + self.np_ub
                        self.cellz = z*self.sidelength
                        #due to triclinic nature, for volume of cell we want the full length of the sides
                        self.cellz_edge = (z - 1)*(self.sidelength) + self.np_ub
                        
                    if y ==1 and x==1 and z == self.zNum:
                        self.cellz_adjusted = self.curr_z - self.baseArray[2] + self.np_ub
                        self.cellz_edge_adj = (z - 1)*(self.sidelength) + self.np_ub
                        
                    
                    self.curr_x += self.a_1[0]
                    self.curr_y += self.a_1[1]
                    self.curr_z += self.a_1[2]
                
                if y == 1 and x == self.xNum:
                    self.cellx = self.curr_x - self.cellx + self.np_ub + self.x_jit
                    
                    #length of side in the x direction
                    self.cellx_edge = (self.xNum-1)*(self.sidelength) + self.np_ub + self.x_jit
                   
                base_z += self.a_2[2]
                
                self.curr_z = base_z
                self.curr_x += self.a_2[0]
                self.curr_y += self.a_2[1]
                
            base_x += self.a_3[0]
            base_z_y += self.a_3[2]
            base_z = base_z_y
            
            self.curr_x = base_x
            self.curr_z = base_z
            
            self.curr_y += self.a_3[1]
            
            base_id += (self.xNum*self.zNum-base_id)
            
        self.celly = self.curr_y - self.a_3[1] - self.celly + self.np_ub + self.y_jit
        
        #length of side in the y direction
        self.celly_edge = (self.yNum-1)*(self.sidelength) + self.np_ub + self.y_jit
        
        return npArray,idArray
        
    #a method to calculate the volume of the nanoparticle solid
    def calculateVolume(self):
        return self.cellx_edge*self.celly_edge*self.cellz_edge*math.sqrt(1-math.cos(self.alpha)**2 - math.cos(self.beta)**2-math.cos(self.gamma)**2 + 2*math.cos(self.alpha)*math.cos(self.beta)*math.cos(self.gamma))
    
    #return which y-layer a nanoparticle belongs to, critical if including location disorder
    def calculateYLayer(self, y):
        return int((y - self.baseArray[1])/self.a_3[1] + 0.5) #note: need to add 0.5 to round (0.5,1) up to 1
    
    #a method to remove individual nanoparticles at random
    def createVacancies(self,np_array,id_array):
        indexArray = [] #to pop nanoparticle ids out when relevant
        numberVacancies = int(self.vDensity*len(np_array))
        for x in range(0,numberVacancies):
            randIndex = np.random.randint(0,len(np_array))
            #make sure that the nanoparticle being removed is not at either end in the z direction (will affect transport code)
            zValue = np_array[randIndex][2]
            while zValue < self.dDiameter or zValue > (self.cellz - self.dDiameter):
                randIndex = np.random.randint(0,len(np_array))
                zValue = np_array[randIndex][2]
            np_array.pop(randIndex)
            id_array.pop(randIndex)
        return np_array, id_array
    
    #a method to define where the grains will seed from in the lattice
    def defineGrains(self):
        self.grainList = []
        
        if(self.gbNumber > 1):
            #as a crude way of spacing out grains, make sure they occur at regular intervals in z
            for x in range(0,self.gbNumber):
                randIndex = np.random.randint(0,len(self.np_array))
                yValue = self.np_array[randIndex][1]
                zValue = self.np_array[randIndex][2]
                grainSpace = self.cellz/(self.gbNumber+1)
                while yValue > self.dDiameter or zValue < (x+1)*grainSpace or zValue > (x+3/2)*(grainSpace):
                    randIndex = np.random.randint(0,len(self.np_array))
                    zValue = self.np_array[randIndex][2]
                    yValue = self.np_array[randIndex][1]

                self.grainList.append(randIndex)
        elif(self.gbNumber == 1):
            randIndex = int(len(self.np_array)/4 - self.zNum/2) #want our nanoparticle to be in the middle of the bottom layer
            self.grainList.append(randIndex)
            
    #a method to define where the twin will seed from in the lattice
    def defineTwinAngle(self,array,currIndex):
        completelyRandom = False
        currNP = array[currIndex]
        x0 = currNP[0]
        z0 = currNP[2]
        
        if completelyRandom:
            #define perfect symmetry angle lines going from one corner to the other. These will not give twins, but
            #regular lattices
            offset = abs(self.gamma*180/math.pi - 90)/2
            sAngle1 = -(45 - offset)*math.pi/180
            sAngle2 = (45 + offset)*math.pi/180
            isSAngle = False

            #want to choose a nanoparticle at random in the top layer
            spread = 11
            indices = []
            for i in range(spread):
                rangeMiddle = (i - spread//2)*self.zNum + currIndex
                if i!= spread//2:
                    indices.extend(range(rangeMiddle - spread//2, rangeMiddle + spread//2 + 1))
                else:
                    indices.extend([currIndex - 1, currIndex + 1])

            randIndex = np.random.choice(indices)

            twinNP = array[randIndex]
            xTwin = twinNP[0]
            zTwin = twinNP[2]

            deltaX = xTwin - x0
            deltaZ = zTwin - z0


            vertical = np.isclose(deltaZ, 0.0, atol = 2e-1) and np.isclose(offset, 0.0)
            horizontal = np.isclose(deltaX, 0.0, atol = 2e-1)
            if not (horizontal):
                #calculate theta and the twin line's slope
                theta = -math.atan(deltaZ/deltaX)
                isSAngle = (np.isclose(theta, sAngle1) or np.isclose(theta, sAngle2))
            else:
                theta = math.pi/2
                isSAngle = True

            if abs(theta) > 80.0*math.pi/180.0:
                    isSAngle = True


            while(randIndex == currIndex or isSAngle):
                randIndex = np.random.choice(indices)
                twinNP = array[randIndex]
                xTwin = twinNP[0]
                zTwin = twinNP[2]

                deltaX = xTwin - x0
                deltaZ = zTwin - z0

                vertical = np.isclose(deltaZ, 0.0, atol = 2e-1) and np.isclose(offset, 0.0)
                horizontal = np.isclose(deltaX, 0.0, atol = 2e-1)
                if not (horizontal):
                    #calculate theta and the twin line's slope
                    theta = -math.atan(deltaZ/deltaX)
                    isSAngle = (np.isclose(theta, sAngle1) or np.isclose(theta, sAngle2))
                else:
                    theta = math.pi/2
                    isSAngle = True

                if abs(theta) > 80.0*math.pi/180.0:
                    isSAngle = True
        else:
            #indices which will return good lattice angles
            #these modifiers don't all bisect
            #modifiers = [77, 84, 81, -198, -200, -37, -121, -203, 116, 35, -75, -119, 118, -115, -204, 43, 155, 85]
            #these modifiers all bisect
            modifiers = [81, -198, -200, -121, -203, -119, 118, -204]
            indexModifier = np.random.choice(modifiers)
            randIndex = currIndex + indexModifier
            twinNP = array[randIndex]
            xTwin = twinNP[0]
            zTwin = twinNP[2]
            theta = -math.atan((zTwin - z0)/(xTwin - x0))
            
            #print("Theta is: ",theta*180/math.pi, " and the modifier was: ",indexModifier)
        
        return theta, randIndex
        
    #a method to create a void between two grains
    def voidSeparation(self):
        voidBoolFloat = np.random.uniform()
        
        if voidBoolFloat > (1-self.voidCutoff):
            sepMean = 2.5*self.dDiameter
            sepSTD = 0.5*self.dDiameter
            gSeparation = np.abs(np.random.normal(sepMean,sepSTD))
        else:
            gSeparation = 0
            
        return gSeparation
     
    #delete points along and over grain line
    #will remove points below or above line in the x-direction depending on rotation angle theta
    #will remove nanoparticles that intersect in the two arrays, from the grainArray (not rotated array)
    def removeNanoparticlesOnGrainLine(self, theta,grainArray,rotated_np_array,lineValues):
        newArray = []
        z1 = lineValues[0]
        x1 = lineValues[1]
        m1 = lineValues[2]
        line = np.asarray([x1,z1])
        lineVector = np.asarray([1, m1])
        indexArray = []
        removalArray = []
        
        #allow for void space between grains
        tSeparation = self.voidSeparation()
        
        compareArray = np.asarray(rotated_np_array)
        
        for index, i in enumerate(grainArray):
            #get all the coordinates as variables x,y,z
            x = i[0]
            y = i[1]
            z = i[2]
            yN = self.calculateYLayer(y) #y layer number, starting with 0
            remove = False
            
            x_boundary = (z-(z1+yN*self.a_3[2]))/m1 + (x1+yN*self.a_3[0])
            
            if theta > 0 and x >= x_boundary + 0.5*self.dDiameter:
                remove = True
            elif theta < 0 and x  <= x_boundary - 0.5*self.dDiameter:
                remove = True   
                
            #to create a void, calculate the distance from nanoparticle to grain boundary
            if tSeparation > 0:
                point = np.asarray([x,z])
                projPointOntoLine = line + lineVector*np.dot(point - line,lineVector)/np.dot(lineVector,lineVector)
                distV = projPointOntoLine - point
                dist = math.sqrt(np.dot(distV,distV))    
                if dist < tSeparation:
                    remove = True

            #remove on line in grainArray, by checking if the nanoparticles overlap directly between two arrays
            if not remove:
                newArray.append(i)
                indexArray.append(index) 
            else:
                removalArray.append(index)

        return newArray, removalArray, indexArray
    
    #delete points along and over twin line
    def removeNanoparticlesOnTwinLine(self, theta,twinArray,lineValues,direction):
        newArray = []
        removalArray = []
        edgeArray = []
        indexArray = []
        
        z1 = lineValues[0]
        x1 = lineValues[1]
        m1 = lineValues[2]
        line = np.asarray([x1,z1])
        lineVector = np.asarray([1, m1])
        
        for index,i in enumerate(twinArray):
            #get all the coordinates as variables x,y,z
            x = i[0]
            y = i[1]
            z = i[2]
            yN = self.calculateYLayer(y) - (self.yNum - 1) #y layer number, ending at 0 (so ..., -2, -1, 0)
            remove = False
            edge = False
            
            #calculate distance to the twin line. Only will remove nanoparticles if the distance 
            #of the nanoparticle to the twin line is greater than a half diameter
            point = np.asarray([x,z])
            projPointOntoLine = line + lineVector*np.dot(point - line,lineVector)/np.dot(lineVector,lineVector)
            distV = projPointOntoLine - point
            dist = np.linalg.norm(distV)
            
            x_boundary = (z-(z1+yN*self.a_3[2]))/m1 + (x1+yN*self.a_3[0])

            if direction == 'right':
                
                if theta > 0:
                    if x >= x_boundary + 0.5*self.dDiameter:
                        remove = True
                    elif dist < 2*self.dDiameter:
                        edge = True
                elif theta < 0:
                    if x <= x_boundary - 0.5*self.dDiameter:
                        remove = True
                    elif dist < 2*self.dDiameter:
                        edge = True
            elif direction == 'left':
                if theta > 0:
                    if x <= x_boundary - 0.5*self.dDiameter:
                        remove = True
                    elif dist < 2*self.dDiameter:
                        edge = True
                elif theta < 0:
                    if x >= x_boundary + 0.5*self.dDiameter:
                        remove = True  
                    elif dist < 2*self.dDiameter:
                        edge = True
            
            if not remove:
                newArray.append(i)
                indexArray.append(index)
                if edge:
                    edgeArray.append(index)
            else:
                removalArray.append(index)

        return newArray, removalArray, indexArray
    
    #remove points outside crystal bounds, for a specified array
    def removeNanoparticlesOutside(self,grainArray):
        newArray = []
        removalArray = []
        
        #right z-edge line
        x2 = self.baseArray[0]
        z2 = self.cellz - self.baseArray[2]
        xEdgeAngle = self.gamma

        if xEdgeAngle != math.pi/2:
            m2 = 1/math.tan(xEdgeAngle)
        else:
            m2 = 0

        #left z-edge line
        x3 = self.baseArray[0]
        z3 = self.baseArray[2]
        m3 = m2
        
        if self.trim == 'tric':
            for i in grainArray:
                x = i[0]
                y = i[1]
                z = i[2]
                radius = i[3]/2
                #print(y)
                yN = self.calculateYLayer(y)
                
                position_disorder = (self.np_ub - i[3])/2.0

                #Standard definition of crystal bounds  as linear lines. yN gives the y layer number
                #so that the layers can be shifted. Using the location disorder upper and lower bounds to allow
                #for jitter if applicable
                leftz = m3*(x-(x3+yN*self.a_3[0])) + yN*self.a_3[2] - .001
                rightz = m2*(x-(x2+yN*self.a_3[0])) + (z2 + z3 +yN*self.a_3[2]) + .001
                topx = self.cellx + yN*self.a_3[0] + .001
                bottomx = yN*self.a_3[0] - .001

                if not((z - radius) < leftz or (z + radius) > rightz or (x + radius) > topx or (x - radius) < bottomx) :
                    newArray.append(i)
                    
        elif self.trim == 'cubic':
            
            for i in grainArray:
                x = i[0]
                y = i[1]
                z = i[2]
                yN = self.calculateYLayer(y)
                
                if self.zNumAdjust <= 0:
                    leftz = self.baseArray[2] + yN*self.a_3[2] + self.npl_lb - .001
                    rightz = self.cellz - self.baseArray[2] +yN*self.a_3[2] + self.npl_ub + .001
                    topx = self.cellx - self.baseArray[0] + yN*self.a_3[0] + self.npl_ub + .001
                    bottomx = self.baseArray[0] + yN*self.a_3[0] +self.npl_lb - .001
                else:
                    leftz = self.baseArray[2] + yN*self.a_3[2] + self.xNum*self.a_2[2] + self.npl_lb - .001
                    rightz = self.cellz - self.baseArray[2] + self.xNum*self.a_2[2] + yN*self.a_3[2] + self.npl_ub + .001
                    topx = self.cellx - self.baseArray[0] + yN*self.a_3[0] + self.npl_ub + .001
                    bottomx = self.baseArray[0] + yN*self.a_3[0] +self.npl_lb - .001
                
                if not(z < leftz or z > rightz or x > topx or x < bottomx) :
                    newArray.append(i)
                    
        return newArray
    
            
    #a method to add unevenness to the left z-boundary of an array
    #used in the formation of grain boundaries
    def add_zBoundary_randomness(self, npList, idList,size, phi):
        
        #add randomness without phi included, so that magnitude of rotation is correct (seems to be off)
        sorted_npArray = [x for _,x in sorted(zip(idList,npList))]
        sorted_idList = sorted(idList)
        npList[:] = sorted_npArray
        
        baseSize = int(len(npList)/self.yNum)
        workingList = npList[0:size] + npList[baseSize:(baseSize+size)]
        
        for index, j in enumerate(workingList):
            
            r = np.random.random_sample()
            if r <= 0.5:
                #add nanoparticle to the left
                np_diameter = np.random.normal(self.dDiameter,self.np_d_std)
                while np_diameter < self.np_lb or np_diameter > self.np_ub:
                    np_diameter = np.random.normal(self.dDiameter,self.np_d_std)
                    
                #to ensure no overlap a priori, radius + jitter cannot exceed desired radius + ligand length
                upperl_bound = (self.np_ub - np_diameter)/2.0
                iterations = 0

                #random location in 3D space
                x_jit = np.random.normal(0,self.lSTD)
                y_jit = np.random.normal(0,self.lSTD)
                z_jit = np.random.normal(0,self.lSTD)
                total_jit = math.sqrt(self.x_jit**2 + self.y_jit**2 + self.z_jit**2)
                while total_jit >= upperl_bound:
                    x_jit = np.random.normal(0,self.lSTD)
                    y_jit = np.random.normal(0,self.lSTD)
                    z_jit = np.random.normal(0,self.lSTD)
                    total_jit = math.sqrt(self.x_jit**2 + self.y_jit**2 + self.z_jit**2)
                    iterations += 1
                    if iterations > 20:
                        x_jit = 0.0
                        y_jit = 0.0
                        z_jit = 0.0
                        break

                x_jit = 0.0
                y_jit = 0.0
                z_jit = 0.0
                b_x = -(self.a_1[0]*math.cos(phi)+self.a_1[2]*math.sin(phi))
                b_z = -(-self.a_1[0]*math.sin(phi)+self.a_1[2]*math.cos(phi))

#                 yLayer = self.calculateYLayer(j[1])
#                 x0 = round((j[0] - self.a_3[0]*(yLayer) - self.baseArray[0])/self.a_2[0])*self.a_2[0] + self.baseArray[0] + self.a_3[0]*yLayer
#                 y0 = (yLayer)*self.a_3[1] + self.baseArray[1]
#                 z0 = round((j[2] - self.a_3[2]*(yLayer) - self.baseArray[2])/self.a_2[2])*self.a_2[2] + self.baseArray[2] + self.a_3[2]*yLayer
                np_x = j[0] - j[4] + x_jit + b_x 
                np_y = j[1] - j[5]+ y_jit
                np_z = j[2] - j[6] + z_jit + b_z
                
                #append nanoparticle
                npList.append([np_x,np_y,np_z,np_diameter])
                    
        return npList, sorted_idList
            
    #the method to create grains inside the nanoparticle samples. 
    def createGrains(self):
        #initialize grain seed points
        grainArray = self.np_array
        grainIDArray = self.np_id_array
        self.defineGrains()
        
        #iterate through all grain seed points. j is the index of the nanoparticle inside np_array
        for index, j in enumerate(self.grainList):
            #the angle is defined outside the class, for repeatability
            theta = self.gAngles[index]
            
            #grain line
            #line equation: z-z1 = m1*(x-x1)
            i = self.np_array[j]
            x1 = i[0]
            y1 = i[1]
            z1 = i[2]
            m1 = -math.tan(theta)#(z/x)
            
            #right x-edge line
            x2 = self.baseArray[0]
            #z2 = self.cellz - self.baseArray[2]
            z2 = self.cellz_adjusted - self.baseArray[2]
            xEdgeAngle = self.gamma
                
            if xEdgeAngle != math.pi/2:
                m2 = 1/math.tan(xEdgeAngle)
            else:
                m2 = 0
            
            #left x-edge line
            x3 = self.baseArray[0]
            z3 = self.baseArray[2]
            m3 = m2
            
            
            #find lattice edge intercepts
            have_top = False
            have_bottom = False
            
            #checking if intersection of grain line with right x-edge line exists
            x_right = (z1-z2-m1*x1+m2*x2)/(m2-m1)
            if x_right >= self.baseArray[0] and x_right <= x1:
                #bottom intersects right x-edge
                z_b = m1*(x_right-x1)+z1
                x_b = x_right
                have_bottom = True
                
            elif x_right >= x1 and x_right <= self.cellx - 2*self.baseArray[0]:
                #top intersects right x-edge
                z_t = m1*(x_right-x1)+z1
                x_t = x_right
                have_top = True
             
            #checking if intersection with left x-edge line exists
            x_left = (z1-z3-m1*x1+m3*x3)/(m3-m1)
            if x_left >= x1 and x_left <= self.cellx-2*self.baseArray[0]:
                #top intersects left x-edge
                z_t = m1*(x_left-x1) + z1
                x_t = x_left
                have_top = True
                
            elif x_left <= x1 and x_left >= self.baseArray[0]:
                #bottom intersects left x-edge
                z_b = m1*(x_left-x1)+z1
                x_b = x_left
                have_bottom=True
            
            #check if intersection instead with bottom x-edge
            if not have_bottom:
                x_b = self.baseArray[0]
                z_b = m1*(x_b-x1)+z1
                
            #check if intersection instead with top x-edge
            if not have_top:
                x_t = self.cellx-self.baseArray[0]
                z_t = m1*(x_t-x1)+z1
                
            #now start to define how much big the rotated lattice needs to be to cover minimally the entire region
            line_length = math.sqrt((x_t-x_b)**2 + (z_t-z_b)**2)
            bottom_overhang = 0
            
            #store the lattice angles and lattice vectors, and generate the new ones for the 011 grain
            #the angles and vectors will be restored after creating the 011 lattice
            angles = [self.gamma, self.beta, self.alpha]
            latticeVectors = [self.a_1, self.a_2, self.a_3]
            self.setSpacing011Plane()
            
            #defining the rotation angle for the lattice. It will be a counterclockwise rotation
            phi = theta + math.pi/2 - self.gamma
            
            #Now we will calculate how much overhang is needed, and adjust bottom_overhang accordingly
            #top x-edge
            if not have_top and phi < 0:
                #need to have more overhang. Relevant lattice angle is math.pi - gamma
                z_distance = self.cellz_edge_adj + self.baseArray[2] - z_t
                crit_angle = 90*math.pi/180 + theta
                x_angle = math.pi - crit_angle - (math.pi - self.gamma) 
                if x_angle > 0:
                    x_overhang = math.sin(x_angle)*z_distance/math.sin(math.pi - self.gamma)
                else:
                    x_overhang = 0
                line_length += x_overhang
            
            #bottom x-edge
            if not have_bottom and phi > 0:
                #need to have more overhang. Relevant lattice angle is gamma
                z_distance = self.cellz_edge_adj + self.baseArray[2] - z_b
                crit_angle = math.pi/2 - theta
                x_angle = math.pi - crit_angle - self.gamma
                if x_angle > 0:
                    x_overhang = math.sin(x_angle)*z_distance/math.sin(self.gamma)
                else:
                    x_overhang = 0
                line_length += x_overhang
                bottom_overhang += x_overhang
               
            #left z-edge
            if have_top and phi > 0:
                #need to have more overhang
                #define the top left z location of the lattice
                topleftz = self.baseArray[2] + (self.xNum-1)*self.a_2[2]
                
                x_distance = math.sqrt((self.cellx_edge - self.baseArray[0] - x_t)**2 + (topleftz - z_t)**2)
                crit_angle = phi
                x_angle = self.gamma - crit_angle
                if x_angle > 0:
                    x_overhang = math.sin(x_angle)*x_distance/math.sin(math.pi-self.gamma)
                else:
                    x_overhang = 0
                line_length += x_overhang
            
            if have_bottom and phi < 0:
                #need to have more overhang
                x_distance = math.sqrt((x_b - self.baseArray[0])**2 + (z_b - self.baseArray[2])**2)
                crit_angle = theta + self.gamma - math.pi/2
                x_angle = math.pi - crit_angle - self.gamma
                if x_angle > 0:
                    x_overhang = math.sin(x_angle)*x_distance/math.sin(self.gamma)
                else:
                    x_overhang = 0
                line_length += x_overhang
                bottom_overhang += x_overhang
                
            #right z-edge
            if have_top and phi < 0:
                #need more overhang
                x_distance = math.sqrt((x_t - self.baseArray[0])**2 + (z_t - self.cellz_edge_adj)**2)
                crit_angle = math.pi - (theta + self.gamma - math.pi/2)
                x_angle = self.gamma - crit_angle
                if x_angle > 0:
                    x_overhang = math.sin(x_angle)*x_distance/math.sin(math.pi - self.gamma)
                else:
                    x_overhang = 0
                line_length += x_overhang

            
            #condition for bisection of the transport direction
            if not have_top and not have_bottom:
                bisect = 1
            else: 
                bisect = (x_t - x_b)/(self.cellx - self.dDiameter)
                
            #determine if new lattice will be large enough in the x-direction
            #Note: the 011 lattice cellx_edge will remain the same, as the lattice constant in the x-direction does not change
            if line_length <= self.cellx_edge:
                old_xNum = self.xNum
            else:
                #create new lattice of appropriate size
                old_xNum = self.xNum
                self.xNum = int(math.ceil((line_length + bottom_overhang - self.dDiameter)/self.sidelength + 2)) 
            
            #store the original cell values
            cellArray = [self.cellx,self.celly,self.cellz]
            baseArray = self.baseArray
            zNumAdjust = self.zNumAdjust
             
            #center the new array in the y-direction
            currYHeight = latticeVectors[2][1] #the stored a_3[1]
            yLayerCorrection = (currYHeight - self.a_3[1])/2
            self.baseArray[1] = self.baseArray[1] + yLayerCorrection
            
            #create the new array
            self.zNumAdjust = 0
            self.initiatePositions()
            new_npArray, new_idArray = self.generateLattice()
            new_npArray[:], new_idArray[:] = self.createVacancies(new_npArray, new_idArray)
            
            #add randomness to the left z-edge of the nanoparticle array (to be rotated)
            npArrayCopy, npArrayCopyIDs = self.add_zBoundary_randomness(new_npArray, new_idArray,self.xNum, 0)
            #npArrayCopy = new_npArray[:]
                
            rotate_length = math.sqrt((x1-x_b)**2 + (z1-z_b)**2)
            rotate_length += bottom_overhang
            
            if self.gamma < math.pi/2:
                x_rotate = self.baseArray[0] + rotate_length*math.sin(self.gamma)
                z_rotate = self.baseArray[2] + rotate_length*math.cos(self.gamma)
            else:
                x_rotate = self.baseArray[0] + rotate_length*math.cos(self.gamma-math.pi/2)
                z_rotate = self.baseArray[2] - rotate_length*math.sin(self.gamma-math.pi/2)
                     
            rotated_np_array = [ [i[0] - x_rotate, i[1], i[2] - z_rotate,i[3]] for i in npArrayCopy]
            
            #perform rotation in place
            rotated_np_array[:] = [[i[0]*math.cos(phi)+i[2]*math.sin(phi),i[1],-i[0]*math.sin(phi)+i[2]*math.cos(phi),i[3]] for i in rotated_np_array]
            rotated_np_array[:] = [ [i[0] + x_rotate, i[1], i[2] + z_rotate,i[3]] for i in rotated_np_array] 
                
            #move rotated array to correct position
            x_add = x1 - x_rotate
            z_add = (z1 - z_rotate)
            rotated_np_array[:] = [[i[0] + x_add, i[1], i[2] + z_add,i[3]] for i in rotated_np_array]
            
            #recover the original cell values, latticeVectors, and xNum
            [self.cellx,self.celly,self.cellz] = cellArray
            [self.a_1, self.a_2, self.a_3] = latticeVectors
            [self.gamma, self.beta, self.alpha] = angles
            self.baseArray = baseArray
            self.xNum = old_xNum
            self.zNumAdjust = zNumAdjust
            
            #delete points along and over line, and add rotated array
            lineValues = [z1,x1,m1]
            newgrainArray, grainRemovalArray, grainArrayIndices = self.removeNanoparticlesOnGrainLine(theta,grainArray, rotated_np_array, lineValues)
            
            rotatedA = np.asarray(rotated_np_array)
            grainA = np.asarray(newgrainArray)
            min_sep_tol = self.min_compression_thr_tol
            overlap_thresh = self.NP_grain_compression_overlap_thr
            compressed = True
            while compressed:
                compressed = False
                for ind,i in enumerate(grainA):
                    index = grainArrayIndices[ind]

                    cc = i[0:3] - rotatedA[:,:3]
                    dist = np.linalg.norm(cc,axis = 1)
                    overlap = np.subtract(i[3]/2 + rotatedA[:,3]/2, dist)

                    n_array = np.where(overlap > -min_sep_tol/4)[0]
                    if len(n_array) != 0:
                        if len(n_array) == 1:
                            n = n_array[0]
                            total_overlap = overlap[n]
                            total_cc = i[0:3] - rotatedA[n,:3]

                            if total_overlap > overlap_thresh:
                                grainRemovalArray.append(index)
                                #print("overlap too large so removed ",index)

                            #if overlap less than threshold, compress nanoparticles rather than remove them
                            else:
                                total_overlap += min_sep_tol
                                cc_hat = np.divide(total_cc,np.linalg.norm(total_cc))
                                overlap_vector = np.multiply(cc_hat,total_overlap)
                                #print("originally compressing: ",index)
                                grainA = self.compressNanoparticlesAfterRemoval(grainA, grainRemovalArray,grainArrayIndices, index, overlap_vector)
                                compressed = True
                        else:
                            grainRemovalArray.append(index)
                            #print("overlapped with two so removed ",index)
                        
            remainingGrainArray = [i for index,i in enumerate(grainA) if grainArrayIndices[index] not in grainRemovalArray]

            #update the final crystal
            totalArray = remainingGrainArray + rotated_np_array
            
            #remove points outside crystal bounds
            grainArray = self.removeNanoparticlesOutside(totalArray)
            self.boundaryInfo = (theta, bisect)
         
        #update the final crystal
        self.np_array = grainArray
        
    def compressNanoparticlesRaw(self, array, index, overlapVector, direction): #add direction? e.g. left or right
        array[index][0:3] = np.add(array[index][0:3],overlapVector)

        if direction == 'left':
            leftSide = (index%self.zNum == 0)
            bottomEdge = (index%(self.zNum*self.xNum) < self.zNum)
            topEdge = ((index + self.zNum)%(self.zNum*self.xNum) < self.zNum)

            if(leftSide and not bottomEdge and not topEdge):
                neighbors = [index + self.zNum, index-self.zNum]
            elif(leftSide and not bottomEdge):
                neighbors = [index-self.zNum]
            elif(leftSide and not topEdge):
                neighbors = [index + self.zNum]
            elif(not leftSide and topEdge):
                neighbors = [index-1,index-self.zNum-1, index-self.zNum]
            elif(not leftSide and bottomEdge):
                neighbors = [index+self.zNum, index+self.zNum-1, index-1]
            else:
                neighbors = [index+self.zNum, index+self.zNum-1, index-1,index-self.zNum-1, index-self.zNum]
        elif direction == 'right':
            rightSide = (index%(self.zNum) == self.zNum - 1)
            bottomEdge = (index%(self.zNum*self.xNum) < self.zNum)
            topEdge = ((index + self.zNum)%(self.zNum*self.xNum) < self.zNum)

            if(rightSide and not bottomEdge and not topEdge):
                neighbors = [index + self.zNum, index-self.zNum]
            elif(rightSide and not bottomEdge):
                neighbors = [index-self.zNum]
            elif(rightSide and not topEdge):
                neighbors = [index + self.zNum]
            elif(not rightSide and topEdge):
                neighbors = [index+1,index-self.zNum+1, index-self.zNum]
            elif(not rightSide and bottomEdge):
                neighbors = [index+self.zNum, index+self.zNum+1, index+1]
            else:
                neighbors = [index+self.zNum, index+self.zNum+1, index+1,index-self.zNum + 1, index-self.zNum]

        min_sep_tol = self.min_compression_thr_tol

        for i in neighbors:
            cc = np.subtract(array[i][0:3],array[index][0:3])
            dist = np.linalg.norm(cc)
            overlap = array[index][3]/2 + array[i][3]/2 - dist
            
            if overlap > -min_sep_tol/4:
                #if overlap < min_sep_tol: overlap += min_sep_tol
                overlap += min_sep_tol
                cc_hat = np.divide(cc,dist)
                overlap_vector = np.multiply(cc_hat,overlap)
                self.compressNanoparticlesRaw(array, i, overlap_vector, direction)

        return array
        
    def compressNanoparticlesAfterRemoval(self, array, removalArray, indexArray, index, overlapVector): #add direction? e.g. left or right
        n = indexArray.index(index)
        array[n][0:3] = np.add(array[n][0:3],overlapVector)

        leftSide = (index%self.zNum == 0)
        rightSide = (index%(self.zNum) == self.zNum - 1)
        bottomEdge = (index%(self.zNum*self.xNum) < self.zNum)
        topEdge = ((index + self.zNum)%(self.zNum*self.xNum) < self.zNum)
        bottomY = index < self.zNum*self.xNum
        topY = index >= self.zNum*self.xNum*(self.yNum - 1)

        #on the left side but not the top or bottom (in x)
        if(leftSide and not (bottomEdge or topEdge)):
            neighbors = [index + self.zNum, index-self.zNum, index + 1, index + self.zNum + 1, index - self.zNum + 1]
        #on the left side and on the top in x
        elif(leftSide and topEdge):
            neighbors = [index-self.zNum, index + 1, index - self.zNum + 1]
        #on the left side and on the bottom in x
        elif(leftSide and bottomEdge):
            neighbors = [index + self.zNum, index + 1, index + self.zNum + 1]
        elif(rightSide and not (bottomEdge or topEdge)):
            neighbors = [index + self.zNum, index-self.zNum, index - 1, index + self.zNum - 1, index - self.zNum - 1]
        #on the side and on the top in x
        elif(rightSide and topEdge):
            neighbors = [index-self.zNum, index - 1, index - self.zNum - 1]
        #on the side and on the bottom in x
        elif(rightSide and bottomEdge):
            neighbors = [index + self.zNum, index - 1, index + self.zNum - 1]
        #on the top in the middle
        elif(not (leftSide or rightSide) and topEdge):
            neighbors = [index-1,index+1, index-self.zNum-1, index-self.zNum, index - self.zNum + 1]
        #on the bottom in the middle
        elif(not (leftSide or rightSide) and bottomEdge):
            neighbors = [index+self.zNum, index+self.zNum-1, index + self.zNum + 1, index-1, index + 1]
        #in the middle
        else:
            neighbors = [index+self.zNum, index+self.zNum-1, index + self.zNum + 1, index-1, index + 1, index-self.zNum-1, index-self.zNum, index - self.zNum + 1]
        
        if bottomY:
            topCenter = index + self.zNum*self.xNum
            neighbors.extend([topCenter, topCenter + self.zNum, topCenter - self.zNum, topCenter + 1, topCenter - 1])
            if rightSide:
                neighbors.remove(topCenter + 1)
            if leftSide:
                neighbors.remove(topCenter - 1)
            if bottomEdge:
                neighbors.remove(topCenter - self.zNum)
            if topEdge:
                neighbors.remove(topCenter + self.zNum)
        elif topY:
            bottomCenter = index - self.zNum*self.xNum
            neighbors.extend([bottomCenter, bottomCenter + self.zNum, bottomCenter - self.zNum, bottomCenter -1, bottomCenter + 1])
            
            if rightSide:
                neighbors.remove(bottomCenter + 1)
            if leftSide:
                neighbors.remove(bottomCenter - 1)
            if bottomEdge:
                neighbors.remove(bottomCenter - self.zNum)
            if topEdge:
                neighbors.remove(bottomCenter + self.zNum)
        else:
            topCenter = index + self.zNum*self.xNum
            bottomCenter = index - self.zNum*self.xNum
            
            neighbors.extend([bottomCenter, bottomCenter + self.zNum, bottomCenter - self.zNum, bottomCenter + 1, bottomCenter - 1])
            neighbors.extend([topCenter, topCenter + self.zNum, topCenter - self.zNum, topCenter + 1, topCenter - 1])
            
            if rightSide:
                neighbors.remove(topCenter + 1)
                neighbors.remove(bottomCenter + 1)
            if leftSide:
                neighbors.remove(topCenter - 1)
                neighbors.remove(bottomCenter - 1)
            if bottomEdge:
                neighbors.remove(topCenter - self.zNum)
                neighbors.remove(bottomCenter - self.zNum)
            if topEdge:
                neighbors.remove(topCenter + self.zNum)
                neighbors.remove(bottomCenter + self.zNum)
            
        min_sep_tol = self.min_compression_thr_tol
        
        for j in neighbors:
            #print("neighbor: ",j)
            if not j in removalArray:
                i = indexArray.index(j)
                cc = np.subtract(array[i][0:3],array[n][0:3])
                dist = np.linalg.norm(cc)
                overlap = array[i][3]/2 + array[n][3]/2 - dist
                
                #if overlap > self.NP_grain_compression_overlap_thr:
                    #print("large overlap")
                    #print(array[indexArray.index(index)])
                    #print(array[indexArray.index(j)])
                
                if overlap > -min_sep_tol/4:
                    overlap += min_sep_tol
                    cc_hat = np.divide(cc,dist)
                    overlap_vector = np.multiply(cc_hat,overlap)
                    self.compressNanoparticlesAfterRemoval(array, removalArray,indexArray, j, overlap_vector)
        
        return array
        
    #the method to create twins inside the nanoparticle samples. 
    def createTwin(self):
        #initialize twin boundary seed points, recyling the grain boundary method
        self.twinList = []
        
        for x in range(0,self.tbNumber):
            #want to target the top right corner of the top layer
            self.twinList.append(int(len(self.np_array)-1))
            
            
        for index, j in enumerate(self.twinList):
            i = self.np_array[j]

            #twin line
            x0 = i[0]
            y0 = i[1]
            z0 = i[2]
            #m0 = -math.tan(theta)#(z/x)

            #right z-edge line
            x1 = self.baseArray[0]
            z1 = self.cellz - self.baseArray[2]
            xEdgeAngle = self.gamma

            if xEdgeAngle != math.pi/2:
                m1 = 1/math.tan(xEdgeAngle)
            else:
                m1 = 0

            #left z-edge line
            x2 = self.baseArray[0]
            z2 = self.baseArray[2]
            m2 = m1

            #create a lattice twice as large
            old_xNum = self.xNum
            old_zNum = self.zNum
            zNumAdjust = self.zNumAdjust
            
            self.xNum = 2*old_xNum
            self.zNum = 2*old_zNum
            self.zNumAdjust = 0

            #store the original cell values
            cellArray = [self.cellx,self.celly,self.cellz]
            
            #initiate lattice placement
            self.initiatePositions()

            #generate lattices
            twinArray, new_idArray = self.generateLattice()
            
            #define the twin nanoparticle to calculate theta with respect to
            #currIndex = [i[:3] for i in twinArray].index([x0, y0, z0])
            currIndex = int(self.xNum*self.zNum*(self.yNum - 0.5)) + int(self.zNum/2) #target middle of top layer
            
            #calculate theta and the twin line's slope
            theta, twinIndex = self.defineTwinAngle(twinArray, currIndex)
            m0 = -math.tan(theta)
            
            #create vacancies in the lattices
            twinArray[:], new_idArray[:] = self.createVacancies(twinArray, new_idArray)

            #recover the original cell values
            [self.cellx,self.celly,self.cellz] = cellArray
            
            #perform reflection via projection, e.g. Ref(v) = 2*dot(v,l)/dot(l,l)*l - v
            mirroredArray = []
            
            #define the twin boundary line vector
            lineV = np.asarray([1, m0])
            for i in twinArray:
                #make line seed the origin
                l = np.asarray([x0,z0])
                point = np.asarray([i[0],i[2]])
                
                projLine = l + lineV*np.dot(point - l,lineV)/np.dot(lineV,lineV)
                refV = (2*projLine - point).tolist()
                
                #make sure nanoparticles are different diameters, to preserve statistics
                np_jitter = math.sqrt(i[4]**2 + i[5]**2 + i[6]**2)
                new_ub = self.np_ub - 2*np_jitter
                np_diam = np.random.normal(self.dDiameter,self.np_d_std)
                while np_diam < self.np_lb or np_diam > new_ub:
                    np_diam = np.random.normal(self.dDiameter,self.np_d_std)
                nano = [refV[0],i[1],refV[1],np_diam]
                mirroredArray.append(nano)
            
            #shift the two arrays by 10 nanoparticles along the diagonal direction, so that the seed is 
            #approximately in the center of the sample
            x_add = -(old_xNum/2)*(self.a_2[0])
            z_add = -(old_zNum/2)*(self.a_1[2] + self.a_2[2])
            twinArray[:] = [[i[0] + x_add, i[1], i[2] + z_add,i[3]] for i in twinArray]
            mirroredArray[:] = [[i[0] + x_add, i[1], i[2] + z_add,i[3]] for i in mirroredArray]
            
            #adjust twin location to center
            x0 = x0 + x_add
            z0 = z0 + z_add
            
            ##################################
            #remove overlapping nanoparticles#
            ##################################
            
            #remove nanoparticles on and past the twin line, for both arrays
            lineValues = [z0,x0,m0]
            [newTwinArray, twinRemovalArray, twinArrayIndices] = self.removeNanoparticlesOnTwinLine(theta,twinArray,lineValues,'right')[0:]
            [newMirroredArray, mirrorRemovalArray, mirroredArrayIndices] = self.removeNanoparticlesOnTwinLine(theta,mirroredArray,lineValues,'left')[0:]
            
            mirrorA = np.asarray(newMirroredArray)
            twinA = np.asarray(newTwinArray)
            min_sep_tol = self.min_compression_thr_tol
            overlap_thresh = self.NP_twin_compression_overlap_thr
            compressed = True
                
            while compressed:
                compressed = False
                for ind,i in enumerate(twinA):
                    index = twinArrayIndices[ind]

                    cc = i[0:3] - mirrorA[:,:3]
                    dist = np.linalg.norm(cc,axis = 1)
                    overlap = np.subtract(i[3]/2 + mirrorA[:,3]/2, dist)

                    n_array = np.where(overlap > -min_sep_tol/4)[0]
                    if len(n_array) != 0 and index not in twinRemovalArray:
                        if len(n_array) == 1:
                            n = n_array[0]
                            total_overlap = overlap[n]
                            total_cc = i[0:3] - mirrorA[n,:3]

                            #if overlap greater than threshhold, remove NP
                            if total_overlap > overlap_thresh:
                                twinRemovalArray.append(index)

                            #if overlap less than threshold, compress nanoparticles rather than remove them
                            elif total_overlap > -min_sep_tol/4:
                                total_overlap += min_sep_tol
                                cc_hat = np.divide(total_cc,np.linalg.norm(total_cc))
                                overlap_vector = np.multiply(cc_hat,total_overlap)
                                twinA = self.compressNanoparticlesAfterRemoval(twinA, twinRemovalArray, twinArrayIndices, index, overlap_vector)
                                compressed = True
                        elif len(n_array) > 1:
                            twinRemovalArray.append(index)
                    
            remainingTwinArray = [i for index,i in enumerate(twinA) if twinArrayIndices[index] not in twinRemovalArray]
            
            #create the fnew crystal
            totalArray = newMirroredArray + remainingTwinArray
            
            #################################################################################
            #Rotate the lattice a random angle, so that twins can be in different directions#
            #################################################################################
            
            #rotate the crystal a random angle
            #phi = np.random.uniform(-math.pi/2, math.pi/2)
            phi = 0.0
            
            #adjust the lattice so that the twin line origin point is at the actual origin
            rotated_np_array = [ [i[0] + x_add, i[1], i[2] + z_add,i[3]] for i in totalArray]

            #perform rotation in place
            rotated_np_array[:] = [[i[0]*math.cos(phi)+i[2]*math.sin(phi),i[1],-i[0]*math.sin(phi)+i[2]*math.cos(phi),i[3]] for i in rotated_np_array]

            #move rotated array to be centered over old array
            rotated_np_array[:] = [ [i[0] - x_add, i[1], i[2] - z_add,i[3]] for i in rotated_np_array] 
            
            #redefine m0 now that the twin line has rotated an additional ϕ
            newAngle = theta + phi
            m0 = -math.tan(newAngle)
            
            ################################################
            #calculate bisection of the transport direction#
            ################################################
            
            #define what the cell-limits are for bisection purposes
            x_topEdge = self.cellx - self.baseArray[0]
            x_botEdge = self.baseArray[0]
            
            #calculate what the bottom x and top x intersections are
            if(m0 != m1):
                #intersection with right x-edge
                x_right = (z1 - m1*x1 - z0 + m0*x0)/(m0-m1)
                #intersection with left x-edge
                x_left = (z2 - m2*x2 - z0 + m0*x0)/(m0-m2)
                
                #set x_t (x-top) and x_b (x-bottom)
                if(x_right > x_left):
                    x_t = x_right
                    x_b = x_left
                else:
                    x_t = x_left
                    x_b = x_right

                if(x_t > x_topEdge):
                    x_t = x_topEdge
                if(x_b < x_botEdge):
                    x_b = x_botEdge
            else:
                x_t = x_topEdge
                x_b = x_botEdge
            
            bisect = (x_t - x_b)/(self.cellx - 2*self.baseArray[0])
            
            #############################################
            #Remove nanoparticles outside crystal bounds#
            #############################################
                            
            #trim edges
            self.np_array = self.removeNanoparticlesOutside(rotated_np_array)
            self.boundaryInfo = (theta, bisect, phi*180/math.pi)
            
            #recover the original xNum and zNum
            self.xNum = old_xNum
            self.zNum = old_zNum
            self.zNumAdjust = zNumAdjust
            
    def createCrack(self):
        bot_orgInd = int(len(self.np_array)/4 - self.zNum/2) #middle of bottom layer
        top_orgInd = int(bot_orgInd + len(self.np_array)/2)
        
        #crackLength = 11 #number nanoparticles
        crackWidth = 6.6 #nm
        min_sep_tol = self.min_compression_thr_tol
        
        npArray = np.asarray(self.np_array)
        
        orgXLayer = bot_orgInd//self.zNum
        x_ideal = int((self.crack_length-1)/2)
        x_upper = x_ideal if x_ideal <= self.zNum - orgXLayer - 1 else self.zNum - orgXLayer - 1
        cM_upper = x_upper if x_ideal <= self.zNum - orgXLayer - 1 else x_upper + 1
        x_lower = x_ideal if x_ideal <= orgXLayer else orgXLayer
        cM_lower = x_lower if x_ideal <= orgXLayer else x_lower + 1
    

        for i in it.chain(range(bot_orgInd-x_lower*self.zNum,bot_orgInd + (x_upper)*self.zNum,self.zNum),range(top_orgInd-x_lower*self.zNum,top_orgInd + (x_upper)*self.zNum,self.zNum)):
            crackModifier = int((i-bot_orgInd)/self.zNum if i < len(self.np_array)/2 else (i-top_orgInd)/self.zNum)
            if(crackModifier >= 0):
                currCrackWidth = crackWidth - crackModifier//cM_upper*(0.5*crackWidth)
            else:
                currCrackWidth = crackWidth - abs(crackModifier)//cM_lower*(0.5*crackWidth)

            #apply crack in the x-direction
            #leftIndex = originIndex + i*self.zNum
            leftIndex = i
            npArray[leftIndex][2] -= currCrackWidth/2
            leftNP = npArray[leftIndex]

            #rightIndex = originIndex + i*self.zNum + 1
            rightIndex = i +1
            npArray[rightIndex][2] += currCrackWidth/2
            rightNP = npArray[rightIndex]

            #left
            #cc = leftNP[0:3] - npArray[:,:3]
            leftArray = np.concatenate((npArray[:leftIndex],npArray[(leftIndex+1):]))
            cc = leftArray[:,:3] - leftNP[0:3]
            dist = np.linalg.norm(cc,axis = 1)
            overlap = np.subtract(leftNP[3]/2 + leftArray[:,3]/2, dist)

            n = np.argmax(overlap)
            total_overlap = overlap[n]
            total_cc = leftArray[n,:3] - leftNP[0:3]

            oIndex = n
            if(oIndex >= leftIndex):
                oIndex += 1

            if total_overlap > -min_sep_tol/4:
                total_overlap += min_sep_tol
                cc_hat = np.divide(total_cc,np.linalg.norm(total_cc))
                overlap_vector = np.multiply(cc_hat,total_overlap)
                npArray[:] = self.compressNanoparticlesRaw(npArray,oIndex, overlap_vector,'left')

            #right
            rightArray = np.concatenate((npArray[:rightIndex],npArray[(rightIndex+1):]))
            cc = rightArray[:,:3] - rightNP[0:3]
            dist = np.linalg.norm(cc,axis = 1)
            overlap = np.subtract(rightNP[3]/2 + rightArray[:,3]/2, dist)

            n = np.argmax(overlap)
            total_overlap = overlap[n]
            total_cc = rightArray[n,:3] - rightNP[0:3]

            oIndex = n
            if(oIndex >= rightIndex):
                oIndex += 1

            if total_overlap > -min_sep_tol/4:
                total_overlap += min_sep_tol
                cc_hat = np.divide(total_cc,np.linalg.norm(total_cc))
                overlap_vector = np.multiply(cc_hat,total_overlap)
                npArray[:] = self.compressNanoparticlesRaw(npArray,oIndex, overlap_vector,'right')
        
        newArray = self.removeNanoparticlesOutside(npArray)
        self.np_array = newArray
        
    def rotateLattice(self):
        #define rotation angle as a random angle between -90 and +90 degrees
        phi = np.random.uniform(-math.pi/2,math.pi/2)
        
        #create a lattice twice as large
        old_xNum = self.xNum
        old_zNum = self.zNum
        self.xNum = 2*old_xNum
        self.zNum = 2*old_zNum

        #store the original cell values
        cellArray = [self.cellx,self.celly,self.cellz]

        #initiate lattice placement
        self.initiatePositions()

        #generate lattices
        largerArray, new_idArray = self.generateLattice()
        
        centerIndex = int(self.zNum*(self.xNum/2) + 0.5*self.zNum)
        centerNP = largerArray[centerIndex]
        x0 = centerNP[0]
        z0 = centerNP[2]
        
        #create vacancies in the lattices
        largerArray[:], new_idArray[:] = self.createVacancies(largerArray, new_idArray)

        #recover the original cell values, xNum and zNum
        [self.cellx,self.celly,self.cellz] = cellArray
        self.xNum = old_xNum
        self.zNum = old_zNum
        
        #rotate by it's center
        rotated_np_array = [ [i[0] - x0, i[1], i[2] - z0,i[3]] for i in largerArray]
            
        #perform rotation in place
        rotated_np_array[:] = [[i[0]*math.cos(phi)+i[2]*math.sin(phi),i[1],-i[0]*math.sin(phi)+i[2]*math.cos(phi),i[3]] for i in rotated_np_array]
        
        #move rotated array to be centered over old array
        rotated_np_array[:] = [ [i[0] + x0/2, i[1], i[2] + z0/2,i[3]] for i in rotated_np_array] 
        
        self.np_array = self.removeNanoparticlesOutside(rotated_np_array)
        self.boundaryInfo = (phi, 0)
                                    
    def returnLattice(self):
        self.setSpacingLattice()
        self.initiatePositions()
        self.np_array, self.np_id_array = self.generateLattice()
        self.np_array[:], self.np_id_array[:] = self.createVacancies(self.np_array, self.np_id_array)
        self.cellV = self.calculateVolume()
        self.cellsize_array = [self.cellx, self.celly, self.cellz]
        self.createGrains()
        self.createTwin()
        
        if self.crack:
            self.createCrack()
        if self.rotate:
            self.rotateLattice()
        
        #create the flipped array
        if self.lType=="cubic" or self.lType == "cpc" or self.trim == "cubic":
            self.np_array_flipped = [[i[0], i[1], self.cellsize_array[2]-i[2], i[3]] for i in self.np_array] #flipped in the transport direction
        elif self.lType=="tric" or self.trim == 'tric':
            #right x-edge line
            x0 = self.baseArray[0]
            z0 = self.cellz #- self.baseArray[2]
            m = 1/math.tan(self.gamma)

            #z_right = z0 + m*(i[0]-x0), add an extra m*(i[0]-x0) so that it is still triclinic
            #self.np_array_flipped = [[i[0], i[1], z0 + 2*m*(i[0]-x0)-i[2], i[3]] for i in self.np_array] #flipped
            #note: this will reverse the orientation of the lattice
            self.np_array_flipped = [[i[0], i[1], z0 - i[2], i[3]] for i in self.np_array] #flipped in the transport direction
        
        return self.np_array, self.np_array_flipped, self.cellsize_array,self.cellV, self.boundaryInfo
        
        
class cubic_lattice(triclinic_lattice):
    #def _init_(self,stdMethod,sizeRange,currDiameter,paramDict):
    #    super()._init_(stdMethod,sizeRange,currDiameter,paramDict)
    def setLatticeAngles(self):
        #define lattice angles
        self.alpha = 90*math.pi/180
        self.beta = 90*math.pi/180
        self.gamma = 90*math.pi/180
    
class cpc_lattice(triclinic_lattice):
    #def _init_(self,stdMethod,sizeRange,currDiameter,paramDict):
    #    super()._init_(stdMethod,sizeRange,currDiameter,paramDict)
    def setLatticeAngles(self):
        #define lattice angles
        self.alpha = 90*math.pi/180
        self.beta = 90*math.pi/180
        self.gamma = 90*math.pi/180
        
    def setLatticeSpacing(self):
        self.x_spacing = self.dDiameter + 2*self.average_ll
        self.y_spacing = layer_spacing = (dDiameter + 2*average_ll)*math.sqrt(3)/2
        self.z_spacing = self.dDiameter + 2*self.average_ll
        
    def generateLattice(self):
        base_id = -1
        id = 0
        for y in range(1,self.yNum+1):
            for x in range(1,self.xNum+1):
                base_id += 1
                id = base_id
                for z in range(1,self.zNum+1):
                    npArray = self.appendNP(npArray)
                    idArray.append(id)
                    id+=self.xNum
                    if y == 2 and x == 1 and z == self.zNum:
                        cellz = curr_z - cellz + dDiameter + self.z_jit
                        self.cellz_edge = (self.zNum-1)*(self.sidelength) + self.dDiameter + self.z_jit
                    curr_z += self.z_spacing
                    
                self.curr_z = self.baseArray[2]
                if y == 2 and x == self.xNum:
                    self.cellx = self.curr_x - self.cellx + self.dDiameter + self.x_jit
                    self.cellx_edge = (self.xNum-1)*(self.sidelength) + self.dDiameter + self.x_jit
                curr_x += bottom_spacing
            if y%2 != 0:
                self.curr_x = self.baseArray[0] + self.x_spacing/2.0
                self.curr_z = self.baseArray[2] + self.z_spacing/2.0
            else:
                self.curr_x = self.baseArray[0]
                self.curr_z = self.baseArray[2]
                
            self.curr_y += self.y_spacing
            
        self.celly = self.curr_y - self.y_spacing - self.celly + self.dDiameter + self.y_jit  
        self.celly_edge = (self.yNum-1)*(self.sidelength) + -self.celly + self.dDiameter + self.y_jit
        return npArray,idArray

