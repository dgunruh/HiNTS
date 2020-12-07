#!/usr/bin/env python
# coding: utf-8
import numpy as np
import math
import os
from IPython.display import clear_output
import timeit
import itertools as it
import create_NP_superlattice as makeLattice

class disordered_NP_simulator():
    def createSample(self, method,currDiameter, paramDict):
        if method == "cubic":
            sample = makeLattice.cubic_lattice(method,currDiameter, paramDict)
        elif method == "cpc":
            sample = makeLattice.cpc_lattice(method,currDiameter, paramDict)
        elif method == "tric":
            sample = makeLattice.triclinic_lattice(method, currDiameter, paramDict)

        return sample.returnLattice()
    
    def generateSamples(self, method,currDiameter, paramDict):
        #generate samples, with inverted pairs
        sampleArray = []
        angleArray = []
        bisectArray = []
        numberSamples = paramDict['numberSamples']
        length = paramDict['length']
        elDensity = paramDict['electronDensity']
        s = 0
        total_cellV = 0
        while s < (numberSamples/2):
            print("\r" + "Creating sample pair %03d and %03d"%(2*s,2*s+1),end="")
            sArray, cellV, grainAngles,bisect = self.grainAngleRepetition(method, currDiameter, paramDict)
            sampleArray.extend(sArray)
            angleArray.append(grainAngles)
            bisectArray.append(bisect)
            total_cellV += cellV
            s += 1

        print("")

        #print out electron numbers for each nanoparticle, will be slightly different for each sample, so average over volume
        average_cellV = 2*total_cellV/(length)
        elNumber = elDensity*average_cellV
        print("NP Diameter: %f,  Number of electrons: %d" %(currDiameter, elNumber))

        return sampleArray, angleArray, bisectArray

    def outputNanoparticleSamples(self, filePath, length, sampleArray):
        i = 0
        while i < (length):
            filename = "nanoparticles" + str(i) + ".inp"
            fullpath = os.path.join(filePath, filename)

            np_array = sampleArray[i][0]
            cellsize_array = sampleArray[i][1]

            with open(fullpath,'w+') as f:
                f.write("BANDS\n")
                f.write("8 8\n") #degeneracy , always 8 8 for PbSe
                f.write("NANOPARTICLES\n")
                f.write("cell(1), %f\n" %cellsize_array[0])#X  center-center distance leftmost-rightmost + desired diameter
                f.write("cell(2), %f\n" %cellsize_array[1])#Y  Defines height (layers)
                f.write("cell(3), %f\n" %cellsize_array[2])#Z  Transport direction

                for j in range(0,len(np_array)):
                    x = str(np_array[j][0])
                    y = str(np_array[j][1])
                    z = str(np_array[j][2])
                    d = str(np_array[j][3])
                    f.write(x + ", " + y + ", " + z + ", " + d + "\n")

            i += 1

    def outputSampleAngles(self, filePath, angleArray, bisectArray):
        filename = "SampleAngles.txt"
        fullpath = os.path.join(filePath, filename)

        #merge angleArray and bisectArray
        for n,i in enumerate(angleArray):
            if isinstance(bisectArray[n],list):
                i.extend(bisectArray[n])
            else:
                i.append(bisectArray[n])

        with open(fullpath,'w+') as f:
            for j in range(1,len(angleArray)+1):
                #angleArray[j-1].append(bisectArray[j-1])
                string_of_angles = str(angleArray[j-1][0])
                for i in range(1,len(angleArray[j-1])):
                    string_of_angles += " "
                    string_of_angles += str(angleArray[j-1][i])    
                f.write(str(j) + " " + string_of_angles + "\n")
                
    def grainAngleRepetition(self, method, currDiameter, paramDict):
        gbNumber = paramDict['grainBoundaryNumber']
        tbNumber = paramDict['twinBoundaryNumber']
        angleRep = paramDict['angleRep']
        bisectOnly = paramDict['bisection']

        #define grain angle as a random angle between -80 and 80 degrees. This will be a rotation counterclockwise
        grainAngles = []
        i = 0
        if gbNumber > 0:
            while i < gbNumber:
                if not bisectOnly:
                    grainAngles.append(np.random.uniform(-80*math.pi/180,80*math.pi/180))
                else:
                    grainAngles.append(np.random.uniform(-40*math.pi/180, 50*math.pi/180))
                i+=1  
        paramDict['grainAngles'] = grainAngles

        sampleArray = []
        bisectArray = []
        angleArray = []
        total_cellV = 0.0
        a = 0
        while a < angleRep:
            np_array, np_array_flipped, cellsize_array, cellV, boundaryInfo= self.createSample(method, currDiameter, paramDict)
            if tbNumber > 0:
                (angle, bisect, phi) = boundaryInfo
                bisectArray.append([bisect, phi])
            else:
                (angle, bisect) = boundaryInfo
                bisectArray.append(bisect)
            sampleArray.append([np_array, cellsize_array])
            sampleArray.append([np_array_flipped, cellsize_array])
            angleArray.append(angle*180/math.pi)
            total_cellV += cellV

            a+=1

        return sampleArray, total_cellV, angleArray, bisectArray
    
    def getParams(self, stdMethod, sizeRange, cutMethod):
        paramDict = {}

        xNum = int(input("Number x nanoparticles: "))
        yNum = int(input("Number y nanoparticles: "))
        zNum = int(input("Number z nanoparticles: "))

        if sizeRange == True:
            dString = input("Desired Diameter range (in nm), use ' ' to separate: ")
            dDiameter = [float(i) for i in dString.split()]
            elDensity = float(input("Density of electrons (#/nm^3): "))
        else:
            dDiameter = float(input("Desired Diameter (in nm): "))
            elDensity = float(input("Density of electrons (#/nm^3): "))

        if stdMethod == "percent":
            dSTD = float(input("Percent Nanoparticle Size Variation (std): "))
        elif stdMethod == "fixed":
            dSTD = float(input("Nanoparticle Size Standard Deviation (in nm): "))

        lSTD = float(input("Nanoparticle Location Standard Deviation (in nm): "))
        vDensity = float(input("Percent Vacancies (in the QDxV(1-x) style): "))
        defectType = input("Type of Defect to create (G = grain, T = twin, C = crack, R = rotated, 110 = 110 lattice, N = none): ")

        gbNumber = 0
        tbNumber = 0
        crack = False
        rotate = False
        oneonezero = False

        if defectType == "G":
            gbNumber = int(input("Number of Grain Boundaries?: "))
        elif defectType == "T":
            tbNumber = int(input("Number of Twin Boundaries?: "))
        elif defectType == "C":
            crack = True
        elif defectType == "R":
            rotate = True
        elif defectType == "110":
            oneonezero = True

        if gbNumber > 0:
            paramDict['voidCutoff'] = voidCut = float(input("What percentage of samples should have a void at the grain boundary?: "))
        else:
            paramDict['voidCutoff'] = voidCut = 0

        if gbNumber > 0 or tbNumber > 0:
            paramDict['angleRep'] = int(input("How many samples of the same angle?: "))
        else:
            paramDict['angleRep'] = 1

        npnumber = xNum*yNum*zNum

        paramDict['defectType'] = defectType
        paramDict['crackLength'] = 7
        paramDict['npnumber'] = npnumber
        paramDict['xNum'] = xNum
        paramDict['yNum'] = yNum
        paramDict['zNum'] = zNum
        paramDict['dDiameter'] = dDiameter
        paramDict['dSTD'] = dSTD
        paramDict['stdMethod'] = stdMethod
        paramDict['lSTD'] = lSTD
        paramDict['electronDensity'] = elDensity
        paramDict['vacancyDensity'] = vDensity
        paramDict['grainBoundaryNumber'] = gbNumber
        paramDict['twinBoundaryNumber'] = tbNumber
        paramDict['crackBool'] = crack
        paramDict['rotateBool'] = rotate
        paramDict['110Grain'] = oneonezero
        paramDict['bisection'] = True #will the grains bisect only or will they be random
        paramDict['grainAngles'] = [np.random.uniform(-80*math.pi/180,80*math.pi/180)] #the default grain angle
        paramDict['trimmingMethod'] = cutMethod
        paramDict['numberSamples'] = int(input("How many samples (must be even number)?: "))  #number of samples (not including angle rep)
        paramDict['length'] = paramDict['numberSamples']*paramDict['angleRep'] #total number of samples including repetition
        
        return paramDict
    
    def run(self):
        stdMethod = "percent"
        #"percent": percent of mean
        #"fixed": user set value

        sizeRange = False
        #True: will generate nanoparticles samples in a size range (for diameter)
        #False: will generate nanoparticles at only one diameter

        method = "cubic"
        #cubic = simple cubic
        #cpc = close packed cubic (body centered)
        #tric = triclinic, angles are defined inside the triclinic sample generation code

        cutmethod = 'cubic'
        #cubic: vertical sides
        #tric: original lattice

        #get user inputs regarding parameters, and read into dictionary
        paramDict = self.getParams(stdMethod, sizeRange, cutmethod)
        #numberSamples = int(input("How many samples (must be even number)?: "))  #number of samples (not including angle rep)
        
        #use user built dictionary to define quantities
        angleRep = paramDict['angleRep'] #number of times grain angle will be repeated between samples
        length = paramDict['length'] #total samples (including angle rep)
        npnumber = paramDict['npnumber'] #number of nanoparticles
        dDiameter = paramDict['dDiameter'] #nanoparticle diameter
        defectType = paramDict['defectType']
        numberSamples = paramDict['numberSamples']

        start = timeit.default_timer()

        #choose to create samples for either a single diameter, or a range of diameters
        if sizeRange: 
            currDiameter = dDiameter[0]
            while currDiameter <= dDiameter[1]:

                sampleArray, angleArray, bisectArray = self.generateSamples(method, currDiameter, paramDict)

                #assign file path
                filePath = 'c:/Users/Davis/Research/data/nanoparticles/' + str(npnumber) + '_' +                    str(currDiameter) + 'nm' + '_' + str(numberSamples) +str(defectType)+'/'
                if not os.path.exists(filePath):
                    os.makedirs(filePath)
                #################
                # output angles #
                #################
                self.outputSampleAngles(filePath, angleArray,bisectArray)

                ##################
                # output samples #
                ##################
                self.outputNanoparticleSamples(filePath, length, sampleArray)

                #np diameter increment
                currDiameter += 0.5

        else:
            filePath = 'data/nanoparticles/' + str(npnumber) + '_' + str(dDiameter) + 'nm'+ '_' +                str(numberSamples) +str(defectType)+'/'
            if not os.path.exists(filePath):
                os.makedirs(filePath)
                
            sampleArray, angleArray, bisectArray = self.generateSamples(method, dDiameter, paramDict)

            #################
            # output angles #
            #################
            self.outputSampleAngles(filePath, angleArray,bisectArray)

            ##################
            # output samples #
            ##################

            self.outputNanoparticleSamples(filePath, length, sampleArray)

        stop = timeit.default_timer()

        print('Time: ', stop - start)
