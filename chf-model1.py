# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 09:04:00 2019

@author: jimmyli2002
"""

import sys
import re
import math

"""
HYPERPARAMETERS
"""

kClosestAtoms = 10
distanceTestRadius = 100
atomicDistanceTolerance = 0.1 #Default: 0.01 for extremely exact match
atomicAngleTolerance = 0.1

"""
SECTION
Name: Helper functions
Descritpion: A collection of functions which can be used to prepare the data
"""

"""
Name: pdbToData
Description: Takes a pdb file as input and converts it to a usable form
Input: input file path, output file path
"""
#match_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s{2}(\S{1,3})\s+([A-Z]{1,4})\s([\s\w])\s+(\d+)\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?")
match_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s+(\S{1,4})\s+([A-Z]{1,4})\s([\s\w])\s+(\d+)\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?")

def pdbToData(in_file_path, out_file_path):
    
    toWriteString = ""
    
    with open(in_file_path, "r") as f:
    
        for line in f.read().splitlines():
            if line.startswith("ENDMDL"):
                break
            match_list=match_pattern.findall(line)
            if match_list:
                final_str = ""
                for i in range(9):
                    if i != 8:
                        final_str += match_list[0][i] + "|"
                    else:
                        final_str += match_list[0][i] + "\n"
                toWriteString += final_str
            elif line.startswith("ATOM"):
                print("! " + line)
                
    with open(out_file_path, "w") as f:
        
        f.write(toWriteString)

"""
SECTION
Name: Get Template
Description: Gets the information of the template
"""

"""
FUNCTION
Name: get3DDistance
Description: Get the distance between two 3D coordinates
"""

def get3DDistance(x1, y1, z1, x2, y2, z2):
    
    return math.sqrt(1.0 * math.pow(float(x2) - float(x1), 2) + math.pow(float(y2) - float(y1), 2) + math.pow(float(z2) - float(z1), 2))

"""
FUNCTION
Name: findCenter
Description: Find the average point of all input locations
"""

def findCenter(input_locations):
    
    n = 0
    average = 0
    
    for location in input_locations:
        n += 1
        average += float(location)
        
    return round(1.0 * average / n, 3)

"""
FUNCTION
Name: getAngle
Description: Get angle between two 3D vectors
Input: v1 (x, y, z), v2 (x, y, z)
"""

def getAngle(v1, v2):
    
    dotProduct = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(math.pow(v1[0],2) + math.pow(v1[1],2) + math.pow(v1[2],2))
    mag2 = math.sqrt(math.pow(v2[0],2) + math.pow(v2[1],2) + math.pow(v2[2],2))
    
    return math.acos(dotProduct/(mag1*mag2))

"""
FUNCTION
Name: getTemplate
Description: Returns a sorted tuple with (center tuple(center atomic name, x, y, z), list of tuples with (atomic name, distance from center, angle between previous and itself))
Input: input file path, [(res_name, chain_id, res_seq)], splitting character
"""

def getTemplate(file_input, list_res, split_char):
    
    atomsToUse = [] #Format: [(atomic name, x, y, z)]
    xLocs = [] #Format: [float]
    yLocs = [] #Format: [float]
    zLocs = [] #Format: [float]

    #Information about the numerical center    
    centerX = 0 #Guarenteed to be filled
    centerY = 0 #Guarenteed to be filled
    centerZ = 0 #Guarenteed to be filled
    
    #Information about the atomic center atom
    centerAtomicName = ""
    centerAtomicDist = sys.maxsize
    centerAtomicX = 0
    centerAtomicY = 0
    centerAtomicZ = 0
    
    atomDistListWithVectors = [] #Format: [(atomic name, distance from center atom, vectorX, vectorY, vectorZ)]
    atomDistListWithAngleToPrev = [] #Format: [(atomic name, distance from center atom, angle to previous {center atom is 0}]
    
    #Fill atomsToUse, xLocs, yLocs, zLocs
    
    with open(file_input, "r") as f:
        
        for line in f:
            
            splitLine = line.split(split_char)
            
            resName = splitLine[1]
            chainId = splitLine[2]
            resSeq = int(splitLine[3])
            
            for residue in list_res:
                
                if resName == residue[0] and chainId == residue[1] and resSeq == residue[2]:
                    
                        #The line we are currently reading matches the requirements for any one of the residues
                        
                        atomTupple = (splitLine[0][0], splitLine[4], splitLine[5], splitLine[6])
                        atomsToUse.append(atomTupple)
                        
                        xLocs.append(float(splitLine[4]))
                        yLocs.append(float(splitLine[4]))
                        zLocs.append(float(splitLine[4]))
                        
    #atomsToUse, xLocs, yLocs, zLocs is now filled properly
    
    #Find numerical center information
    
    centerX = findCenter(xLocs)
    centerY = findCenter(yLocs)
    centerZ = findCenter(zLocs)
    
    #Find center atom information
    
    for atom in atomsToUse:
        
        dist = get3DDistance(centerX, centerY, centerZ, atom[1], atom[2], atom[3])
        
        if dist < centerAtomicDist:

            centerAtomicName = atom[0]
            centerAtomicDist = dist
            centerAtomicX = atom[1]
            centerAtomicY = atom[2]
            centerAtomicZ = atom[3]
            
    #All center atom information is now found
    
    #Find distance information between center and all atoms
    
    for atom in atomsToUse:
        
        dist = get3DDistance(centerAtomicX, centerAtomicY, centerAtomicZ, atom[1], atom[2], atom[3])
        
        atomDistListWithVectors.append((atom[0],
                                        dist,
                                        float(atom[1]) - float(centerAtomicX),
                                        float(atom[2]) - float(centerAtomicY),
                                        float(atom[3]) - float(centerAtomicZ)))
        
    #Sort all the atoms and calculate angle between them
    
    atomDistListWithVectors.sort(key=lambda tup: tup[1])
    
    atomDistListWithAngleToPrev.append((atomDistListWithVectors[0][0], atomDistListWithVectors[0][1], 0.0))
    atomDistListWithAngleToPrev.append((atomDistListWithVectors[1][0], atomDistListWithVectors[1][1], 0.0))
    
    for i in range(2, len(atomDistListWithVectors)):
        
        xPrev = atomDistListWithVectors[i - 1][2]
        yPrev = atomDistListWithVectors[i - 1][3]
        zPrev = atomDistListWithVectors[i - 1][4]
        
        xCurr = atomDistListWithVectors[i][2]
        yCurr = atomDistListWithVectors[i][3]
        zCurr = atomDistListWithVectors[i][4]
        
        angleBetween = getAngle((xPrev, yPrev, zPrev), (xCurr, yCurr, zCurr))
        
        atomDistListWithAngleToPrev.append((atomDistListWithVectors[i][0], atomDistListWithVectors[i][1], angleBetween))
        
    return ((centerAtomicName, centerAtomicX, centerAtomicY, centerAtomicZ), atomDistListWithAngleToPrev)

"""
SECTION
Name: Generating Interesting Regions Through Matching
"""

"""
FUNCTION
Name: getProtein
Description: Gets all pertinent information of a protein and returns it
Input: input file path, splitting character
"""

def getProtein(file_input, split_char):
    
    returnList = []
    
    with open(file_input, "r") as f:
        
        for line in f:
            
            splitLine = line.split(split_char)
            
            atomicName = splitLine[0][0]
            x = splitLine[4]
            y = splitLine[5]
            z = splitLine[6]
            
            returnList.append((atomicName, x, y, z))
            
    return returnList

"""
FUNCTION
Name: getImportantCenterAtomsNamesAndDist
Description: Returns all center atoms which may be a match using names and distance and the confidence of each
Input: protein [(aName, x, y, z)], template (center atom tupple, [atom tupple])
"""

def getImportantCenterAtomsNamesAndDist(protein, template):
    
    possibleCenterAtoms = [] #Format: [((atomic name, x, y, z), confidence)]
    counter = 0
    
    for atom in protein:
        
        counter += 1
        thresholdTracker = 0
        
        if counter % 1000 == 0 or counter == len(protein):
            
            print("Atom: [" + str(counter) + "/" + str(len(protein)) + "]")
        
        if atom[0] == template[0][0]:
            
            #The atomic name of the atom and the center of the template match, continue search
            
            #Get the distance to all other atoms and store them in a list
            
            atomDist = [] #Format: [(atomic name, distance)]
            
            for atom2 in protein:
                
                dist = get3DDistance(atom[1], atom[2], atom[3], atom2[1], atom2[2], atom2[3])
                
                atomDist.append((atom2[0], dist))
            
            atomDist.sort(key=lambda tup: tup[1])
            #DEBUG if(atom == ("N", "-1.074", "4.126", "57.539")): print(atomDist[:100])
            
            #atomDist is filled and sorted with distance from possible atom
            
            proteinPosition = 0
            
            for i in range(kClosestAtoms):
                
                #Check for 10 atoms
                
                requiredNextDistance = template[1][i][1]
                insideLoop = False
                
                if abs(atomDist[proteinPosition][1] - requiredNextDistance) < atomicDistanceTolerance and atomDist[proteinPosition][0] == template[1][i][0]:
                    
                    thresholdTracker += abs(atomDist[proteinPosition][1] - requiredNextDistance) 
                    insideLoop = True
                    
                if not insideLoop: #If the above check was not met, check the next j atoms
                
                    for j in range(distanceTestRadius):
                        
                        #Test the j next atoms in the list to see if distance match
                        
                        proteinPosition += 1
                        #TODO Add difference to threshold tracker if successful and return probabiliy with center
                        if abs(atomDist[proteinPosition][1] - requiredNextDistance) < atomicDistanceTolerance and atomDist[proteinPosition][0] == template[1][i][0]:
                            
                            #Increase the value to use to calculate threshold
                            
                            thresholdTracker += abs(atomDist[proteinPosition][1] - requiredNextDistance)
                            
                            #The atom at proteinPosition matches in distance and atom name
                            
                            insideLoop = True
                            break
                        
                if not insideLoop: break #If inside loop is false, then this atom is not possible
                
                #DEBUG if(atom == ("N", "-1.074", "4.126", "57.539")): print(i)
                
                if i == (kClosestAtoms - 1) and insideLoop: possibleCenterAtoms.append((atom, round(1 - (thresholdTracker / atomicDistanceTolerance / kClosestAtoms), 3)))
                #If the (i-1)th atom is being checked and the resultis that inside loop is true, then add this to possible along with confidence
                
    return possibleCenterAtoms

"""
USAGE AREA
"""

#Generate data file for protein 1aay
"""
pdbToData("C:\\proj\\proteinModelling\\motifSearching\\1aay.pdb", "C:\\proj\\proteinModelling\\motifSearching\\1aay_formatted.data")
"""

#Example with 1aay

#Example 1 residue for 1aay
#list_res = [("CYS", "A", 165), ("CYS", "A", 168), ("HIS", "A", 181), ("HIS", "A", 185)]
#Example 2 residue for 1aay
#list_res = [("CYS", "A", 107), ("CYS", "A", 112), ("HIS", "A", 125), ("HIS", "A", 129)]                        
#Example 3 residue for 1aay
#list_res = [("CYS", "A", 137), ("CYS", "A", 140), ("HIS", "A", 153), ("HIS", "A", 157)]                        

#Generate template for the third zinc finger of protein 1aay

#template = getTemplate("C:\\proj\\proteinModelling\\motifSearching\\1aay_formatted.data", list_res, "|")

#Generate important center atoms using anem and distance check

#protein = getProtein("C:\\proj\\proteinModelling\\motifSearching\\1aay_formatted.data", "|")

#importantCenters = getImportantCenterAtomsNamesAndDist(protein, template)

#Example with 4aiy

#pdbToData("C:\\proj\\proteinModelling\\motifSearching\\4aiy.pdb", "C:\\proj\\proteinModelling\\motifSearching\\4aiy_formatted.data")
list_res = [("HIS", "B", 10), ("HIS", "F", 10), ("HIS", "J", 10)]
template = getTemplate("C:\\proj\\proteinModelling\\motifSearching\\4aiy_formatted.data", list_res, "|")
protein = getProtein("C:\\proj\\proteinModelling\\motifSearching\\4aiy_formatted.data", "|")
importantCenters = getImportantCenterAtomsNamesAndDist(protein, template)
