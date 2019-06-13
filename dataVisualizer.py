# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 09:22:46 2019

@author: jimmyli2002

THIS IS AN EXAMPLE SPECIFIC TO MAN BINDING SITES ON 2CIX and 2CJ2, THE FUNCTIONS AND IMPLEMENTATION CAN BE APPLIED TO OTHER PROTEINS
"""

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

plt.figure()
ax = plt.axes(projection="3d")

"""
SECTION
Name: Helper functions
"""

"""
FUNCTION
Name: get3DCoordsProtein
Description: Get all 3D coordinates for one protein
"""

def get3DCoordsProtein(file_input, split_char):
    
    returnX = []
    returnY = []
    returnZ = []
    
    with open(file_input, "r") as f:
        
        for line in f:
            
            splitLine = line.split(split_char)
            
            returnX.append(splitLine[4])
            returnY.append(splitLine[5])
            returnZ.append(splitLine[6])
            
    return (returnX, returnY, returnZ)

"""
FUNCTION
Name: getTemplateToPlot
Description: Get all 3D coordinates needed to plot a template
Input: file path, list of residue tuples [(resName, chainId, resSeq)], splitting character
"""

def getTemplateToPlot(file_name, list_res, split_char):
    
    xTemplate = []
    yTemplate = []
    zTemplate = []
    
    with open(file_name, "r") as f:
        
        for line in f:
            
            splitLine = line.split(split_char)
            
            resName = splitLine[1]
            chainId = splitLine[2]
            resSeq = int(splitLine[3])
            
            for residue in list_res:
                
                if resName == residue[0] and chainId == residue[1] and resSeq == residue[2]:
                    
                        #The line we are currently reading matches the requirements for any one of the residues
                        
                        xTemplate.append(float(splitLine[4]))
                        yTemplate.append(float(splitLine[5]))
                        zTemplate.append(float(splitLine[6]))
                        
    #xTemplate, yTemplate, zTemplate filled correctly
        
    xTemplateToPlot = [float(i) for i in xTemplate]
    yTemplateToPlot = [float(i) for i in yTemplate]
    zTemplateToPlot = [float(i) for i in zTemplate]
    
    return xTemplateToPlot, yTemplateToPlot, zTemplateToPlot

"""
SECTION
Name: Usage
"""

"""
EXAMPLE
Graphing one protein
"""

xprotein, yprotein, zprotein = get3DCoordsProtein("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", "|")

xProteinToPlot = [float(i) for i in xprotein]
yProteinToPlot = [float(i) for i in yprotein]
zProteinToPlot = [float(i) for i in zprotein]

#ax.scatter(xProteinToPlot, yProteinToPlot, zProteinToPlot, c="r", marker="^")

"""
EXAMPLE
Graphing various templates
"""

#First binding site example with color yellow

list_res = [("THR", "A", 238), ("SER", "A", 239), ("GLU", "A", 266)]
xTemplateToPlot, yTemplateToPlot, zTemplateToPlot = getTemplateToPlot("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", list_res, "|")

ax.scatter(xTemplateToPlot, yTemplateToPlot, zTemplateToPlot, c="y", marker="o")

#Second binding site example with color magenta

list_res = [("ALA", "A", 151), ("ASP", "A", 152), ("GLU", "A", 155), ("THR", "A", 250)]
xTemplateToPlot, yTemplateToPlot, zTemplateToPlot = getTemplateToPlot("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", list_res, "|")

ax.scatter(xTemplateToPlot, yTemplateToPlot, zTemplateToPlot, c="m", marker="o")

#Third binding site example with color chocolate

list_res = [("GLU", "A", 155), ("GLN", "A", 159), ("SER", "A", 251), ("PRO", "A", 253)]
xTemplateToPlot, yTemplateToPlot, zTemplateToPlot = getTemplateToPlot("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", list_res, "|")

ax.scatter(xTemplateToPlot, yTemplateToPlot, zTemplateToPlot, c="chocolate", marker="o")

#Fourth binding site example with color lightcoral
#NOTE: Was not found with atomicDistanceTolerance 0.1

list_res = [("THR", "A", 252), ("ILE", "A", 261)]
xTemplateToPlot, yTemplateToPlot, zTemplateToPlot = getTemplateToPlot("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", list_res, "|")

ax.scatter(xTemplateToPlot, yTemplateToPlot, zTemplateToPlot, c="lightcoral", marker="o")


"""
EXAMPLE
Graphing individual atoms
"""
#ax.scatter(-23.721, 50.792, 41.744, c="g", marker="^")
#ax.scatter(-23.245, 57.127, 30.323, c="y", marker="^")

"""
EXAMPLE
Given a list of centers in the format returned by getImportantCenterAtomsNamesAndDist, graph them
"""

list_important_centers = [(('O', '3.384', '10.402', '-8.056'), 0.594), (('O', '8.293', '8.931', '-7.639'), 0.768), (('O', '9.077', '9.156', '-10.686'), 0.722), (('O', '5.799', '40.470', '4.358'), 0.737), (('O', '1.888', '26.997', '-2.653'), 0.674), (('O', '-1.906', '21.678', '10.735'), 0.602), (('O', '-5.291', '27.351', '-8.109'), 0.647), (('O', '-0.867', '28.576', '-9.963'), 0.595)]

for site in list_important_centers:
    ax.scatter(float(site[0][1]), float(site[0][2]), float(site[0][3]), c="g", marker="^")
    
"""
EXAMPLE
Setting the limits on the axis of the 3D graph
"""

ax.set_xlim3d(-20, 40)
ax.set_ylim3d(0, 50)
ax.set_zlim3d(-30, 30)

plt.show()
