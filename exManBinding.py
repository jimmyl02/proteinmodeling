# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:28:30 2019

@author: jimmyli2002

Example of using template matching on 2CIXBC4 and 2CJ2BC8 MAN binding sites
"""

from proteinModel import pdbToData, getTemplate, getProtein, getImportantCenterAtomsNamesAndDist
from configparser import ConfigParser

"""
HYPERPARAMETERS
"""

config = ConfigParser()
config.read("config.ini")

if len(config.sections()) == 0:
    config.add_section("proteinModel")
    
config.set("proteinModel", "kClosestAtoms", "10")
config.set("proteinModel", "distanceTestRadius", "100")
config.set("proteinModel", "atomicDistanceTolerance", "0.075")
config.set("proteinModel", "atomicAngleTolerance", "0.1")

with open("config.ini", "w") as f:
    config.write(f)

#Generate data file for protein 2cix and 2cj2

pdbToData("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cix.pdb", "C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cix_formatted.data")
pdbToData("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2.pdb", "C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data")

#List important residues in template (2CIXBC4 as template)

list_res = [("SER", "A", 239), ("SER", "A", 242)]

#Generate the template data from 2cix

template = getTemplate("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cix_formatted.data", list_res, "|")

#Generate the protein data from 2cj2

protein = getProtein("C:\\proj\\proteinModelling\\motifSearching\\exManBinding\\2cj2_formatted.data", "|")

#Perform matching operation

importantCenters = getImportantCenterAtomsNamesAndDist(protein, template)
