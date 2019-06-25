# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:10:43 2019

@author: jimmyli2002
"""

from os import path
import glob
import re
import requests
import pandas as pd
import numpy as np

"""
FUNCTION
Name: downlaodRequiredPDB
Description: Downloads all pdb files if they do not exist already
"""

def downloadRequiredPDB(input_file):

    with open(input_file, "r") as input_data:
        
        matchPattern = re.compile("PDBID (\S+)")
        
        for line in input_data:
            
            regexMatch = matchPattern.findall(line)
            
            if regexMatch:
                
                pdbID = regexMatch[0]            
            
                if not path.exists(r"./pdbdata/" + str(pdbID) + ".pdb"):
                    
                    url = "https://files.rcsb.org/download/" + pdbID.upper() + ".pdb"
                    req = requests.get(url)
                    print("Downloading: " + pdbID.upper() + ".pdb")
                    
                    with open(r"./pdbdata/" + pdbID.upper() + ".pdb", "wb") as output_file:
                        
                        output_file.write(req.content)
                
# Each pdb file will get its own training file to make sure we do not lose data during sliding window
                        
"""
FUNCTION
Name: generateBaseAndMetaData
Description: Create the base csv file with all information and moves the most negative to left coordinate. Also creates metadata file
"""

def generateBaseAndMetaData(list_pdbid_files):
    
    for pdbID in list_pdbid_files:
        
        if not path.exists(r"./pdbdata/" + str(pdbID) + ".pdb"):
            
            print("WARNING: The file " + str(pdbID) + ".pdb does not exist, base data not being generated")
            
            continue
        
        completeMol = []
        resName = []
        resNum = []
        atomX = []
        atomY = []
        atomZ = []
        occupancy = []
        tempFac = []
        label = []
        
        # Variables to store the metadata
        xShift = 100000.0
        yShift = 100000.0
        zShift = 100000.0
        
        
        with open(r"./pdbdata/" + str(pdbID) + ".pdb", "r") as f:
            
            match_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s+(\S{1,4})\s+([A-Z]{1,4})\s([\s\w])\s+(\d+)\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d{0,2})?\s*(\-?[0-9]\d{0,2}\.\d*)?")
            
            for line in f.read().splitlines():
                
                if line.startswith("ENDMDL"):
                    
                    break
                
                match_list=match_pattern.findall(line)
                
                if match_list:

                    completeMol.append(match_list[0][0])
                    resName.append(match_list[0][1])
                    resNum.append(int(match_list[0][3]))
                    atomX.append(float(match_list[0][4]))
                    atomY.append(float(match_list[0][5]))
                    atomZ.append(float(match_list[0][6]))
                    occupancy.append(match_list[0][7])
                    tempFac.append(float(match_list[0][8]))
                    label.append("Other")
                    
                    xShift = min(xShift, float(match_list[0][4]))
                    yShift = min(yShift, float(match_list[0][5]))
                    zShift = min(zShift, float(match_list[0][6]))
                    
        xShift *= -1
        yShift *= -1
        zShift *= -1
        
        # Save meta data
        
        with open(r"./training_data/" + str(pdbID) + "_meta.txt", "w") as f:
            
            f.write(str(xShift) + "|" + str(yShift) + "|" + str(zShift))
        
        for i in range(len(atomX)):
            
            atomX[i] += xShift
            atomY[i] += yShift
            atomZ[i] += zShift
                        
        raw_data = {"mol_name": completeMol,
                    "res_name": resName,
                    "res_num": resNum,
                    "atom_x": atomX,
                    "atom_y": atomY,
                    "atom_z": atomZ,
                    "occupancy": occupancy,
                    "temp_fac": tempFac,
                    "label": label
                    }
        
        # Save the data to a dataframe and then to CSV
        df = pd.DataFrame(raw_data, columns = ["mol_name", "res_name", "res_num", "atom_x", "atom_y", "atom_z", "occupancy", "temp_fac", "label"])
        df.to_csv(r"./training_data/" + str(pdbID) + "_data.csv", index=False)
        del df

"""
FUNCTION
Name: generateDataPerPDB
Description: Fill out the known information into the csv, replacing old values
"""
                        
def generateDataPerPDB(data_file, label_name):
    
    match_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s+(\S{1,4})\s+([A-Z]{1,4})\s([\s\w])\s+(\d+)\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d*)?\s*(\-?[0-9]\d{0,2}\.\d{0,2})?\s*(\-?[0-9]\d{0,2}\.\d*)?", re.MULTILINE)
    
    with open(data_file, "r") as data_file_content:
        
        sections = data_file_content.read().split("\n\n")
        
        currIndex = 0
        
        for section in sections:
            
            currIndex += 1
            
            pdbIdRegex = re.compile("^PDBID (.*)$", re.MULTILINE)
            pdbIdSearch = pdbIdRegex.findall(section)
            
            if pdbIdSearch:
                
                pdbID = pdbIdSearch[0]
                
                if not path.exists(r"./training_data/" + str(pdbID) + "_data.csv"):
                    
                    print("WARNING: Need to download " + str(pdbID))
                    
                    continue
            
                # We now have the PDBID of the section and must extract the atom info to exchange. Must use metadata
                
                data = pd.read_csv("./training_data/" + str(pdbID) + "_data.csv")
                
                match_list = match_pattern.findall(section)
                
                with open("./training_data/" + str(pdbID) + "_meta.txt") as metadata_file:
                    
                    metadataInfo = metadata_file.read().split("|")
                    
                    xShift = float(metadataInfo[0])
                    yShift = float(metadataInfo[1])
                    zShift = float(metadataInfo[2])
                    
                    for match in match_list:
                        
                        shiftedX = float(match[4]) + xShift
                        shiftedY = float(match[5]) + yShift
                        shiftedZ = float(match[6]) + zShift
                        
                        data.loc[np.isclose(data["atom_x"], shiftedX) & np.isclose(data["atom_y"], shiftedY) & np.isclose(data["atom_z"], shiftedZ), ["label"]] = label_name
                    
                data.to_csv("./training_data/" + str(pdbID) + "_data.csv", index=False)
                
                if currIndex % 10 == 0:
                    print("Complete with [" + str(currIndex) + "/" + str(len(sections)) + "] sections")
                        
#downloadRequiredPDB("./manBindingSites.txt")

# Steps to get all downloaded pdb files            
list_pdbId = []
pdbFileMatch = re.compile(r"\\(.*).pdb")
for file in glob.glob("./pdbdata/*"):
    list_pdbId.append(pdbFileMatch.findall(file)[0])
    
# Use generated list to call generateMetaData
generateBaseAndMetaData(list_pdbId)

generateDataPerPDB("./manBindingSites.txt", "MAN")
