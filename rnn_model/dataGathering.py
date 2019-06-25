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
Name: generateBaseData
Description: Create the base csv file with all information and moves the most negative to left coordinate
"""

def generateBaseData(list_pdbid_files):
    
    for pdbID in list_pdbid_files:
        
        if not path.exists(r"./pdbdata/" + str(pdbID) + ".pdb"):
            
            print("WARNING: The file " + str(pdbID) + ".pdb does not exist, meta data not being generated")
            
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
        
        
        with open(r"./pdbdata/" + str(pdbID) + ".pdb", "r") as f:
            
            match_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s+(\S{1,4})\s+([A-Z]{1,4})\s([\s\w])\s+(\d+)\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?\s+(\-?[0-9]\d{0,2}\.\d*)?")
            
            for line in f.read().splitlines():
                
                if line.startswith("ENDMDL"):
                    
                    break
                
                match_list=match_pattern.findall(line)
                
                if match_list:
                    
                    completeMol.append(match_list[0][0])
                    resName.append(match_list[0][1])
                    resNum.append(match_list[0][3])
                    atomX.append(match_list[0][4])
                    atomY.append(match_list[0][5])
                    atomZ.append(match_list[0][6])
                    occupancy.append(match_list[0][7])
                    tempFac.append(match_list[0][8])
                    label.append("Other")
                        
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
Name: generateMetaData
Description: Generate meta-data for each pdb file inputted and store it
"""

def generateMetaData(list_pdbid_files):
    
    for pdbId in list_pdbid_files:
        
        #TEMP, NEED TO DECIDE WHAT META DATA IS NEEDED
        #Currently have
        # pdbId, xShift, yShift, zShift
        
        print("TODO")

"""
FUNCTION
Name: generateDataPerPDB
Description: Fill out the known information into the csv, replacing old values
"""
                        
def generateDataPerPDB(data_file, list_pdb_name):
    
    with open(data_file, "r") as data_file_content:
        
        sections = data_file_content.read().split("\n\n")
        
        for section in sections:
            
            pdbIdRegex = re.compile("^PDBID (.*)$", re.MULTILINE)
            pdbIdSearch = pdbIdRegex.findall(section)
            
            if pdbIdSearch:
                
                pdbId = pdbIdSearch[0]
                
                print(pdbId)
                
                if not path.exists(r"./training_data/" + str(pdbId) + "_data.csv"):
                    
                    print("WARNING: Need to download " + str(pdbId))
                    
                    continue
            
            
#downloadRequiredPDB("./manBindingSites.txt")

# Steps to get all downloaded pdb files            
list_pdbId = []
pdbFileMatch = re.compile(r"\\(.*).pdb")
for file in glob.glob("./pdbdata/*"):
    list_pdbId.append(pdbFileMatch.findall(file)[0])
    
# Use generated list to call generateMetaData
#generateBaseData(list_pdbId)
    
#generateMetaData(list_pdbId)
            
generateDataPerPDB("./manBindingSites.txt", list_pdbId)
