# Protein motif searching program

This is a program which uses a known motif on a certain protein and attempts to generalize and identify similar or exact matches of the motif in another protein. The data should be taken from the wwPDB and is based on spatial understanding and 3D template matching

## Current methods

 - Distance and Name - Attempts to find templates using a center atom and finding a section of the protein which has atoms displaced from a selected center atom a similar distance to the template.

![Graphic explaining distance matching](https://i.ibb.co/Lv3Bqsj/distance-And-Name-Matching-Graphical.png)
## Methods to explore

 - 3d template matching wtih eigenspaces
 - Utilizing angle between atoms to reduce search space

## Usage

The main file is **proteinModel<span></span>.py** which contains all methods and code needed to load a template from a protein. The functions are documented within the file. An example of using the MAN binding site on protein with PDBID 2CIX can be found in the **exManBinding<span></span>.py** file and shows how the config options are set and how to load and use data. In the example, the template from 2CIX is attempted to match and find all sites on 2CJ2. 

The following is an explanation of the hyperparameters to be set in the config file

 - kClosestAtoms - Defines how many atoms to check, in example 2, kClosestAtoms = 5
 - distanceTestRadius - Defines the search radius as illustrated in the diagram
 - atomicDistanceTolerance - How different the distance from a possible center atom can be compared to the template information
 - atomicAngleTolerance (To be implemented) - How different the angle from a possible center atom can be compared to the template information

