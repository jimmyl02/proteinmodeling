# Protein motif searching program

This is a program which uses a known motif on a certain protein and attempts to generalize and identify similar or exact matches of the motif in another protein. The data should be taken from the wwPDB and is based on spatial understanding and 3D template matching

## Current methods

 - Distance and Name - Attempts to find templates using a center atom and finding a section of the protein which has atoms displaced from a selected center atom a similar distance to the template.

![Graphic explaining distance matching](https://i.ibb.co/Lv3Bqsj/distance-And-Name-Matching-Graphical.png)
## Methods to explore

 - 3D template matching wtih eigenspaces
 - Utilizing angle between atoms to reduce search space

## Usage

NOTE:The **proteinModel<span></span>.py** file should not be modified unless attempting to change how the program works

The main file is **proteinModel<span></span>.py** which contains all methods and code needed to load a template from a protein. The functions are documented within the file. An example of using the MAN binding site on protein with PDBID 2CIX can be found in the **exManBinding<span></span>.py** file and shows how the config options are set and how to load and use data. In the example, the template from 2CIX is attempted to match and find all sites on 2CJ2. 

The following is an explanation of the hyperparameters to be set in the config file

 - kClosestAtoms - Defines how many atoms to check, in example 2, kClosestAtoms = 5
 - distanceTestRadius - Defines the search radius as illustrated in the diagram
 - atomicDistanceTolerance - How different the distance from a possible center atom can be compared to the template information
 - atomicAngleTolerance (To be implemented) - How different the angle from a possible center atom can be compared to the template information

## Examples

### 2CIX and 2CJ2 for MAN binding sites

The code for the visualization can be found in **dataVisualizer<span></span>.py**

The following is the data extracted from a comparison between proteins with pdbID 2CIX and 2CJ2 with the template extracted from 2CIX.

Here is the protein with each individual atom graphed

![Protein atoms plotted](https://i.ibb.co/j4Ghpts/protein2-CJ2.png)

Here are the important center atoms deemed with an atomic distance tolerance of 0.1 graphed with matplot3D

![Plotted important centers of 2CJ2](https://i.ibb.co/0tkXTpD/important-Centers2-CJ2.png)

Here are the known point clouds of functional sites for MAN binding plotted in colors along with the predicted important center atoms
NOTE: There are more MAN binding sites, included are only four examples

![Point clouds and predicted sites plotted](https://i.ibb.co/nkhm5qq/important-Centers-And-Point-Clouds.png)
