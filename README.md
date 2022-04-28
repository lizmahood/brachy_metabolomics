# Welcome to brachy_metabolomics
***
## A set of scripts for processing and visualizing metabolomics data, used in Mahood et al, 2022.
_While these scripts were tested for my data, they have not been tested against any other data. Please let me know if you encounter errors._

### NOTES ON USAGE:
These scripts assume that you are working on a Windows computer and have processed your metabolomic data with [MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html "MSDIAL Home").

These scripts launch the programs CANOPUS/Zodiac (both part of the [SIRIUS4 package](https://bio.informatik.uni-jena.de/software/sirius/ "SIRIUS Home")). 
They assume you have downloaded and set up SIRIUS4.

They launch [MS-FINDER's](http://prime.psc.riken.jp/compms/msfinder/main.html "MSFINDER Home") molecular networing functionality. The network is visualized using [Cytoscape](https://cytoscape.org/ "Cytoscape home"), and some Cytoscape output files are needed for full functionality of the molecular networking scripts. To run all molecular networking scripts, you must have MSFINDER and Cytoscape installed.

**All of the above softwares are free!**

### NOTES ON COLUMN NAMES/SAMPLE NAMES:
All of my samples were named as such: Hydro.Ctrl.Leaf_5, Sym.SporeW.Root_2, etc. The nomenclature is: group.condition.organ_replicate. This nomenclature is essential to run certain scripts properly (anything comparing a condition to its control). Control conditions must be labelled as Ctrl (see above).

For most of these scripts, if you are wondering what arguments they take, just launch the script without any arguments at all. It will tell you what arguments to provide, in what order.

I should also note that these scripts are not comprehensive of everything I've done with the Brachypodium data. Some things (file renaming/moving, removing the "dead volume", etc.) are very hard-coded for my data and thus are not included in this repo. Check **metabolomics_processing** (the original repo) if you need these.

### HOW TO RUN THESE SCRIPTS:

##### 1_making_TICS.R
Inputs are: folder of mzml files, file with name conversions and desired output file. This script can be run before MSDIAL.

Then, you must do MSDIAL processing. The parameters I used for MSDIAL are found in ./MSDIAL_parameters. 
From MSDIAL, you should output an MGF of the aligned results and the PeakHeight, PeakArea, and SN_matrix files. THEN: REMOVE THE FIRST 4 ROWS OF THESE 3 FILES.