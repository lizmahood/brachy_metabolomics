# Welcome to brachy_metabolomics
***
## A set of scripts for processing and visualizing metabolomics data, used in Mahood et al, 2022.
_While these scripts were tested for my data, they have not been tested against any other data. Please let me know if you encounter errors._

### NOTES ON USAGE:
These scripts assume that you are working on a Windows computer and have processed your metabolomic data with [MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html "MSDIAL Home")

These scripts launch the programs CANOPUS/Zodiac (both part of the [SIRIUS4 package](https://bio.informatik.uni-jena.de/software/sirius/ "SIRIUS Home")). 
They assume you have downloaded and set up SIRIUS4.

Finally, they launch [MS-FINDER's](http://prime.psc.riken.jp/compms/msfinder/main.html "MSFINDER Home") molecular networing functionality. You must have MS-FINDER installed for this.

**All of the above softwares are free!**

From MSDIAL, you should output an MGF of the aligned results and a PeakHeight or PeakArea file. THEN: REMOVE THE FIRST 4 ROWS OF THIS FILE.