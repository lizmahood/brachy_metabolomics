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
All of my samples were named as such: Hydro.Ctrl.Leaf_5, Sym.SporeW.Root_2, etc. The nomenclature is: group.condition.organ_replicate. This nomenclature is essential to run certain scripts properly (anything comparing a condition to its control). Control conditions must be labelled as Ctrl (see above). If you have different genotypes, put them in "group".

For most of these scripts, if you are wondering what arguments they take, just launch the script without any arguments at all. It will tell you what arguments to provide, in what order.

I should also note that these scripts are not comprehensive of everything I've done with the Brachypodium data. Some things (file renaming/moving, removing the "dead volume", etc.) are very hard-coded for my data and thus are not included in this repo. Check **metabolomics_processing** (the original repo) if you need these.

### HOW TO RUN THESE SCRIPTS:

##### 1_making_TICS.R
Inputs are: folder of mzml files, file with name conversions and desired output folder. This script can be run before MSDIAL.

Please find example file with name conversions in example_files.

Then, you must do MSDIAL processing. The parameters I used for MSDIAL are found in example_files. These may not work for everyone. Remember to assign your samples to groups!! From MSDIAL, you should output an MGF of the aligned results and the PeakHeight, PeakArea, and SN_matrix files. THEN: REMOVE THE FIRST 4 ROWS OF THESE 3 FILES.

##### 2_normalize_metabolites.R
Inputs are: too many to list here. There are 20 inputs. Run the file without any inputs/arguments to see them. 
_Explaination of selected inputs:_
2) Index of "Alignment.ID" column
3) If you don't have any treatment columns in your data, put "0"
17) If you only had one extraction, and thus have one set of blanks to use, you will only have one blank averaged column (which is the average of all metabolite values in your blank samples). If this is your case, just put the index of that column.
19) If you are interested in a certain metabolite, put its Alignment ID in here -- this will tell you if your metabolite was removed at a certain step. If you have multiple, separate with a comma.

This outputs a file with noisy/blank metabolites removed, and metabolite values normalized however you wanted to. It also outputs replicate correlation files. The script removes lowly correlated replicates -- but if you find replicates that passed my correlation filter, yet are still outliers -- those should be removed.

##### graphing/pcas.R
I recommend doing this next so that you can find any additional outliers based on PCAs. These should be removed before you find differentially abundant metabolites (done next).
Inputs: 1) Path to filtered/normalized peak file 2) Output directory 3) All experiments (put yes unless your normalized peak file only contains samples of one group) 4) Do you want your tissues/organs to be graphed together. If "no" -- right now it's hard-coded to look for leaf and root samples, you may have to change that. Additionally, there are only so many different groups you can have. If you have more than 10 conditions, you will have to put in more colors. 5) Put "metab"

##### 3_change_vsn_negative
If you used VSN in step 2 -- it changes 0 values to negative values. This script changes them back to 0 (but only if they were 0 before VSN). To check this, you run step 2 again, but change argument 13 to "no", but keep everything else the same. You can skip this if you did NOT do VSN.
Inputs: 1) VSN normalized file 2) file with everything else done except VSN

Outputs another file with negative values changed to 0.01

##### 4_find_normalized_fold_change_cutoffs.R
If you did not use VSN, you can skip this one also. 
Our method for finding differentially abundant peaks is complicated by VSN. The cutoff for fold change is based on 
normalized and non-normalized values. 
For normalized cutoffs: the script gets all peak area values with fold change between -2 and 2. Then, it gets the fold changes of these peak areas after normalization. It finds the 99th and 1st percentile. These are the top and bottom cuttoffs, respectively, for normalized fold change. Non-normalized fold change had to be at least 2 or -2. Normalized fold change must be above top cutoff or below bottom cutoff. Both fold changes must be in the same direction.

Inputs: 1) File of normalized peaks 2) File of non-normalized peaks 3) Put "yes"
This doesn't output anything, but prints the normalized fold change cutoffs

##### 5_get_diff_abund_metabolites.R
Inputs: 1) Normalized (with VSN) peak file 2) Input file with KNN done (and all filtering) but not VSN 3) output name for volcano plots 4) Put "yes" 5) top cutoff (output by previous script) 6) Bottom cutoff 

This outputs: 1) Differentially abundant metabolites for each condition (if there are any) 2) p-values of all metabolites changing in that condition 3) IDs of any/all metabolites that are differentially abundant 4) volcano plots (one for each condition).

Although I've never done it without VSN, if you're in that boat, I'd suggest: Put the same input file in twice, put 2 for top cutoff, and -2 for bottom cutoff. _This hasn't been tested._

##### Launch Canopus
Go through the scripts in the launching_canopus folder. First the mgf is filtered so that metabolites filtered out by script 2 (and any adducts) are removed. Then launch the 2nd script (which actually does Canopus). The third script gets the output and associates "confidence scores" with superclass/class/subclass/level.5/most.specific_class predictions. The fourth script removes predictions if their confidence score is <0.5 Finally it makes a couple of graphs.

After all this, the most important output file you'll have is canopus_summary.tsv_filtered.tsv

##### 6_split_DAMs_into_updown_merge_canopus.R
As this script requires the canopus_summary.tsv_filtered.tsv, Canopus needs to be completed before it. 
Inputs: 1) Folder containing DAM files (made by script 5) 2) canopus summary filtered file 3) desired output file name
Outputs: One file with DAMs across all conditions, whether or not they are up or down accumulated, and their Canopus class.

##### 7_new_information_theory.py
This script calculates global information theory metrics, and metrics per canopus class. 
Inputs: 1) Input normalized peak area alignment file, with negative values turned to 0 2) output name 3) metab OR trans (put metab) 4) Do you want RDPI broken down by superclasses? none OR path to canopus_summary file

Run it with the canopus_summary.tsv_filtered.tsv file FIRST, then change the name of "RDPI_per_stress_replicate.tsv" to "RDPI_per_stress_replicate_class_averaged.tsv". Then run without the canopus ID file. This will output another "RDPI_per_stress_replicate.tsv" file (global RDPI)

##### graphing/plotting_info_theory.R
The inputs are like this:
1) XX_shannon_di_per_sample.tsv 2) XX_RDPI_per_stress_replicate.tsv 3) XX_RDPI_per_stress_replicate_class_averaged.tsv 4) XX_pk_area_change_per_metabolite.tsv 5) XX 6) XX_RDPI_per_class_and_stress_replicate.tsv

##### 8_get_metabolites_condition_specific.py
Inputs: 1) input file (alignment file) to find condition specific metabolites in 2) Desired output filename (and directory) 3) What omics type? trans OR metab (put metab)

##### MS/MS Networking & identification
This section is qualitative. It assumes you have done MS/MS Networking through MSFINDER, and imported the resulting node and edge files into Cytoscape for network creating. Once this is done, select the nodes in the main network and output a node file of exclusively these. This network_node file is used as an input for some of these scripts. Other scripts assume GNPS FBMN has been performed on your data -- the identification scripts require the .graphml file that GNPS creates. Open this file in Cytoscape, and export the node table. Again, this node table is an input to some of these scripts. 

###### molecular_networking/1_update_nodetable_with_canopus.R
Inputs: 1) mgf file used in molecular networking 2) the network node table output by MSFINDER 3) canopus_summary.tsv_filtered.tsv (this must have come from the mgf input here as argument 1)

Outputs: A tsv node table that can be imported into cytoscape to visualize how superclasses cluster together.

###### identification/1_identification_gnps_msdial_canopus.R
This script assumes you have done identification with MSDIAL (using some msp files when you processed the data) and GNPS.
Inputs: 1) Fully processed and filtered peak area file, with identification 2) Canopus file 3) GNPS' .graphml cytoscape node file 4) folder of DAMs 5) file of our MSFINDER network nodes 6) output folder and file name prefix

Outputs: A file with "high-confidence" identifications (with identification scores >0.8) from both MSDIAL and GNPS. The conditions that each identified metabolite is differentially abundant in (if any) are noted.

###### molecular_networking/2_update_nodetable_with_ids.R
Inputs: 1) Node file updated with canopus superclass (output of script 1_update_nodetable_with_canopus.R) 2) File of identified network compounds (output of 1_identification_gnps_msdial_canopus.R) 3) Name and prefix of output file

Outputs: a tsv node table that can be imported into cytoscape to visualize how superclasses and identifications cluster together.

###### identification/2making_final_identification_files_for_paper.R
This is a purely optional script. It just reformats the identifications file.
Inputs: ARGS: 1) Identification GNPS/MSDIAL and canopus file 2) canopus file (all canopus annotations) 3) output file prefix

Outputs: reformatted identifications file.

###### identification/3_check_overlap_canopus_cf.R
This script is also purely optional, and it requires you to have the ChemOnt classifications of your identifications. I obtained these by: selecting the INCHIs from GNPS' graphml file (you must first import this into Cytoscape and export the node table), getting the INCHIKEY from this using [chembl_ikey](https://github.com/mnowotka/chembl_ikey), and getting the ChemOnt classifications through [ClassyFire Batch](https://cfb.fiehnlab.ucdavis.edu/). 

Once you get your annotations, paste them into your identifications file, and rename to columns to be: CF_superclass, CF_class, CF_subclass, and CF_level.5. Then you can run this script.

Input: The identification file with ClassyFire annotations pasted into it
Outputs: Nothing, but prints % overlap between ClassyFire annotations and CANOPUS annotations. SO -- if 100 compounds were annotated at the class level by both CANOPUS and ClassyFire, and 86 of those got the same class annotation, the overlap is 86%.

##### WGCNA
###### 1WGCNA_clustering_script.R
Most of the code in here came from the [WGCNA Tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/). The traits file (lines 26-32) is hard-coded for my data, you may have to change those lines.
Inputs: 1) Peak Area file 2) Path and prefix of output

###### 2_make_wgcna_module_plots.R
This script is hard-coded to my data again. Lines for new users to change are 141-157.
Many inputs, please run the script with no inputs to see the necessary inputs.

###### 3_calculate_wgcna_mod_enrichment.py
This script determines if any metabolite class is enriched in any of your WGCNA modules.
Inputs: ARGS: 1) WGCNA module output file (output of 1WGCNA_clustering_script) 2) canopus_summary.tsv_filtered.tsv 3) peak area file 4) output path and prefix

### This takes you through all major scripts
Some scripts haven't been covered here -- these are probably going to be used less-often. The 2 scripts in the main directory that were left out were for making a cosine score matrix and for making the input data for our BAR eFP browser. There are other scripts in the graphing folder that are also not on this readme. Let me know if you have any questions about these: ehm79 AT cornell DOT edu

### Thank You!