# Mean species responses predict effects of environmental change on coexistence
---

These data include: 
1/code (.R files and .py files) to reproduce the results in https://doi.org/10.1111/ele.14278. 
2/macroinvertebrate data (csv file) and accompanying metadata (xlsx file) that were used to produce Fig.5.

## Description of the data and file structure

Root.R: Please open this file first and run it in R to load all packages and the self-made functions to do all simulations and analyses. 
Bugs_tools.R and Theory_tools.R: the self-made functions to do all simulations and analyses, ran by running Root.R. 
MVN.py: Python script to work with multivariate normal distributions. Ran from witihin R, by running Root.R. To do so, you need to have the R reticulate package installed.

Once Root.R has been run, please go to:
Theory_computations.R: Line by line, you can execute the script in R to reproduce the theoretical results. This file is commented throughout.
Bug_computations.R: Line by line, you can execute the script in R to reproduce Fig.5. This file is commented throughout.

Bugs_Data.csv: Contains the macroinvertebrate data used to produce Fig.5. It contains the following variables:
lat: Latitude (decimal degrees)
long: Longitude (decimal degrees)
period: Five year time span in which location was sampled
stress; Classification of sampling location environmental stress, based on land use
chem_ph; pH (standard pH units)
chem_cond; Conductivity (uS/cm)
chem_turb; Turbidity (NTU)
chem_doc: Dissolved Organic Carbon (DOC) (mg/L)
chem_calcium: Calcium (mg/L)
chem_ammonia_n: Amount of N as Ammonia (mg/L)
chem_nitrate_n	: Amount of N as Nitrate (mg/L)
chem_silica: Silica (mg/L)
Ablabesmyia - Thienemannimyia: Genera with 0/1 absence/presence indicating community composition for given location
community: 0/1 absence/presence of members within the given community of the entire dataset, without taxonomic identifiers
community_size; Number of unique genera found within given community

Community_Metadata: Contains the metadata just described above, plus the variable types (e.g. continuous, character)

Some of the cells in Bugs_Data.csv are NA. This means that the variable has not been measured at that specific location and time (i.e. Not Available).

## Sharing/Access information

The contents of this Dryad submission is also available on https://github.com/fdelaend/De_Laender_et_al_feasibility. 

## Code/Software

Everything is done from R, although . For a description of the workflow, please read Description of the data and file structure.