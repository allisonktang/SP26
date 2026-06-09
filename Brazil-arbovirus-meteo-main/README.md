# Brazil-arbovirus-meteo

This repository contains code to recreate the spatiotemporal regression modelling study in the paper 'Spatiotemporal relationships between extreme weather events and arbovirus transmission across Brazil'. 

## Code
All analyses can be run from the main.R file, which calls scripts with functions to run the analysis (contained in the R folder). 

## Data
You need to download the data separately from Zenodo (doi 10.5281/zenodo.13286671) and save in the data folder. Please cite both the paper and the data if including them in your work. 

Key for how the municipality was assigned to a case in the data files:


 * _resid files. The municipality assigned to a case was the municipality of residence for all data entries.
 * _report files. The municipality assigned to a case was the municipality where a case was reported or hospitalised.
 * _exp_resid files. The municipality assigned to a case was the municipality of suspected infection for 64%, 58%, and 66% of the data entries in the CHIKV, DENV and ZIKV datasets respectively which reported this information, and for the remaining data entries we used the municipality of residence. This is the data used in the main analysis.
 * _exp_report files. The municipality assigned to a case was the municipality of suspected infection for 64%, 58%, and 66% of the data entries in the CHIKV, DENV and ZIKV datasets respectively which reported this information, and for the remaining data entries we used the municipality where a case was reported or hospitalised.

