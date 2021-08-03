# Ultrasound_Liver_Regeneration_Simulations

This repository contains all code and data needed to reproduce the results of the paper:

"Longitudinal ultrasound imaging and computational modeling reveal sex-dependent suppression of liver regeneration after resection in alcoholic liver disease"

The repository contains two folders, "Auxilliary_Code_and_data" and "Main_code," which should be downloaded into a common repository.
Check the headings of all "Main_code" files for lists of packages that need to be installed prior to use.

For example: 
>library(package_name) 

indicates that

>install.packages("package_name") 

should be run, if that package has not already been installed to your local verison of R. 



Also, each code contains a line to set the working directory 

>setwd("~/")

which must be edited to indicate the path of the downloaded files from this repository. For the Rmarkdowns, this is substituted for the corresponding line

>opts_chunk$set(root.dir ="~/")


To reproduce the results of Barnhart et al., run the following R scripts in order:

Collect_Parameter_Scan_Volume_Data_windows.R
Weighted_Parameter_Scans_windows.R
Unweighted_Parameter_Scans_windows.R *optional*

Copies of these are included for linux users which allow better functionality for parallel processing. The first R script will run in a few hours. The weighted and unweighted parameter scans can take up to a week to complete on an 8-core windows PC. 

Rmarkdown files are included to plot the results of the simulations, which should be run within RStudio. 
