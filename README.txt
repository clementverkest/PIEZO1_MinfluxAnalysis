# PIEZO1 MINFLUX analysis
Matlab scripts and datasets for "Cluster nanoarchitecture and structural diversity of PIEZO1 in intact cells"

by Stefan Lechner, Clement Verkest, Lucas Roettger and Nadja Zeitzschel 
Contact: s.lechner@uke.de

This set of Matlab scripts are used to visualize and explore PIEZO1 DNA-PAINT Minflux data in Verkest et al., 2024 (https://www.biorxiv.org/content/10.1101/2024.11.26.625366v1.full)
The scripts and data are provided for academic and visualization purposes only. Commercial usage of the provided Minflux data (reproduction outside of the publication, etc. ) is forbidden.
Please contact the authors for further inquiries.


# Requirements

Scripts generated and tested with Matlab R2023b and R2024a. Require additional Matlab features such as the Signal_Processing_Toolbox or Image_Processing_Toolbox

Required Main scripts:
-PIEZO1_GFP_ClusterAnalysis_1
-PIEZO1_GFP_ClusterAnalysis_2
-PIEZO1_GFP_ClusterAnalysis_3
-PIEZO1_ALFA_analysis_1
-PIEZO1_ALFA_TrimerInPlaneProjection

and associated subroutines:
-DBSCAN
-dbscan2
-PlotClusterAnalysisResult
-PlotRawData
-SetGFPanalysisParameters_V2
-CalcTrimerAngle1
-CalculateTraceMean
-PIEZO1Superparticle
-PlotTrimerAnalysisResult
-SetALFAanalysisParameters


Raw data provided:

-confocal images as tif files, with region selected for Minflux scan highlighted


-Matlab structure files with raw, unfiltered minflux localizations and traces
	-PIEZO1mGL_CTL_RawData, PIEZO1mGL_OSMO_RawData, PIEZO1_ALFAmGL_SOMA_RawData, PIEZO1_ALFAmGL_NEURITE_RawData and PIEZO1_ALFAmGL_CYTOD_RawData contain data array (X rows, 8 cols) with raw loc and traces, organized as follows:
	col1: X coord /col2: Y coord / col3: Z coord / col4: TID (trace identification number) / col5: efo (Hz) / col6: cfr / col7: background (Hz) / col8: time (s)
	
	
	
-Matlab structure files with selected PIEZO clusters
	-GFP_CTL_all_selected_clusters and GFP_OSMO_all_selected_clusters contain data array (X rows and 4 cols) with the cluster coordinates (X,Y,Z) and TID




# How to use 

-Download folders (Fig1-2 and Fig3-4) and their content
-Add them to your Matlab path

-Run PIEZO1_GFP_ClusterAnalysis_1.m
	-This script is used to visualize Minflux unfiltered and filtered data from DNA-PAINT GFP
	-Select options in the 'CHOOSE OPTIONS' section (set 'true' or 'false')
	-Select the data source (Control or hypo-osmotic) in 'Load data' section, (set 'DataSource = 0', '1', '2', '3' or '4')

-Run PIEZO1_GFP_ClusterAnalysis_2.m
	-This script is used to plot panels in Fig 1f-g, as well as 3D plot of selected PIEZO cluster
	-Select options in the 'set options' section !  Warning ! Plotting the full cluster dataset might slow down things. 
	-Select the data source ('DataSource = 1', '2', '3' or '4')

-Run PIEZO1_GFP_ClusterAnalysis_3.m
	-This script is used to plot panels in Fig 2c-f
	-Select options in the 'set options' section.
 
-Run PIEZO1_ALFA_Analysis_1.m	
	-This script is used to plot most panels in figure 3 and 4, and Extended Data Fig 9.
	-Select options in 'CHOOSE OPTIONS' section
	-Select the data source in the 'load data' section (0 = soma, 1 = neurite, 2 = cytochalasin-D)
	-Possibility to run PIEZO1_ALFA_TrimerInPlaneProjection.m after running PIEZO1_ALFA_Analysis_1.m (it generates 'IndivTrimersRAW', that is required for the script to work)




# Citation
preprint: https://www.biorxiv.org/content/10.1101/2024.11.26.625366v1.full

