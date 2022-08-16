Folder CeDAR_reproduction contains all R scripts used in simulation and real data analysis in CeDAR paper. 

There are three folders: src, data, analysis. 
1. Folder src store all functions used to summarize result or generate figures.
   Specifically, if the newest cedar cannot be accessed through TOAST package (due to a slow update process), please source the file functions_cedar_2.0.R    in folder '/src' or download      the latest version through following command line: devtools::install_github("ziyili20/TOAST", build_vignettes=TRUE).
   
2. Folder data store datasets used in simulation and partial data used in realdata analysis. For data sets not available in this folder, please check the instructions in corresponding script to get access to the data. (You may also look for the data through the following Dropbox link: https://www.dropbox.com/sh/nvbcyr8sebgh9a3/AACnupEOGyNheQwHm7W_j9rra?dl=0) 

3. Folder analysis contains all R scripts that perform simulation and real data analysis. 
   Before running provided code:
   a. Please make sure working path has been correctly set (a.k.a setting the variable main.path at beginning of code to the folder CeDAR_reproduction).
   b. Please make sure all required packages have been installed and loaded

Some tips:
1. For the simulation, reported results are average of 50 simulations. However, it could take long to run. Here we suggest to set the simulation number 
   smaller, if computation time is a concern.
2. For real data anlaysis, like RA anlaysis and smoke status analysis, the data size is too large to run a local computer. We suggest to run it on a cluster.
3. If multiple R scripts in same folder, please run the code in following order: preprocessing -> run_csDM -> result_summary (these are keywords in script name)
4. Results maybe slightly different even with same random seed under differnt systems, but the trend and conclusion should be same.
5. All analyses were performed on R 4.1.3 (Mac), 4.0.2 (Linux). For R 4.2.0, some packages (like CellMix) installation may confront error. We are still working on this to figure it out.

Map between result in Manuscript and folder
/analysis/simulation/computation_time: Table S9
/analysis/simulation/data_noise: Table S5, S6
/analysis/simulation/de_pattern: Figure 4; Figure S2, S3, S4; Table S2, S3
/analysis/simulation/est_prop: Figure S8; Table S8
/analysis/simulation/est_tree_prob: Figure S5; Table S4
/analysis/simulation/mis_tree: Figure S6, S7; Table S7
/analysis/simulation/sample_size: Figure 3; Table S1
/analysis/real-data/blood_gender: Figure 6; Figure S9
/analysis/real-data/blood_ra: Figure 7
/analysis/real-data/blood_sle: Figure S11
/analysis/real-data/blood_smoke: Figure S12
/analysis/real-data/brain_ds: Figure S10
/analysis/real-data/brain_gender: Figure 5
/analysis/real-data/pval_corr: Figure 1; Figure S1


Please email luxiao.chen@emory.edu, if there is any question about the CeDAR reproduction source.


