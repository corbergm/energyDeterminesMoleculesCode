**readme for the code associated to the manuscript 
"How energy determines spatial localisation and copy number of molecules in neurons"
by Cornelius Bergmann, Kanaan Mousaei, Silvio Rizzoli and Tatjana Tchumatchenko.**

# 1. System requirements 

The only required program is MATLAB 2022b with the "Parallel Computing Toolbox" installed. It can also be used without the latter with decreased performance by replacing every instance of "parfor" by "for".

Software has successfully been tested on various machines with Windows as operating system.

No specific hardware is necessary to run the code.

# 2. Installation guide 

Ensure that MATLAB2022b and the aforementioned toolboxes are installed. 
First, change the directory to the desired location of the code, then open a terminal and clone the repository or download it.

First, change the directory to the repository folder and run "initPath_MATLAB". 
It adds the whole project folder with subfolders to the MATLAB path. 

# 3. Demo 

An exemplary simulation of mRNA and protein distributions and the associated costs can be found in "master_exampleSimulation.m". 
The expected run time is less than a minute on an average PC (depending on the chosen parameters).

# 4. Instructions for use

Data downloads
------------------------------------------------------------

For Figure 3, download and store the following datasets in the "data" folder after renaming them as follows:
> "Table S1" of Tushev et al., 2018 (https://doi.org/10.1016/j.neuron.2018.03.030) as "Tushev_2018.xls"
> "Supplementary Data 2" from Zappulo et al., 2017 (https://doi.org/10.1038/s41467-017-00690-6) as "Zappulo_2017.xslx"
> "Supplementary Table 2" from Helm et al., 2021 (https://doi.org/10.1038/s41593-021-00874-w) as "Helm_2021.xslx"
> "Supplementary Data 1" from Fornasiero et al., 2018 (https://doi.org/10.1038/s41467-018-06519-0) as "Fornasiero_2018.xlsx"

For Figure 4, download and store the following datasets in the "data" folder after renaming them as follows:
> "Supplementary Data 2" from Zappulo et al., 2017 (https://doi.org/10.1038/s41467-017-00690-6) as "Zappulo_2017.xslx"
> "Supplementary file 3" from Perez et al., 2021 (https://doi.org/10.7554/elife.63092) as "Perez_2021.xslx"
> Raw data for mRNA from http://linnarssonlab.org/cortex/ associated with Zeisel et al., 2015 (https://doi.org/10.1126/science.aaa1934) as "Zeisel_2015.txt"

For Figures 1B, 1C
------------------------------------------------------------

To create a panel, execute the corresponding script in the folder "figure_1". 

For Figure 3
------------------------------------------------------------

Before creating a figure, make sure that all required datasets have been put in the "data" folder with the correct names (see section "Data downloads" above).

To create the figure, execute the corresponding script in the folder "figure_3". 

For Figures 1E, 1F, 2, 4
--------------------------------------------------------------------

To run the simulations on the full parameter space as described in the manuscript, execute "master_list" for a single dendrite length at a time. Parallelization is established within the function "run_list". 
Do not execute the script for two dendrite lengths in one run, because the output formats of the "costPerSeg" variable (encoding the cost per 10 micron of dendrite) depend on the respective dendrite length and their sizes are hence pairwise incompatible.

The run duration for one dendrite length and the parameters used throughout the manuscript is some minutes on an average PC when parallelized using the "Parallel Computing Toolbox".

For rapid access, the simulation output for 250, 500, 750, and 1000 micron dendrites with the parameters used throughout the manuscript is precomputed and can be found in the "files" folder under the date "2024_02_03".

Before creating a figure, make sure that all required datasets have been put in the "data" folder with the correct names (see section "Data downloads" above).

Finally, to create a figure, execute the corresponding script in the folder "figure_X".

For Figure 5
------------

To simulate the temporal evolution of somatically labelled Shank3 protein in dendrites and spines, run the script "master_Shank3Dynamics". The actual implementation of the explicit and implicit finite difference scheme is in the function "run_Shank3Dynamics".

On a single computing node on an HPC, the run duration is some hours for both the explicit and the implicit numerical schemes.

For rapid access, the simulation output with the parameters used in the manuscript is precomputed and can be found in the "files" folder under the date "2023_05_04".

Finally, to create the figure, execute the corresponding script in the folder "figure_5". 

# Project structure   

"initPath_MATLAB" adds the whole project folder with subfolders to the MATLAB path. 

"/data" contains all experimental data we retrieved from online databases and publications (data sources see below)

"/files" contains precomputed simulation outputs to recreate the manuscript figures.

"/functions" contains all scripts and functions used to perform simulations, analyze data or create figures.
  
> "/figure_X" contains all scripts that create, when run, at least one of the panels of figure X.
  
> "/model_core" contains all the core functions (i.e., those actually implementing or solving equations). 

>> "get_distribution" solves the main mRNA and protein equations

>> "get_cost" computes the associated metabolic cost,

>> "get_ensembleDiffConstFit", "get_steadyStateSol3stateModel", and "get_analyticalSteadyStateSol1StateModel" are needed to fit an ensemble diffusion constant to the 3-state transport model as explained in the manuscript.
  
>> "master_multipleDiffConstFit", "master_list", and "run_list" are needed to run the simulations for a whole list of parameter combinations.
  
# Data sources

Let us clarify that this code contains data from previously published manuscripts that were manually derived from these:

"/data/cumulativeMitochondriaDensitySTD_Lopez-Domenech_2016"
> https://doi.org/10.1016/j.celrep.2016.09.004

"/data/mRNADensity100Microns_Fonkeu_2019"
> https://doi.org/10.1016/j.neuron.2019.06.022

"/data/proteinDensity100Microns_Fonkeu_2019"
> https://doi.org/10.1016/j.neuron.2019.06.022

"/data/Shank3Dynamics_binned_Tsuriel_2006"
> https://doi.org/10.1371/journal.pbio.0040271
