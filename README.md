# NuclearFRAP
R package to run model fits and profile likelihood analysis on nuclear FRAP experiments

This package is part of a method to obtain single cell estimations of nuclear transport parameters (including import rate, export rate and nuclear fixed fraction of molecules), called scFRAPs. A description of the method and an application case can be found in: 
L. Durrieu et al., “Characterization of cell-to-cell variation in nuclear transport rates and identification of its sources,” iScience, vol. 26, no. 1, 2023.
The full original data and scripts from this work are available at:
https://data.mendeley.com/datasets/wfhgv3hxn6/1

# Overview of the scFRAPs method

Following, is a brief description of the main steps of the protocol, intended to give a sense of the type of data and use of the package.

* Tag your protein of interest with a fluorescent protein

* Perform a train of 4 sequential FRAPs (Fluorecence Recovery After Photobleaching) in the nucleus of a cell. Each FRAP consist in 5 images of the pre-perturbation state, an intensional photobleaching pulse, and 30 more images to asess tyhe fluorecence recovery. The numbers 5, 4 and 30 are at the moment hardcoded in the package, but of course can be modified. This train can be performed in either independent cells or in mother-daughter pairs.

* Create a directory for the experiment and copy the images obtained from the microscope to the subdirectory /img/

* Process and quantify the images in imageJ, using the Time Series Analizer v2 macro. Save the quantifications as tab delimited text files with the name lev??-time-series.xls. Save a reference image, with the cell and the quantified regions drawn on it, with the name lev??.stk.jpg.

* Create the following files:

- levData (or posData): a tab delimited text file, with the description of each position. This file should have columns 'lev' (or 'pos') with the position number, 'pb' with the photobleaching order, and other columns describing the strain, cytokinesis state, etc.
- OIF-date.txt: text file with the image capture date of each OIF file. This can be created by running:
sfk filter -ls+"ImageCaputre" -file .txt > OIF-date.txt
There is a script for creating this file available at the NuclearFRAP extras project
- lev??-tif-fname.txt: text file with the file names of all the tif files of the position
There is a script for creating this file available at the darksideoftheshmoo/NuclearFRAP extras project

* Run the fitting and analysis of the data using NuclearFRAP
There is an example script for running the analysis, and a sample folder with all the required files at the darksideoftheshmoo/NuclearFRAP extras project

* Analise the fitted parameters
There are some analysis scripts at LuDurrieu/NuclearTransportRates.



