# NuclearFRAP
R package to run model fits and profile likelihood analysis on nuclear FRAP experiments.

This package is part of a method to obtain single cell estimations of nuclear transport parameters (including import rate, export rate and nuclear fixed fraction of molecules), called scFRAPs. A description of the method and an application case can be found in reference (1).

The full original data and scripts from this work are available at:
https://data.mendeley.com/datasets/wfhgv3hxn6/1

A step-by-step protocol on the experiments and fittings can be found at reference (2).

# Overview of the scFRAPs method

Following is a brief description of the main steps of the protocol, intended to give a sense of the type of data and use of the package.

* Tag your protein of interest with a fluorescent protein

* Perform a train of 4 sequential FRAPs (Fluorecence Recovery After Photobleaching) in the nucleus of a cell. Each FRAP consist in 5 images of the pre-perturbation state, an intensional photobleaching pulse, and 30 more images to asess tyhe fluorecence recovery. The numbers 5, 4 and 30 are at the moment hardcoded in the package, but can be modified in the script. This train can be performed in either independent cells or in mother-daughter pairs.

* Create a directory for the experiment and copy the images obtained from the microscope to the subdirectory /img/

* Process and quantify the images in imageJ, using the Time Series Analizer v2 macro. Save the quantifications as tab delimited text files with the name lev??-time-series.xls. Save a reference image, with the cell and the quantified regions drawn on it, with the name lev??.stk.jpg.

* Create the following files:

** levData (or posData): a tab delimited text file, with the description of each position. This file should have columns 'lev' (or 'pos') with the position number, 'pb' with the photobleaching order, and other columns describing the strain, cytokinesis state, etc.
** OIF-date.txt: text file with the image capture date of each OIF file. This can be created by running:
sfk filter -ls+"ImageCaputre" -file .txt > OIF-date.txt
There is a script for creating this file available at the NuclearFRAP extras project
** lev??-tif-fname.txt: text file with the file names of all the tif files of the position
There is a script for creating this file available at the darksideoftheshmoo/NuclearFRAP extras project

* Run the fitting and analysis of the data using NuclearFRAP

* Analyse the fitted parameters
There are some analysis scripts available at LuDurrieu/NuclearTransportRates.

# Installation

1.	Download and install R (https://cran.r-project.org). 
2.	Although not strictly necessary, it is recommended to download and install RStudio (https://posit.co).
3.	Install the following package dependencies (available at the CRAN repository): R.utils, R.oo, R.methodsS3, R.plyr and R.doBy
4.	Install the package Rcell2 (https://github.com/darksideoftheshmoo/rcell2) (3).
5.	Download the NuclearFRAP R package (here!).
6.	Install NuclearFRAP from a local folder.

# Usage

There is an example script for running the analysis and a sample folder with all the required files at the darksideoftheshmoo/NuclearFRAP extras project

## Load the NuclearFRAP package

If you already have the package installed, load it with

> library(NuclearFRAP)

Load the GeomVFrac database (this is a single database containing geometrical information for the cells in all the experiments)

> geomVfrac.db <-read.delim("/path/to/analysis/geomVfrac-db.dat")

## Load a nuclear FRAP experiment to R

After  finishing all the manual analysis of an experiment, you should end with a directory containing the
following  files (see the paper (2) for more details):

a.	levData (or posData): a tab delimited text file, with the description of each position. 
b.	OIF-date.txt: text file with the image capture date of each OIF file. 
c.	lev??-time-series.xls: tab-delimited text file created by ImageJ's plugin Time Series Analyzer v2
d.	lev??-tif-fname.txt: text file with the file names of all the .tif files of the position
e.	lev??.stk.jpg: reference image, with the cell and the quantified regions drawn on it.

To load all this data into R, use the loadFrapPairs function of NuclearFRAP. This function returns a
cell.data object that can be manipulated with Rcell's functions.

> d<-loadFrapPairs("path/to/my/data")

If you have the original tif images of the experiment in a di erent folder than the  les described above,
you can specify an img.path parameter to loadFrapPairs.

> d<-loadFrapPairs("path/to/my/data", img.path="path/to/my/images")

In some experiments, because of the way the were analyzed, there is one cell per position. This means
that the mother might be in pos=1 and the daughter in pos=2. In order to deal with all datasets in a
transparent manner, you can restructure these datasets. In the restructured data mother and daughter will
have the same position number. In order to do the restructuring, the 'levData.txt'  le must have a 'pairs'
columns, with the pair number.

> d<-restructure.pairs(d)

If you haven't got a dataset to load, you can still go through this vignette using the example dataset.

> data(YFP)

Select the specific geometric data for the experiment. Note that the strain and date information must be identical to the one used in the GeomVfrac database. For example:

> d$experiment <- "Ace2_S288C[2023-05-01]"
> geom.db<-subset(geomVfrac.db,experiment=="Ace2_S288C/2023-05-01")

## Raw data visualization

Once the data has been loaded (and restructured) you can look at it using standard methods, such as ggplot and Rcell2's functions.

## Fitting the model

Optional: Define the parameters of the experiment and the model. Here can be defined the type (fitted or fixed) of the parameters, and their values and search boundaries.

> param.FRAP.pairs <- list(VERSION.STR="51"
                         ,MODEL.PARAMETERS=c("log10_kEVfrac","log10_kI","tot0","Fn0","Fc0","PBn","PBfrac","log2_geomVfrac","auto0n","auto0frac","an1","an2","ac0","ac1","ac2","log2_geomVfrac_mean","log2_geomVfrac_sd")
                         ,FIXED.PARAMETERS=c(log10_kEVfrac=0,log10_kI=0,tot0=0,Fn0=0,Fc0=0,PBn=0,PBfrac=0,log2_geomVfrac=0,auto0n=0,auto0frac=1,an1=0,an2=0,ac0=1,ac1=1,ac2=1,log2_geomVfrac_mean=1,log2_geomVfrac_sd=1)
                         ,VALUE.PARAMETERS=c(log10_kEVfrac=NA,log10_kI=NA,tot0=NA,Fn0=0,Fc0=0,PBn=NA,PBfrac=0,log2_geomVfrac=4,auto0n=0,auto0frac=1,an1=0,an2=0,ac0=1,ac1=0,ac2=0,log2_geomVfrac_mean=3.1,log2_geomVfrac_sd=0.434)
                         ,LOWER.PARAMETERS=c(log10_kEVfrac=NA,log10_kI=NA,tot0=NA,Fn0=0,Fc0=0,PBn=0,PBfrac=0,log2_geomVfrac=0,auto0n=0,auto0frac=0.1,an1=-0.1,an2=-0.001,ac0=0.5,ac1=-0.1,ac2=-0.001,log2_geomVfrac_mean=-3,log2_geomVfrac_sd=NA)
                         ,UPPER.PARAMETERS=c(log10_kEVfrac=NA,log10_kI=NA,tot0=NA,Fn0=1e4,Fc0=1e4,PBn=1,PBfrac=1,log2_geomVfrac=8,auto0n=1e4,auto0frac=10,an1=0.1,an2=0.001,ac0=1.5,ac1=0.1,ac2=0.001,log2_geomVfrac_mean=6,log2_geomVfrac_sd=NA)
                         ,PRE.ACTIVATION.FRAMES=5 
                         ,TOTAL.FRAMES.PER.FRAP=35 
                         ,AUTO.F=0 
                         ,F.COR.CYT=1 
                         ,SEC.PER.FRAME=0.22 #time between consecutive frames
                         ,PB.CORRECTION.FUN=function(fl,time) fl 
                         ,SIGMA.FUN=function(fl) fl
                         ,SIGMA.METHOD="SINGLE.CELL.PAIR"
                         #,SIGMA.METHOD="SIGMA.FUN"
                         #,SIGMA.METHOD="FROM.DATA"
                         ,FRAP.ORDER.SEP="-"
                         ,FRAP.ORDER.DAUGHTER.STR="daughter"
                         ,FRAP.ORDER.MOTHER.STR="mother"
                         ,MOTHER.FRAP.NUMBER=4
                         ,DAUGHTER.FRAP.NUMBER=4
                         ,TOT0.SIGMA.FACTOR=3
                         ,GEOMVFRAC.SIGMA.FACTOR=3
                         ,GEOMVFRAC.REL.SD=0.2
                         
                         #parameters used to load the experiment
                         ,ref.img="lev*.stk.jpg"
                         ,time.series="lev*-time-*.xls"
                         ,fname="lev*-tif-fname.txt"
                         ,desc="levData.txt"
                         ,OIF.date="OIF-date.txt"
                         ,channel.name="YFP"
)

> d$parameters <- updateParameters(d$parameters,geom.db)


To fit the model to the data, you can use the function runModelFit.
As most of NuclearFRAP's functions, there is also a 'pair' version of this function, called runModelFitPair.

> param <- runModelFitPair(d[[pos==1,]],p=d$parameters)

The first argument for runModelFitPair is a data.frame, with the data for the pair. This can be easily done with d[[pos==1,]], that returns a data.frame for pos 1. The second parameter is a list containing a whole bunch of parameters for the  t. Check the help page for ?FRAP.options. This list (with default values) is saved in the parameters slot of the cell.data object. The function returns a data.frame with columns type (indicating if it is the mother or daughter cell),  t.index (the index of the  t, for each cell), the cost of the fit and the values of all the parameters of the model.
Once we have the set of optimal parameters, we'll probably want to see how well these parameters (and model structure)  t the data. For this we can use the function getSimDataPair.

> sim.db<-getSimDataPair(d[[pos==1,]],param,p=d$parameters)
> sim.db<-restructureSimDataPair(sim.db)

The getSimDataPair function returns a data.frame with the nuclear and cytoplasmic fluorescence variables for the mother and daughter cell (and the expected standard deviation), for each time point. Because of the way the returned data.frame is organized, there are a lot of NAs. You can easily restructure this data.frame with the function restructureSimDataPair, to a format more convenient for plotting.

## Profile Likelihood Analysis

To get an idea of the uncertainty in the parameter estimation, we can do a pro le likelihood analysis (4). The profile likelihood can be done in one dimension,  xing the value of a single parameter and leaving the others free to minimize the cost function. More interestingly for us, it can be done (and easily visualized) in two dimensions,  xing the value of two parameters, and leaving the others free.
We are particularly interested in the kinetic parameters kI and kEVfrac (see model NuclearFrapModelParameterization.pdf). The region of uncertainty will be the region where the profile likelihood cost is equal to or lower than the minimal cost, plus a threshold. If we use a chi-square cost function, and the residuals are normally distributed, we can get this threshold from the chi-square distribution with an alpha = 0.05 and
the degrees of freedom equal to the total number of free parameters of the model (4).
NuclearFRAP has two functions to calculate the pro le likelihood uncertainty region. The  first one isrunProfileLikelihoodContour (and runProfileLikelihoodContourPair). This function works only for 2D profile likelihood analysis. Starting from the best  t point, it draws 'rays' in different directions (theta) and calculates where the increased profile likelihood cost equals the threshold (boundaryCost). It returns a data.frame with one point per direction. These points lie in the boundary of the uncertainty region. This method works well if the region is sort of O shaped, but fails if it is U shaped.
The second method is implemented in runProfileLikelihoodGrid. This method evaluates a profile like-likelihood cost in a regular grid in the parameter space. It works for both, 1D and 2D profile likelihood analysis. If the regular grid is not passed as an argument, runProfileLikelihoodGrid calls runProfileLikelihoodContour to get an estimation of the region it should sample. From the regular grid, the 2D con dence region can
be calculated. This method can deal with U shaped regions but is more computationally intense. The following command can take several minutes to complete.

> PL.db<-runProfileLikelihoodGridPair(d[[pos==2,]],param,p=d$parameters,scan=c("kEVfrac","kI")
+ ,boundaryCost=qchisq(p=c(0.05), df=9))
  
You can plot the 2D profile likelihood region with the plotPL2D function.

> plotPL2D(PL.db,x="kEVfrac",y="kI",param=param,signifCost=c(qchisq(p=c(0.05), df=9)))

For each cell, two contours are plotted. The solid line corresponds to the contour calculated from the regular grid of points, calculated by runProfileLikelihoodGrid. The dashed line is the contour as calculated by runProfileLikelihoodContour. The points indicate the best  t for the mother 'm1' or daughter 'd1'. If several local minimum are found they are also shown in this plot as 'm2', 'd3', etc. The chi-square cost of each point is indicated.

## Analyze FRAP pair

Some 'high level' functions are included in the package. One of these is analyzeFrapPair. This function calculates the  t and pro le likelihood analysis on kI and kEVfrac, for a pair of cells.

> fp<-analyzeFrapPair(d[[pos==1,]],d$parameters)

analyzeFrapPair returns a frapPair object, that contains all the data from the analysis. The package   also provides a plot method for this class (Figure 5).

> plot(fp)

![image](https://github.com/darksideoftheshmoo/NuclearFRAP/assets/79978782/199d1f6e-07b0-45e3-b2ad-5cf19fe2a0b7)

The created figure has several panels. The top left panels have the traces (points and thin lines) for the nuclear and cytoplasmic fluorescence, for mother and daughter. The panel title specifies on which cell the FRAP was done. These panels also have the model predictions, as solid lines of the color of the correspondent variable. When several local minima are  nd, the different predictions are plotted on these same panels.
The gray-shaded region corresponds to the model prediction of the standard deviation of the data. See the section on the determination of sigma for more details. The dashed black lines are the swapped predictions.

This means that the model is fitted to the daughter's data, using the optimal values of kI and kEVfrac for the mother (and leaving all other parameters free). If this line looks like a good  t to the data, this strongly suggests that the daughter's kinetic parameters are NOT significantly different to the mother's. The same analysis is done for the mother, using the fixing the kinetic parameters to the daughter's optimal value.
The following panel to the right is the 2D pro le likelihood analysis, for kI and kEVfrac. The local minima for the mother and daughter and the profile likelihood regions are represented as explained before (see the Pro le Likelihood Analysis section). The orange curves represent the instrumental boundaries for the kinetic parameters. The dashed line in the upper right is the points that satisfy tau = Tsampling (where tau=1/(kI + kEVfrac), and Tsampling is the time gap between consecutive data points). The dashed curve in the lower left corner is that satisfy tau = Texperiemnt (where Texperiement is the length of the FRAP experiment). The dotted curves are analogous, but use Tsampling/2 and 2*Texperiemnt as boundaries.
The table in the higher right panel gives the values of all parameters for the different local minima found.
The histograms in the lower left panels are the distribution of the residuals of the model (region shaded in §gray), and the normal distribution using the mean and standard deviation of the residuals (green curve). The violet curve is the normal distribution, using mean zero and the standard deviation predicted by the model (see the section on the determination of sigma for more details). If the prediction of the sigma prediction is
working, then the violet and green curve should coincide. These panels can also be used to determine if the residuals are normally distributed.
In the lower right corner the time course for the residuals are shown, for each variable. The black dotted lines are the predicted standard deviation for the data. These panels can be used to see if there is any tendency in the residuals, indicative that the fit is not good.

## Estimating Sigma

The package has two methods to estimate sigma (the expected standard deviation of the data). The first one called "SIGMA.FUN", uses a user de ned function to calculate sigma from the fluorescent level of the data. To use this method, modify the parameters slot of the cell.data object.

> d$parameters$SIGMA.METHOD <- "SIGMA.FUN"
> d$parameters$SIGMA.FUN <- function(fl) 25 + 32 * fl

An alternative method, and the one used by default, is to estimate the noise from the "unused" data.
When doing the FRAP on the daughter, the mother variables are being acquired, but they are not used for the model  fitting. So the idea is to use them for the estimation of sigma. Because we want the residuals, we need some "prediction" to calculate them. The first 5 points of each FRAP are assumed to be constant and the remaining points are fitted to a loess curve. From this coarse model, the residual distribution is calculated. The sigma level is assumed to be constant for each cell.
For the second cell to be photobleached, this method works well. For the first cell to be photobleached, it tends to underestimate the standard deviation, probably because after the photobleaching there is a lower signal and consequently a lower noise level. To correct for this, the package assumes the noise is proportional to the fluorescence level and corrects the sigma by the ratio of the fluorescence levels. This correction is only applied to the first cell to be photobleached. To use this method set SIGMA.METHOD to "SINGLE.CELL.PAIR".

> d$parameters$SIGMA.METHOD <- "SINGLE.CELL.PAIR"


## Profile Likelihood Matrix

A way to check if all parameters are identifiable is to do a pro le likelihood matrix, i.e. a matrix of plots where in each panel a 2D pro le likelihood contour is shown, for the respective variables. In the diagonal,1D pro le likelihood cost curves are shown. To calculate all these profile likelihood analyses, you can use the high-level function profileLikelihoodFrapPair and the correspondent plot method (Figure 6). Beware this is a very computationally expensive thing to do. It takes 12hs to complete a single experiment (with 10 pairs) in my Intel i5 quad-core.

> PL <- profileLikelihoodFrapPair(d[[pos==1,]], d$parameters,
+ vars=c("kEVfrac","kI","tot0","PBn","PBc","auto1n"), pairwise = TRUE)
> plot(PL)

# Contributions

Alan Bush (@abush84) developed the model, wrote the initial version of the NuclearFRAP package, and wrote the documentation. Andrea Katz updated the package to run fully in R and with current packages and R version. Lucía Durrieu (@ludurrieu) developed the model, tested the package and wrote the documentation.

# References

1.	Durrieu, L., Bush, A., Grande, A., Johansson, R., Janzén, D., Katz, A., Cedersund, G., and Colman-Lerner, A. (2023). Characterization of cell-to-cell variation in nuclear transport rates and identification of its sources. iScience 26. 10.1016/j.isci.2022.105906.
2.	Durrieu, L., Bush, A.,and Colman-Lerner, A. (2023). FRAP-based estimation of nuclear import and export rates in single yeast cells.
3.	Méndez, N.A., Beldorati, G., Constantinou, A., and Colman‐Lerner, A. (2023). rcell2: Microscopy‐Based Cytometry in R. Curr. Protoc. 3. 10.1002/cpz1.726.
4.	Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmüller, U., and Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. Bioinformatics 25, 1923–1929. 10.1093/bioinformatics/btp358.


