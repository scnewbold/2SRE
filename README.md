# 2SRE
This repository contains a collection of R scripts sufficient to apply the two-stage random-effects estimator used in Newbold et al 2024.

vslmeta-sim.R and vslmeta-app.R will replicate all results reported in Newbold et al 2024. Results will be written to files named vslmeta-sim.out and vslmeta-app.out in the same folder as the scripts. On a normally configured desktop PC, vslmeta-sim.R requires ~3 hours to run, and vslmeta-app.R requires ~3 minutes to run.

Most *Fun.R scripts are functions used by other scripts. vslmeta-app.R uses pauseFun.R, tictocFun.R, ghFun.R, seroFun.R, qtestFun.R, trimfillFun.R, ptepeeseFun.R, twosremaFun.R, and twosremrFun.R. vslmeta-sim.R uses tictocFun.R, ghFun.R, and seroFun.R

yhatFun.R is a streamlined script that will apply the 2SRE estimator to a correctly formatted meta-dataset, useful for small scale testing or as a starter script users can adapt for applications to other meta-datasets.

getdataFun.R will simulate a small meta-dataset in the correct format or will download the demonstration meta-dataset used in Newbold et al 2024, useful for small scale testing.

Original VSLs and weights.xlsx is an Excel file that contains the original study VSL estimates and weights for the weighted VSL estimates and standard errors used in the meta-analysis.

Send questions and suggestions for bug fixes to Steve Newbold at snewbold@uwyo.edu
