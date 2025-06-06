# Evaluating-native-cytomegaloviruses-as-vectors
This repository contains source code and original data in support of the paper "Evaluating native cytomegaloviruses as vectors for a spillover-disrupting Lassa virus transmissible vaccine"

## The repository contains five folders
1. CPP_Source. This folder contains the C++ source code used to perform Bayesian MCMC. The code is set up to analyze either individual data sets or a numerically indxed sequence of data sets.
2. Chains. This folder contains the thinned chains produced by applying the Bayesian MCMC algorithm. The folder contains two files, one with the chains for MnatCMV2 and another with the chains for MnatCMV3. Because each virus was analyzed independently, each has its own set of chains.
3. Mathematica. This folder contains a Mathematica notebook which derives the likelihood function and other ancillary results.
4. Real_Data. This folder contains the original data, a data key describing the original data column headers, and a stripped down data file for each virus. The stripped down data files for each virus include only age and infection status and are in the format used by the MCMC algorithm.
5. Simulated data. This folder contains the 100 simulated data sets used to test our method. 
