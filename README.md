# belleII-calibration
Repository for a bachelor thesis
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
This repository contains code used in my bachelor thesis.
It also contains data on which the code was demonstrated and some python scripts 
that were used to create some of the graphs in the thesis.
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
1. First we need to cut out the data we want to work with - muon_peak.cpp
2. Divide the data into time intervals with split.
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
3. Because of a big time demand of the the unbinned invariant mass fitting, there 
was a need for paralelizing the process of computing the matrices of the loss.
Running the script run.sh launches the whole process - it needs two parameters as input
-first, one needs to specify number of cores used - should be a multiple of 11.
-second, number of divisions - determined by the number of intervals from dividing the time.
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
4.Run the algorithm by calling ./alg
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
5.Run pyhton scripts for graphs.
