# belleII-calibration

Repository for a bachelor thesis:

This repository contains code used in my bachelor thesis.
It also contains data on which the code was demonstrated and some python scripts 
that were used to create some of the graphs in the thesis.

#------------------------------------------------------------------------#

1. First we need to cut out the data we want to work with - muon_peak.cpp


3. Divide the data into time intervals with split.


5. Because of a big time demand of the the unbinned invariant mass fitting, there 


was a need for paralelizing the process of computing the matrices of the loss.
Running the script run.sh launches the whole process - it needs two parameters as input
-first, one needs to specify number of cores used - should be a multiple of 11.
-second, number of divisions - determined by the number of intervals from dividing the time.

4.Compile the algorithm using g++ aloritghm.cpp -o out -fopenmp 
(it will not work without -fopenmp, it is needed for the paralelization of the process of fitting - in case of n>1000 the algortihm starts to take some time)

5.Run algorithm with "./alg". Inputs have to be stored in the folder results(because they are results of previous calculations. The form and names of the matrices can be seen there. Results also get stored in the folder "results".

5.Run pyhton scripts for graphs.
