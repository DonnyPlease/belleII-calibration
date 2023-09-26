# belleII-calibration

Repository for a bachelor thesis:

The code in the folder src is not mine, but it was given to me for the project. Without it the crucial part of the code for this thesis does not work.

This repository contains code used in my bachelor thesis.
It also contains data on which the code was demonstrated and some python scripts 
that were used to create some of the graphs in the thesis.

For most of the programs to work, one need to have some of the ROOT libraries.

#------------------------------------------------------------------------#

1. First we need to cut out the data we want to work with - cutData.cpp

2.Result of the initial division made by a program not included is stored in timeSplits.txt.

3. Matrices are created and saved by running run.sh:
Because of a big time demand of the the unbinned invariant mass fitting, there 
was a need for paralelizing the process of computing the matrices of the loss.
Running the script run.sh launches the whole process - it needs two parameters as input
-first, one needs to specify number of cores used - should be a multiple of 11.
-second, number of divisions - determined by the number of intervals from dividing the time.
Several files - mostly those associated with the process of fitting Cheb. polynomials was an existing code
"borrowed" for this thesis.


#------------------------------------------------------------------------#
Next we can start algorithm.
#------------------------------------------------------------------------#

4.Compile the algorithm using g++ aloritghm.cpp -o out -fopenmp 
(it will not work without -fopenmp, it is needed for the paralelization of the process of fitting - in case of n>1000 the algortihm starts to take some time)

5.Run algorithm with "./alg". Inputs have to be stored in the folder results(because they are results of previous calculations. The form and names of the matrices can be seen there. Results also get stored in the folder "results".

6.Pyhton scripts are used for graphs.
