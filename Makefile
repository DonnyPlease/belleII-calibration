	INC=include/InvariantMassMuMuIntegrator.h  include/InvariantMassMuMuStandAlone.h  include/Splitter.h  include/tools.h include/ChebFitter.h include/nodes.h
SRC=src/InvariantMassMuMuIntegrator.cc	src/InvariantMassMuMuStandAlone.cc src/main.cc src/ChebFitter.cc src/nodes.cc
OBJ=obj/InvariantMassMuMuIntegrator.o	obj/InvariantMassMuMuStandAlone.o obj/main.o obj/ChebFitter.o obj/nodes.o

# location of the eigen library for linear algebra
EIGEN=/cvmfs/belle.cern.ch/el7/externals/v01-10-02/include

# link object files together
exec: $(OBJ) 
	g++ -g -O3  $^   `root-config  --glibs` -fopenmp  -o exec 

# compile source files to object files
obj/%.o: src/%.cc $(INC)
	g++ -g -c -O3 -Iinclude -I$(EIGEN) $<   `root-config --cflags `  -fopenmp   -o  $@ 
