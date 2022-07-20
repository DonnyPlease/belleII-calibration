#include <iostream>

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include <fstream>
#include <tuple>

#include <Eigen/Dense>
#include <chrono>


using namespace std;
using namespace Eigen;

void save_matrix_to_csv(const MatrixXd &matrix, string name)
{
    IOFormat Precision(FullPrecision, DontAlignCols, ",", "", "", ",\n");
    ofstream file("matrices/" + name);
    file << matrix.format(Precision) << endl;
    file.close();
}

MatrixXd readMatrixFromCsv(string name, int numberOfCols, int numberOfRows) 
{
    string x;
    MatrixXd matrix(numberOfRows,numberOfCols);
    
    ifstream dataObj("matrices/" + name,ios::in); //opening the file.
    if (dataObj.is_open()) { //if the file is open
        for (int i = 0; i < numberOfRows; ++i) {
            for (int j = 0; j < numberOfCols; ++j) {
                getline(dataObj, x, ',');
                matrix(i,j) = stod(x);
            }
	    getline(dataObj, x, '\n');
        }
    }
    else cout << "Unable to open file"; //if the file is not open output
    return matrix;
}


void mergeMatrices(int numOfCores, int size) 
{
    const int numberOfParts = numOfCores/10;  //number of parts in one 
    const int sizeOfOnePart = size / numberOfParts; //nuber of rows to calculate in train and test matrix in one epoch
    
    int j = 0;


    MatrixXd trainMatrix(0,size), testMatrix(0,size);

    for (int i = 0; i < numOfCores; ++i){
        MatrixXd trainMatrixToConcat = readMatrixFromCsv("train_matrix" + to_string(i) + ".csv", size, sizeOfOnePart);
        MatrixXd testMatrixToConcat = readMatrixFromCsv("test_matrix" + to_string(i) + ".csv", size, sizeOfOnePart);
        
        MatrixXd mergeTrain(trainMatrix.rows() + trainMatrixToConcat.rows(), trainMatrix.cols());
        mergeTrain << trainMatrix, trainMatrixToConcat;
        MatrixXd mergeTest(testMatrix.rows() + testMatrixToConcat.rows(), testMatrix.cols());
        mergeTest << testMatrix, testMatrixToConcat;


        trainMatrix = mergeTrain;
        testMatrix = mergeTest;

        if ((i+1) == (sizeOfOnePart + j*sizeOfOnePart)) {
            save_matrix_to_csv(trainMatrix, "train_matrix" + to_string(j) + ".csv");
            save_matrix_to_csv(testMatrix, "test_matrix" + to_string(j) + ".csv");
            ++j;
            trainMatrix.resize(0,size);
            testMatrix.resize(0,size);
        }
    }
}

int main(int argc, char *argv[]) {
    mergeMatrices(stoi(argv[1]), stoi(argv[2]));
    return 0;
}
