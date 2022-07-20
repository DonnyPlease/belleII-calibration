#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <thread>

//#include <omp.h>


//may return 0 when not able to detect
//const auto processor_count = 4;//std::thread::hardware_concurrency();

using namespace std;



//We define a class named Data_fit_pars. The variables inside this class are used throughout the entire program.
class Data_fit_pars {
public:
    vector<double>*get_stack_fits() {
        return &stack_fits;
    }

    vector<vector<int>>*get_ind_of_fit() {
        return &ind_of_fit;
    }

    vector<vector<double>> *get_mat_of_errors() {
        return &mat_of_errors;
    }

    vector<double> get_last_row(){
        vector<double>result;
        for (int i = 0; i < mat_of_errors.size(); ++i) {
            result.push_back(mat_of_errors[i].back());
        }
        return result;
    }

    void set_stack_fits(vector<double> new_stack_fits) {
        stack_fits = std::move(new_stack_fits);
    }

    void set_ind_of_fit(vector<vector<int>> new_ind_of_fit) {
        ind_of_fit = std::move(new_ind_of_fit);
    }

    void set_mat_of_errors(vector<vector<double>> new_mat_of_errors) {
        mat_of_errors = std::move(new_mat_of_errors);
    }

    int size() {
        return mat_of_errors.size();
    }

    void clear() {
        stack_fits.clear();
        ind_of_fit.clear();
        mat_of_errors.clear();

    }

    void set_mat_of_errors_index(int i, int j, double d) {
        mat_of_errors[i][j] = d;
    }

    void set_stack_fits_index(int i, double d) {
        stack_fits[i] = d;
    }

    Data_fit_pars() = default;

private:
    vector<double> stack_fits; //stack_fits stands for memory of stacking errors, which is the core of this entire program
    vector<vector<int>> ind_of_fit; //this stores indices corresponding to stacks of errors from above
    vector<vector<double>> mat_of_errors;    //this stores the huge matrix of errors, which is essential for the algorithm to be quick, however, it is quite time consuming to create it

};


/*Function min_of_e(args) gets as input two arrays (a and b) and a class of Data_fit_pars, from which only ind_of_fit will be used. of errors and an array of arrays of breaking points
that belong to those errors. It adds the two array together and finds index with minimum overall error. It takes this
index and appends it to the array of indices which corresponds to minimum error. It returns the error and the new array of indices
length of which is now, by the way, larger by 1.*/
pair<double,int> min_of_e(vector<double>a, vector<double>b, Data_fit_pars *data_fit_pars) {

    vector<double> c(a.size(),0); //Declare
    transform(a.begin(), a.end(), b.begin() , c.begin(), plus<double>()); //sum vectors a and b and store it in c

    int min_ind = (int) distance(begin(c), min_element(begin(c), end(c))); //find index with the smallest element
    double min_e = c[min_ind];  //store the error, which is minimum of the vector

    return make_pair(min_e,min_ind);
}

/*-This is the most important function in the entire algorithm.
  -It works with the number of steps n we are looking for and "global" variables data_fit_pars.
  -It does not return anything, but it edits variables stack_fits and ind_of_fit in class of Data_fit_pars.*/
void error_of_n_minus_1(int n, Data_fit_pars *data_fit_pars) {

    //Declare what is needed.
    vector<double> new_stack_fits; //this is vector where we stack cumulative errors
    vector<vector<int>>new_ind_of_fit = *data_fit_pars->get_ind_of_fit();
    //(data_fit_pars->get_mat_of_errors()->size()-n,vector<int>(n));
    //for(int i = 0; i < data_fit_pars->get_mat_of_errors()->size()-n; i++) {
    //    new_ind_of_fit[i] = {};
    //}
    vector<int> p1;
    double p2;

    //We need to distinct the case of n = 1, because in this case we do not work with the matrix of errors_in_the_middle.
    if (n==1) {
        if (data_fit_pars->get_stack_fits()->operator[](0) == 0) {
            for (int i = 1; i < data_fit_pars->get_mat_of_errors()->size(); i++) {
                data_fit_pars->set_stack_fits_index(i,data_fit_pars->get_mat_of_errors()->operator[](0)[i-1]);   //This is analogous to the
                                                                                // function errors of the end,
                                                                                //but this time it is reversed.
            }
            data_fit_pars->set_stack_fits_index(0,1); //data_fit_pars.stack_fits[0] = 1 , so we are storing the information whether we calculated this step or not.
        }

    } else {
        //Check whether we calculated this case or not.
        if (data_fit_pars->get_stack_fits()->operator[](0) == n) {
            return;

            //If not, start calculating.
        } else {

            error_of_n_minus_1(n - 1, data_fit_pars); // Call the same function but with (n-1) as second parameter.

            new_stack_fits.push_back(data_fit_pars->get_stack_fits()->operator[](0));

            //We do, however, want to edit the rest, therefore we start a for cycle.
            for (int i = n; i < (data_fit_pars->get_mat_of_errors()->size()); i++){

                //Take sub_stack_fits as part of the original stack_fits with size depending on i.
                auto first = data_fit_pars->get_stack_fits()->begin() + 1;
                auto last = data_fit_pars->get_stack_fits()->begin() + i+2-n;
                vector<double> sub_stack_fits(first,last);

                //Take sub_middle as part of mat_of_errors corresponding to the sub_stack_fits.
                vector<double>sub_middle;
                for (int j = n-1; j < i; j++) {
                    sub_middle.push_back((data_fit_pars->get_mat_of_errors()->operator[](j)[i-1]));
                }

                //Here we do something similar to the function min_of_e.
                vector<double> c(sub_stack_fits.size(),0);
                transform(sub_stack_fits.begin(), sub_stack_fits.end(), sub_middle.begin(), c.begin(), plus<double>()); //sum two vectors together and store it in c.
                int min_ind = (int) distance(begin(c), min_element(begin(c), end(c)));
                double min_e = c[min_ind];
                if (n > 2) {
                    vector<int>i_ind = (data_fit_pars->get_ind_of_fit()->operator[](min_ind));
                    i_ind.push_back(min_ind+n-1);
                    new_ind_of_fit[i-n] = i_ind;
                } else {
                    new_ind_of_fit[i-n].push_back(min_ind+1);

                }
                new_stack_fits.push_back(min_e);

            }
            //Replace the original stack_fits and ind_of_fit with the new_stack_fits and new_ind_of_fit.
            data_fit_pars->set_stack_fits(new_stack_fits);
            data_fit_pars->set_ind_of_fit(new_ind_of_fit);
            data_fit_pars->set_stack_fits_index(0,n); //Set stack_fits[0] = n, so we can easily know whether this calculation has been done or not.
        }
    }
}


//
pair<double, vector<int>> fit_steps(int n, Data_fit_pars *data_fit_pars) {
    
    vector<int> indices;
    
    error_of_n_minus_1(n, data_fit_pars);


    auto first = data_fit_pars->get_stack_fits()->begin() + 1;
    auto last = data_fit_pars->get_stack_fits()->begin() + data_fit_pars->get_stack_fits()->size();
    vector<double> sub_stack_fits(first,last);

    vector<double> end_error = data_fit_pars->get_last_row();
    first = end_error.begin() + n;
    last = end_error.begin() + end_error.size();
    vector<double> end_fits(first,last);

    pair<double,int> e_and_ind = min_of_e(sub_stack_fits,end_fits,data_fit_pars);
    vector<int> ind_vec = data_fit_pars->get_ind_of_fit()->operator[](e_and_ind.second);
    ind_vec.push_back(e_and_ind.second+n);

    return make_pair(e_and_ind.first,ind_vec);
}


//This function reads matrices from csv file with comma (,) as delimeter of elements and (,\n) as delimeters of rows.
//Input is the name (or path to) of the csv file and N which is dimension of matrix shaped (N x N).
vector<vector<double>>read_the_matrix_from_csv(string name, const int N) {

    string x; //Here we store temporarily what we read from file.
    vector<vector<double>>matrix; //Declaration of result.
    
    ifstream dataObj(name,ios::in); //opening the file.
    if (dataObj.is_open()) { //if the file is open 
        for (int i = 0; i < N; ++i) {
            vector<double>row; //Declare one row.
            row.reserve(N);   //Reserve size of the row.
            for (int i = 0; i < N; ++i) {
                getline(dataObj, x, ','); //Read from file until encountering a comma and store it in variable x.
                row.push_back(stod(x));   //stod = string to double - converts string to double and pushes it back to the row we are filling in this cycle.
            }
            getline(dataObj, x, '\n'); //new line after
            
            matrix.push_back(row); // Push back the entire row.
        }
    }
    else cout << "Unable to open file"; //if the file is not open output

    return matrix;
}



pair<double,double> train_test_error(pair<vector<vector<double>>,vector<vector<double>>>matrices,
                                     vector<int> breaking_points) {

    double train_error = 0.0;
    double test_error = 0.0;
    
    breaking_points.insert(breaking_points.begin(),0);
    breaking_points.push_back(matrices.first.size());


    for (int i = 0; i < breaking_points.size()-1; ++i) {
        train_error += matrices.first[breaking_points[i]][breaking_points[i+1]-1];
        test_error += matrices.second[breaking_points[i]][breaking_points[i+1]-1];
    }


    return make_pair(train_error,test_error); //to be corrected !!!!!!
}


pair<vector<int>,vector<int>> read_fold_sizes(){
    string x;
    pair<vector<int>,vector<int>>result;


    ifstream dataObj("results/sizes.txt",ios::in); //opening the file.
    if (dataObj.is_open()) //if the file is open
    {

        for (int i = 0; i < 10; ++i) {
            getline(dataObj, x, ',');
            result.first.push_back(stoi(x));
            cout << "Yes" << endl;
            getline(dataObj, x, '\n'); //new line after
            result.second.push_back(stoi(x));
        }
    }
    else cout << "Unable to open file"; //if the file is not open output
    return result;
}

void save_matrix(vector<vector<double>>matrix, string name, int N)
{
    ofstream result(name);
    
    result << setprecision(9);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 10; j++) {
            result << matrix[i][j];
            if (j!=9) result << ","; else result << endl;
        }
    }
    result.close();
    return;
}

//This function looks for the best number of steps the data are supposed to be fitted with.
//Uses K-fold crossvalidation - input K tells it how many folds there are.
int best_number_of_steps(const int N, const int K, const vector<pair<vector<vector<double>>,vector<vector<double>>>>&matrix) {

    vector<vector<double>>train_errors_of_folds;
    vector<vector<double>>test_errors_of_folds;
    vector<double>zeros;

    zeros.reserve(K);
    for (int j = 0; j < K; j++) {
        zeros.push_back(0.0);
    }
    train_errors_of_folds.reserve(N);
    for (int i = 0; i < N; i++) {
        train_errors_of_folds.push_back(zeros);
    }
    test_errors_of_folds.reserve(N);
    for (int i = 0; i < N; i++) {
        test_errors_of_folds.push_back(zeros);
    }



// This is where K-Fold CV starts


    auto  (*globalVariables) = new Data_fit_pars [K];

//    if (K > processor_count) {
//        delete [] globalVariables;
//        globalVariables = new Data_fit_pars [processor_count];
//    }

//#pragma omp parallel for default(none) shared(cout,globalVariables,train_errors_of_folds,test_errors_of_folds,matrix)
    for (int fold_i = 0; fold_i < K; fold_i++) {
        Data_fit_pars data_fit_pars = globalVariables[fold_i];

        //Data_fit_pars data_fit_pars = globalVariables[fold_i % processor_count];



        data_fit_pars.set_mat_of_errors(matrix[fold_i].first);

        vector<double> stack_fits_local(matrix[fold_i].first.size(), 0);
        data_fit_pars.set_stack_fits(stack_fits_local);
        vector<vector<int>> aloc_ind_of_fit(matrix[fold_i].first.size(), vector<int>(N));
        data_fit_pars.set_ind_of_fit(aloc_ind_of_fit);

        for (int i = 0; i < matrix[fold_i].first.size(); i++) {
            data_fit_pars.get_ind_of_fit()->operator[](i) = {};
        }

        for (int number_of_steps = 1; number_of_steps <= N; number_of_steps += 1) {


            //Calling the algorithm calculating best step function
            pair<double, vector<int>> result_of_fit = fit_steps(number_of_steps, &data_fit_pars);

            //Assigning result to vector of training breaking points.
            vector<int> train_breaking_points = result_of_fit.second;


            //Calculate error.
            pair<double, double> train_test_e = train_test_error(matrix[fold_i], train_breaking_points);


#pragma omp critical
            {
                train_errors_of_folds[number_of_steps - 1][fold_i] = train_test_e.first;
                test_errors_of_folds[number_of_steps - 1][fold_i] = train_test_e.second;
                std::cout << setprecision(9) <<  "Fold number:" << fold_i + 1 << "  Number of steps:" << number_of_steps << "      "
                          << train_test_e.first << "........" << train_test_e.second << std::endl;
                vector<int> breaking_points = result_of_fit.second;

                for (int breaking_point : breaking_points) {
                    cout << breaking_point << "  ";
                }
                cout << endl;
            }
            //<< " Thread number: " << omp_get_thread_num()

        }
        cout << endl;
    }
    //Read sizes of folds data as we need this information for comparing the average loss functions.
    pair<vector<int>,vector<int>>sizes = read_fold_sizes();

    //Average of loss
    vector<double> train_errors_steps;
    vector<double> test_errors_steps;

    //save results to file
    ofstream result("results/result.csv");
    
    result << setprecision(9) << "n,train_error,test_error" << endl;
    for (int i = 0; i < N; i++) {
        double tr_e = 0.0; //cumulative train error
        double te_e = 0.0; //cumulative test error
        for (int j = 0; j < K; j++) {
            tr_e += train_errors_of_folds[i][j] / sizes.first[j];
            te_e += test_errors_of_folds[i][j] / sizes.second[j];
        }
        train_errors_steps.push_back(tr_e / K);
        test_errors_steps.push_back(te_e / K);

        if (tr_e != 0)
        {
            cout << "Number of steps: " << (i + 1) << "  Train error: " << tr_e / K << "  Test error: " << te_e / K << endl;
            result << (i+1) << "," << tr_e / K << "," << te_e / K << endl;
        }
    }
    result.close();
    save_matrix(train_errors_of_folds,"results/train-errors-of-folds.csv",N);
    int min_ind = 1 + (int) distance(begin(test_errors_steps), min_element(begin(test_errors_steps), end(test_errors_steps)));
    cout << min_ind << endl;
    return min_ind;
}


//This function saves the breaking points together with the fitted peaks (and its error) of the data on particular intervals.
void save_steps(vector<vector<double>> & means, vector<int> breaking_points){

    //size = number of steps
    int size = breaking_points.size()+1;

    //add zero and number of divisions we calulated breaking points on so we include the zeroth and last interval.
    breaking_points.insert(breaking_points.begin(),0);
    breaking_points.push_back(means.size());
    
    //Read errors from csv file.
    vector<vector<double>>std_errors = read_the_matrix_from_csv("results/stdErrorMatrix.csv", means.size());

    //Open file and store there the results.
    ofstream file("results/"+to_string(size)+"steps.csv");

    file << "means,start,std_errors" << endl; //Header of csv file.

    for (int i = 0; i < breaking_points.size()-1; ++i) {

        double a = means[breaking_points[i]][breaking_points[i+1]-1]; //Extract from matrix the correct element.
        int b = breaking_points[i]; 
        double c = std_errors[breaking_points[i]][breaking_points[i+1]-1];

        file << setprecision(9) << a  << "," << b << "," << c << endl;
    }

    //add the last one the the file as well because the graph needs it. (could be probably done in python script as well)
    file << means[breaking_points[breaking_points.size()-2]][breaking_points.back()-1] << "," 
         << (breaking_points.back()) << "," 
         << std_errors[breaking_points[breaking_points.size()-2]][breaking_points.back()-1] << endl;
}


//After we find what is the optimal number of steps we fit all the data with that number of steps.
void  fitting_ideal_steps(const int N, pair<vector<vector<double>>,vector<vector<double>>> & matrix) {
    
    auto  (*globalVariables) = new Data_fit_pars[1];

    Data_fit_pars data_fit_pars = globalVariables[0];

    data_fit_pars.set_mat_of_errors(matrix.first);

    vector<double> stack_fits_local(matrix.first.size(), 0);
    data_fit_pars.set_stack_fits(stack_fits_local);
    
    vector<vector<int>> aloc_ind_of_fit(matrix.first.size(), vector<int>(N));
    data_fit_pars.set_ind_of_fit(aloc_ind_of_fit);

    for (int i = 0; i < matrix.first.size(); i++) {
        data_fit_pars.get_ind_of_fit()->operator[](i) = {};
    }

    
    pair<double, vector<int>> result_of_fit = fit_steps(N, &data_fit_pars);
    vector<int> final_breaking_points = result_of_fit.second;

    for (int final_breaking_point : final_breaking_points) {
        cout << final_breaking_point << endl;
    }
    
    save_steps(matrix.second,final_breaking_points);
    
}

void print_matrix(const  vector<vector<double>>&matrix){

    int size_y = matrix.size();
    int size_x = matrix[0].size();
    cout << "Matrix (" << size_x << " x " << size_y << ")" << endl;
    if (size_x <10){ 
        for (int i = 0; i< size_y; ++i) {
            for (int j = 0; j < size_x; ++j) {
                cout << "  " << matrix[i][j];
            }
            cout << endl;
        }
    } else {
        for (int i = 0; i< 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                cout << "  " << matrix[i][j];
            }
            cout << endl;
        }
    }
}



int main() {
    const int divisions = 47;

    pair<vector<int>,vector<int>>sizes = read_fold_sizes();

    vector<pair<vector<vector<double>>,vector<vector<double>>>>matrices_folds;
    for (int i = 0; i < 10 ; ++i) {
        vector<vector<double>> train_matrix = read_the_matrix_from_csv("results/trainLossMatrix"+to_string(i)+".csv", divisions);
        print_matrix(train_matrix);
        vector<vector<double>> test_matrix = read_the_matrix_from_csv("results/testLossMatrix"+to_string(i)+".csv", divisions);
        pair<vector<vector<double>>, vector<vector<double>>> matrices = make_pair(train_matrix, test_matrix);
        matrices_folds.push_back(matrices);
    }
    cout << "reading matrices done" << endl;

    int max_number_of_steps = divisions - 1; // number of steps
    int number_of_folds = 10; // number of folds
    int ideal_number_of_steps =  best_number_of_steps(max_number_of_steps,number_of_folds,matrices_folds);
    cout << "best number of steps: done" << endl;
    
    pair<vector<vector<double>>,vector<vector<double>>>matrices;
    matrices.first = read_the_matrix_from_csv("results/dataLossMatrix.csv", divisions);
    matrices.second = read_the_matrix_from_csv("results/dataPeakMatrix.csv", divisions);
    fitting_ideal_steps(ideal_number_of_steps, matrices);
    
    

    //vector<int>breaking_points = fitting_ideal_steps(number_of_steps,points);
    cout << endl;
    cout << endl;
    
    return 0;
}
