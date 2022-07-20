#include <iostream>

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <TTreeReader.h>
#include <TFile.h>
#include <fstream>
#include <tuple>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <Eigen/Dense>
#include <chrono>

using namespace std;
using namespace Eigen;


//This struct is used to make the work with events easier.
struct MyEvent
{
    double minv; // invariant mass
    double t;    // time
    int id;
};


// comapares times of two structs Event, it is needed for sorting data
bool compareTimeOfEvents(MyEvent e1, MyEvent e2)
{
    return (e1.t < e2.t);
}

// Splits data into two vectors of Events with respect to the "id" parameter of the event. As input, the testing "id" must be stated.
pair<vector<MyEvent>, vector<MyEvent>> train_test_split(const vector<MyEvent> &data, int K)
{
    vector<MyEvent> train_set, test_set;
    train_set.reserve(data.size()); //reserve for time save
    test_set.reserve(data.size() / 7); //reserve a little bit more than needed, on average final size of test_set should be 1/10 of size of data
    for (const MyEvent &e : data) {
        if (e.id != K) {
            train_set.push_back(e);
        } else {
            test_set.push_back(e);
        }
    }
    return make_pair(train_set, test_set);
}

//Read and sorts the data.
vector<MyEvent> read_and_sort()
{
    TFile *f = new TFile("data/muons_peak.root");
    TTreeReader reader("tree", f);
    TTreeReaderValue<int> id(reader, "id");
    TTreeReaderValue<double> time(reader, "time");
    TTreeReaderValue<double> minv(reader, "minv");

    vector<MyEvent> data;              // declare vector of strucures Event for data
    data.reserve(reader.GetEntries()); // reserve lenght for time savings

    // Store data in vector called data
    while (reader.Next())
    {
        MyEvent element; // declare temporary element which will be later pushed back repeatedly
        element.minv = (*minv) * 1000;
        element.t = *time;
        element.id = *id;
        data.push_back(element);
    }

    // Sort data with respect to time (t)
    sort(data.begin(), data.end(), compareTimeOfEvents);
    return data;
}

//Generates basic testing data.
vector<MyEvent> testingData() {
    MyEvent e;
    vector<MyEvent> result;
    for (int i = 0; i < 10000; ++i){
        e.minv = i;
        e.t = i;
        e.id = i % 10;
        result.push_back(e);
    }
    return result;
}

void save_sizes_of_splits()
{
    vector<MyEvent> data = read_and_sort();
    ofstream file("results/sizes.txt");
    for (int i = 0; i < 10; ++i){
        pair<vector<MyEvent>, vector<MyEvent>>split = train_test_split(data, i);
        file << split.first.size() << "," << split.second.size() << endl;
    }
    file.close();
}

int main()
{
    save_sizes_of_splits();
    return 0;
}