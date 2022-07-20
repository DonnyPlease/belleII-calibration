#include <iostream>
#include <tuple>
#include <vector>
#include <filesystem>
#include <numeric>
#include <iomanip>
#include <tuple>
#include <fstream>

#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <Math/Functor.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/DistFunc.h>
#include <Eigen/Dense>

#include <InvariantMassMuMuStandAlone.h>
#include <ChebFitter.h>

using namespace std;
using namespace Eigen;

using Belle2::Pars;
using Belle2::InvariantMassMuMuCalib::Event;
using Belle2::InvariantMassMuMuCalib::getEvents;
using Belle2::InvariantMassMuMuCalib::getInvMassPars;
using Belle2::InvariantMassMuMuCalib::getLikelihood;
using Belle2::InvariantMassMuMuCalib::mainFunction;
using Belle2::InvariantMassMuMuCalib::readEvents;


//This struct is used to make the work with events easier.
struct MyEvent
{
    double minv; // invariant mass
    double t;    // time
    int id;
};


//This struct is used to store information about a job - especially which parts of matrices are calculated in the particular job.
struct JobInfo
{
    int beginRow = 0;
    int beginCol = 0;
    int endRow   = 0;
    int endCol   = 0;
    int split    = 0;
    int size     = 0;
    int jobN     = 0;
};


// comapares times of two structs Event, it is needed for sorting data
bool compareTimeOfEvents(MyEvent e1, MyEvent e2)
{
    return (e1.t < e2.t);
}


//compares times with some other value - used for cutting data.
bool compareTimeOfEventsWithVal(MyEvent ev, double val)
{
    return (ev.t < val);
}


//Prints vector in line (up to 10 first elements) - used for debugging
void print(vector<double> vec)
{
    double bound = min(10, (int)vec.size());
    for (int i = 0; i < bound; ++i)
        cout << vec[i] << ' ';
    cout << endl;
}

/// plots the result of the fit to the Mmumu, i.e. data and the fitted curve, the base function
  static void plotMuMuFitBase(TH1D* hData, TGraph* gr, TH1D* hPull, Pars pars, Eigen::MatrixXd mat, int time)
  {
    bool isBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);

    gStyle->SetOptStat(0);

    TCanvas* can = new TCanvas(Form("canMuMu_%d", time), "");

    TPad* pad1 = new TPad(Form("pad1_%d", time), "", 0, 0.3, 1, 1.0);
    TPad* pad2 = new TPad(Form("pad2_%d", time), "", 0, 0,   1, 0.3);

    pad1->SetBottomMargin(0.05);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);

    pad1->Draw();
    pad2->Draw();

    ///////////////////
    // Main plot
    ///////////////////

    pad1->cd();

    hData->SetMarkerStyle(kFullCircle);
    hData->Draw();
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->Draw("same");
    hData->GetXaxis()->SetLabelSize(0.0001);
    hData->GetYaxis()->SetLabelSize(0.05);
    hData->GetYaxis()->SetTitle("Number of events");
    hData->GetYaxis()->SetTitleSize(0.05);
    hData->GetYaxis()->SetTitleOffset(0.9);

    double mY4S = 10579.4;
    double y = gr->Eval(mY4S);
    TLine* line = new TLine(mY4S, 0, mY4S, y);
    line->SetLineColor(kGreen);
    line->SetLineWidth(2);
    line->Draw();



    TLegend* leg = new TLegend(.15, .4, .35, .87);
    int i = 0, nPars = 0;
    for (auto p : pars) {
      double err = sqrt(mat(i, i));
      if (err != 0) {
        int nDig = log10(p.second / err) + 2;

        TString s   = "%s = %." + TString(Form("%d", nDig)) + "g";
        TString dig = "%." + TString(Form("%d", nDig)) + "g";
        TString digE = "%.2g";
        leg->AddEntry((TObject*)0, Form("%s = " + dig + " #pm " + digE, p.first.c_str(), p.second, err), "h");
        ++nPars;
      }
      ++i;
    }
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();


    double chi2 = 0;
    for (int j = 1; j <= hPull->GetNbinsX(); ++j)
      chi2 += pow(hPull->GetBinContent(j), 2);
    int ndf = hPull->GetNbinsX() - nPars - 1;


    TLegend* leg2 = new TLegend(.73, .75, .93, .87);
    leg2->AddEntry((TObject*)0, Form("chi2/ndf = %.2f", chi2 / ndf), "h");
    leg2->AddEntry((TObject*)0, Form("p = %.2g",   TMath::Prob(chi2, ndf)), "h");

    std::ofstream outfile;
    outfile.open("results/chiSquared.txt", std::ios_base::app); // append instead of overwrite
    outfile << time << "," << (chi2 / ndf) << endl;


    leg2->SetTextSize(0.05);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();


    double mFit = pars.at("m0");
    double yF = gr->Eval(mFit);
    TLine* lineR = new TLine(mFit, 0, mFit, yF);
    lineR->SetLineColor(kRed);
    lineR->SetLineWidth(2);
    lineR->Draw();



    ///////////////////
    // Ratio plot
    ///////////////////

    pad2->cd();
    hPull->SetMarkerStyle(kFullCircle);
    hPull->Draw("p");

    hPull->GetXaxis()->SetTitle("M (#mu#mu) [MeV]");
    hPull->GetYaxis()->SetTitle("pull");
    hPull->GetXaxis()->SetTitleSize(0.13);
    hPull->GetXaxis()->SetTitleOffset(1.25);
    hPull->GetXaxis()->SetLabelSize(0.13);
    hPull->GetXaxis()->SetLabelOffset(0.05);
    hPull->GetXaxis()->SetTickSize(0.07);


    hPull->GetYaxis()->SetTitleSize(0.13);
    hPull->GetYaxis()->SetLabelSize(0.13);
    hPull->GetYaxis()->SetTitleOffset(0.2);
    hPull->GetYaxis()->CenterTitle();


    hPull->GetYaxis()->SetNdivisions(404);

    hPull->GetYaxis()->SetRangeUser(-5, 5);

    TGraph* grLine = new TGraph(2);
    grLine->SetPoint(0, hPull->GetBinLowEdge(1), 0);
    grLine->SetPoint(1, hPull->GetBinLowEdge(hPull->GetNbinsX()) + hPull->GetBinWidth(hPull->GetNbinsX()), 0);
    grLine->SetLineWidth(2);
    grLine->SetLineColor(kRed);
    grLine->Draw("same");


    can->SaveAs(Form("PDFs/mumu%d.pdf", time));


    delete leg;
    delete leg2;
    delete line;
    delete lineR;
    delete grLine;

    delete pad1;
    delete pad2;
    delete can;


    gROOT->SetBatch(isBatch);
  }


/// plots the result of the fit to the Mmumu, i.e. data and the fitted curve
static void plotMuMuFit(const vector<double>& data, const Pars& pars, Eigen::MatrixXd mat, double mMin, double mMax, int time)
  {
    const int nBins = 200;

    // Fill the data histogram
    TH1D::SetDefaultSumw2();
    TH1D* hData = new TH1D("hData", "", nBins, mMin, mMax);
    TH1D* hFit  = new TH1D("hFit", "", nBins, mMin, mMax);
    TH1D* hPull = new TH1D("hPull", "", nBins, mMin, mMax);
    hData->SetDirectory(nullptr);
    hFit->SetDirectory(nullptr);
    hPull->SetDirectory(nullptr);

    // fill histogram with data
    for (auto d : data)
      hData->Fill(d);


    // construct the fitted function
    TGraph* gr = new TGraph();
    const double step = (mMax - mMin) / (nBins);

    for (int i = 0; i <= 2 * nBins; ++i) {
      double m = mMin + 0.5 * step * i;
      double V = mainFunction(m, pars);
      gr->SetPoint(gr->GetN(), m, V);
    }


    // Calculate integrals of the fitted function within each bin
    for (int i = 0; i < nBins; ++i) {
      double lV = gr->GetPointY(2 * i + 0);
      double cV = gr->GetPointY(2 * i + 1);
      double rV = gr->GetPointY(2 * i + 2);

      double I = step / 6 * (lV + 4 * cV + rV);
      hFit->SetBinContent(i + 1, I);
    }

    //Normalization factor
    double F = hData->Integral() / hFit->Integral();

    hFit->Scale(F);

    // Normalize the curve
    for (int i = 0; i < gr->GetN(); ++i)
      gr->SetPointY(i, gr->GetPointY(i) * F * step);


    // calculate pulls
    for (int i = 1; i <= nBins; ++i) {
      double pull = (hData->GetBinContent(i) - hFit->GetBinContent(i)) / sqrt(hFit->GetBinContent(i));
      hPull->SetBinContent(i, pull);
    }


    plotMuMuFitBase(hData, gr, hPull, pars, mat, time);

    delete hData;
    delete hFit;
    delete hPull;
    delete gr;
  }


void histogramPDF(Pars pars, Eigen::MatrixXd mat, const vector<double> & data, int jobN, int jobI)
{   

    int time = stoi(to_string(jobN) + "000" + to_string(jobI));

    plotMuMuFit(data, pars, mat, 10.2e3, 10.8e3, time);

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




// This finction creates pair of means and times after dividing data to n(input) parts. The vector of means has size n, and corresponds to mean of every part.
// The vector of times has size (n-1) and contains times where the data is split.
vector<double> times_of_n_divisions(const vector<MyEvent> &data, long long n)
{
    long long size = data.size(); // get size of data

    vector<double> vector_of_time_slices; // declare second part of the ouput

    long long iSlice = 0; // declare index of slice

    for (int eventNum = 0; eventNum < data.size(); ++eventNum)
    { // go through all the data
        if (eventNum < (((iSlice + 1) * size) / n))
            continue;
        if (iSlice < n - 1)
        {                                                                   // check if we made enough slices already
            double slice_t = (data[eventNum - 1].t + data[eventNum].t) / 2; // store time of slice
            vector_of_time_slices.push_back(slice_t);                       // push the time of slice to the output vector
        }
        ++iSlice; // push back the first "im" of the new slice
    }
    return vector_of_time_slices;
}

pair<vector<double>, int> readTimes()
{
    string x1, x2;
    int N;
    vector<double> vectorOfTimeSlices;

    ifstream dataObj("timeSplits.txt", ios::in); //opening the file.
    if (dataObj.is_open()) //if the file is open
    {
        getline(dataObj, x1, '\n'); //new line after
        N = stoi(x1);
        getline(dataObj, x1, ',');
        vectorOfTimeSlices.push_back(stod(x1));
        for (int i = 0; i < (N-1); ++i) {
            getline(dataObj, x2, '\n');
            getline(dataObj, x1, ','); //new line after
            double time = (stod(x2) + stod(x1)) / 2.0;
            vectorOfTimeSlices.push_back(time);
        }
        getline(dataObj, x1, '\n');
        vectorOfTimeSlices.push_back(stod(x1));
    }
    else cout << "Unable to open file: timeSplits.txt"; //if the file is not open output

    return make_pair(vectorOfTimeSlices, N);
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


//Name explains what the function does.
vector<double> get_data_slice(double t_begin, double t_end, const vector<MyEvent> &data)
{
    //get indices of starting and ending point of the slice we are extracting
    auto startIndex = lower_bound(data.begin(), data.end(), t_begin, compareTimeOfEventsWithVal);
    auto endIndex = lower_bound(startIndex, data.end(), t_end, compareTimeOfEventsWithVal);

    vector<MyEvent> slice(startIndex, endIndex); //perform slice
    vector<double> sliceMinv = {};

    //from the slice we neet to extract minv and put those values in vector
    for (MyEvent &e : slice) {
        sliceMinv.push_back(e.minv);
    }
    return sliceMinv;
}


//
pair<double, double> performFit(const vector<double> & time_divisions, 
                                const vector<MyEvent> & trainData, 
                                const vector<MyEvent> & testData, 
                                int row, int col,
                                JobInfo jobInfo, int jobI)
{

    double t_begin = time_divisions[row];
    double t_end = time_divisions[col + 1];

    vector<double> trainDataSlice = get_data_slice(t_begin, t_end, trainData);
    vector<double> testDataSlice = get_data_slice(t_begin, t_end, testData);

    // Perform maximum likelihood fit

    auto start_time = std::chrono::high_resolution_clock::now();

    Pars pars, pars0;
    MatrixXd mat;

    cout << "fitting..." << endl;

    tie(pars, mat) = getInvMassPars(trainDataSlice, pars0);//Fit

    histogramPDF(pars, mat, trainDataSlice, jobInfo.jobN, jobI);

    auto end_time = chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    cout << "Fitting time:   " << time/chrono::milliseconds(1) << endl;


    start_time = std::chrono::high_resolution_clock::now();
    // Get the - log likelihood, where "pars" define the model, and data are data points
    double trainLoss = getLikelihood(trainDataSlice, pars);
    // Get the - log likelihood, where "pars" define the model, and data are data points
    double testLoss = getLikelihood(testDataSlice, pars);

    end_time = chrono::high_resolution_clock::now();
    time = end_time - start_time;
    cout << "-Log likelihood time :    " << time/chrono::milliseconds(1) << endl;
    
    return make_pair(trainLoss, testLoss);
}


//Is called when a particular job works on some part of one of the 10 train-test sets of matrices.
pair<vector<double>,vector<double>> calculateTrainTestJob(const vector<MyEvent> & data, const vector<double> & time_divisions, const JobInfo & jobInfo)
{
    const int size = time_divisions.size() - 1; //size of one line in matrix
    vector<MyEvent> trainData, testData; //Declare.
    vector<double> trainLossVector, testLossVector; //Declare.

    tie(trainData, testData) = train_test_split(data, jobInfo.split); // Split data.

    int row = jobInfo.beginRow; // Save the value(s) for readability.
    int col = jobInfo.beginCol;
    
    for (int i = 0; i < jobInfo.size; ++i)
    {
        double trainLoss, testLoss;
        tie(trainLoss, testLoss) = performFit(time_divisions, trainData, testData, row, col, jobInfo, i);
        trainLossVector.push_back(trainLoss); //only for testing
        testLossVector.push_back(testLoss); //only for testing

        cout << ((100.0*i)/jobInfo.size) << "% done in job: " << jobInfo.jobN << endl;

        if (row == (size-1) && col == (size-1)) break; //If we reached end of matrix - brake.
        
        //If the last condition didn't end the cycle, increment col and if it is too big, increment row and set col so that we are back on the diagonal.
        ++col;
        if (col < size) continue;

        ++row;
        col = row;


    }

    return make_pair(trainLossVector, testLossVector);
}

//
void save_matrix_to_csv(const MatrixXd &matrix, string name)
{

    IOFormat Precision(FullPrecision, DontAlignCols, ",", "", "", ",\n");
    ofstream file("matrices/" + name);
    file << matrix.format(Precision) << endl;
    file.close();
}


//This function is called when a job works on the 11th set of matrices - principle is the same as in calculation of train-test sets, but there are some differences - no split, more matrices in the set...
tuple<vector<double>, vector<double>, vector<double>> calculateOverallLoss(const vector<MyEvent> &data, const vector<double> &time_divisions, JobInfo jobInfo)
{
    int size = time_divisions.size() - 1;
    vector<double> peaks, stdErrors, loss;
    int row = jobInfo.beginRow;
    int col = jobInfo.beginCol;

    for (int i = 0; i < jobInfo.size; ++i)
    {
        double t_begin = time_divisions[row];
        double t_end = time_divisions[col+1];

        vector<double> data_slice = get_data_slice(t_begin, t_end, data);

        // Perform maximum likelihood fit
        Pars pars, pars0;
        MatrixXd mat;

        auto start_time = chrono::high_resolution_clock::now();
        tie(pars, mat) = getInvMassPars(data_slice, pars0);
        
        histogramPDF(pars, mat, data_slice, jobInfo.jobN, i);

        auto end_time = chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        cout << "Fitting time:   " << time/chrono::milliseconds(1) << endl;


        start_time = chrono::high_resolution_clock::now();
        // Get the - log likelihood, where "pars" define the model, and data are data points
        double logL = getLikelihood(data_slice, pars);
        end_time = chrono::high_resolution_clock::now();
        time = end_time - start_time;
        cout << "-log likelihood and others time:   " << time/chrono::milliseconds(1) << endl;

        loss.push_back(logL);
        peaks.push_back(pars.at("m0"));
        stdErrors.push_back(sqrt(mat(4,4)));
        
        ++col;
        if (col < size) continue;

        ++row;
        col = row;
    }
    return make_tuple(peaks, stdErrors, loss);
}

//Save vector of result with each element in its own new line.
void saveVectorOfResult(vector<double> vec, string name)
{
    ofstream file("results/" + name + ".txt");
    for(double d : vec) {
        file << setprecision(12) << d << endl;
    }
    file.close();
}


JobInfo initJob(const int job, const int numOfJobs, const int size) 
{   

    int sizeOfMatrices = size; // eventually this number should be 1000

    int numOfJobsForOneSet = numOfJobs / 11; // we divide the overall number of jobs with 11 so we know how many jobs correspond to one set
    int jobNumberForThisSet = job % numOfJobsForOneSet; //we want to know which part of the set of matrices is this job working on.

    int numberOfElementsOnOrAboveTheDiagonal = (sizeOfMatrices*(sizeOfMatrices+1)) / 2; //Formula comes from the sum of aritmetic sequence.
    int numberOfElementsForOneJob = (numberOfElementsOnOrAboveTheDiagonal + numOfJobsForOneSet - 1) / numOfJobsForOneSet; //this weird equation (x+y-1)/y guarantees that division x/y is rounded up
    

   
    JobInfo result;
    result.split = job / numOfJobsForOneSet;
    result.jobN = job;

    int elementIndex = 0;
    int indexOfStart = jobNumberForThisSet*numberOfElementsForOneJob; //included
    int indexOfEnding = (jobNumberForThisSet+1)*numberOfElementsForOneJob - 1; //included
    result.size = numberOfElementsForOneJob;

    if (indexOfStart >= numberOfElementsOnOrAboveTheDiagonal) {
        result.size = 0;
        return result;
    }

    //For each element above the diagonal - chceck whether we start there or end there, or nothing.
    for (int row = 0; row < sizeOfMatrices; ++row) {
        for (int col = row; col < sizeOfMatrices; ++col) {
            if (elementIndex == indexOfStart) {
                result.beginRow = row;
                result.beginCol = col;
            } else if (elementIndex == indexOfEnding) {
                result.endRow = row;
                result.endCol = col;
                return result;
            }
            ++elementIndex;
        }
    }
    //This happens only when the last job has less elements than the rest.
    result.endRow = sizeOfMatrices - 1;
    result.endCol = sizeOfMatrices - 1;
    return result;
}

void clearChiSquared()
{
    ofstream file("results/chiSquared.txt");
    file.close();
}

void saveVectorOfTimes(vector<double> times)
{
    ofstream file("results/times.txt");
    file << "times" << endl;
    for(double d : times) {
        file << setprecision(12) << d << endl;
    }
    file.close();
}

//
void create_and_save_matrices(int job, int numOfJobs, int N)
{
    vector<MyEvent> data = read_and_sort(); //Read and sort data.
    cout << "Reading and sorting: Done" << endl; //Send the message.
    //vector<MyEvent> data(all_data.begin(), all_data.begin() + 10000);
    vector<double> times = times_of_n_divisions(data, N); //Divide data by the times. Size is (N-1)

    times.insert(times.begin(), data[0].t); //Add time of the first element in data to the beginning.
    times.push_back(data.back().t + 1); //Add time that is bigger that the time of the last element to the end.
    //This two lines increas the size of the vector of times to (N+1)

    tie(times,N) = readTimes();
    cout << N << endl;
    clearChiSquared();
    

    JobInfo jobInfo = initJob(job, numOfJobs, N); //Initialize the struct contatining info about the job.
 
    if (jobInfo.split < 10) { //if the job is calculating matrices for the fold.
        vector<double> trainLossVector, testLossVector; //Declare/
        tie(trainLossVector, testLossVector) = calculateTrainTestJob(data, times, jobInfo); //Calculate the elements corresponding to this job and save them in vectors.
        
        saveVectorOfResult(trainLossVector, "trainLossVector" + to_string(job)); //Save the vector of train loss to file.
        saveVectorOfResult(testLossVector, "testLossVector" + to_string(job)); //Save the vector of test loss to file.

    } 
    else { // if the job is calculating matrices of overall loss, etc.
        vector<double> peaks, stdErrors, loss;
        tie(peaks, stdErrors, loss) = calculateOverallLoss(data, times, jobInfo);
        
        saveVectorOfResult(peaks, "peaks" + to_string(job));
        saveVectorOfResult(stdErrors, "stdErrors" + to_string(job));
        saveVectorOfResult(loss, "loss" + to_string(job));
    }
    
    saveVectorOfTimes(times); //Save times of divisions.
}


//Main function - receives two arguments.
int main(int argc, char *argv[])
{
    int job = stoi(argv[1]);
    int numOfJobs = stoi(argv[2]);
    int N = stoi(argv[3]);
    create_and_save_matrices(job-1, numOfJobs, N);

    return 0;
}
