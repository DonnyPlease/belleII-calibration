import matplotlib.pyplot as plt
import csv;
import pandas as pd;
import numpy as np;
from datetime import datetime

#get input on which fit we want to import
print("write number of steps you would like to plot")
N = input()



#import results - times and values of m0 with error (called std_error even though it is not accurate)
times = pd.read_csv("results/times.txt")    #import times from csv file
data = pd.read_csv("results/"+str(N) +"steps.csv")    #import data - the format is (means,start(time),std_errors)

#convert the imported pandas frames to numpy
times_array_t = times.T                 #transpose the pandas dataframe
times_array = times_array_t.to_numpy()  #convert it to numpy array
times_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(times_array[0])   #convert unix time in hours to classic date and time
errors = data.std_errors.to_numpy()     #convert standard errors to numpy array
means = data.means.to_numpy()           #convert ...
	


#plot the fitted steps
plt.step(times_normal[data.start.T], data.means, where = "post", label = str(data.means.size-1) + " intervals")    #plot steps of m0 as blue(default) line
plt.fill_between(times_normal[data.start.T], means-errors, means+errors, step='post', color='blue', alpha=0.25)    #plot the errors as grey rectangles  



#import another table of  start and end times of the initial calibration intervals and store them "twice" - borders of intervals and borders of pauses
timePauses = pd.read_csv("results/timeSplits.txt")  
intervals_start = timePauses.start.T.to_numpy()
intervals_end = timePauses.end.T.to_numpy()
pauses_start = timePauses.end.T.to_numpy()[:-1]
pauses_end = timePauses.start.T.to_numpy()[1:]

#perform conversion from unix time for all 4 arrays
p_start_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(pauses_start)   #convert unix time in hours to classic date and time
p_end_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(pauses_end)   #convert unix time in hours to classic date and time
i_s_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(intervals_start)   #convert unix time in hours to classic date and time
i_e_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(intervals_end)   #convert unix time in hours to classic date and time


#import the matrices m0 and its error and convert into numpy
errorsmatrixpd2 = pd.read_csv("results/stdErrorMatrix.csv", header = None)
dataM0pd = pd.read_csv("results/dataPeakMatrix.csv", header = None)
errorsmatrix2 = errorsmatrixpd2.to_numpy()
dataM0 = dataM0pd.to_numpy()



#extract the diagonal elements to the 1D arrays
N = times_array.size-1
errors2 = np.zeros(N+1)
m0s = np.zeros(N+1)
for i in range(N):
    errors2[i] = (errorsmatrix2[i,i])
    m0s[i] = (dataM0[i,i])
#append the last element once again so we can use the method plt.step()
errors2[N] = (errorsmatrix2[N-1,N-1])
m0s[N] = (dataM0[N-1,N-1])




#create ticks for x axis 
Xticks = np.linspace(times_array[0][0],times_array[0][-1],num = 10)
Xticks_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(Xticks)   #convert unix time in hours to classic date and time



#plot lines that connect start of pause with the function value of the previous interval and end of pause with the function value of the next interval
for i in range(p_start_normal.size):
    plt.plot([p_start_normal[i], p_end_normal[i]],[m0s[i],m0s[i+1]], color = 'coral', linewidth = 0.5)



#plot intervals with the function values and their error - the condition inside the loop ensures that there is only one red line in the legend - probably not the most efficient way of doing this
for i in range(i_s_normal.size):
    if (i == 0):
        plt.plot([i_s_normal[i],i_e_normal[i]], [m0s[i],m0s[i]], color = "red", linewidth = 2, label = "47 intervals")
    else:
        plt.plot([i_s_normal[i],i_e_normal[i]], [m0s[i],m0s[i]], color = "red", linewidth = 2)
    plt.fill_between([i_s_normal[i],i_e_normal[i]], m0s[i]-errors2[i], m0s[i]+errors2[i], step= None, color='red', alpha=0.1)    #plot the errorsas grey rectangles   



#modify plot properties
plt.title("")
plt.xlabel('Date (and time)', fontsize=15)                  #set x-axis label
plt.ylabel('$E_{cms}$ [MeV]', fontsize=15)                  #set y-axis label
plt.xticks(Xticks_normal, rotation =35)                     #set x ticks we prepared earlier
plt.grid(True,"both","both",)
plt.show()      #show plot                                                                                                                       
