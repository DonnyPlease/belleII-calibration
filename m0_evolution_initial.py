from errno import EROFS
import matplotlib.pyplot as plt
import csv;
import pandas as pd;
import numpy as np;
from datetime import datetime

#import times from the result of the fit of the step function
times = pd.read_csv("results/times.txt")    #import times from csv file
times_array_t = times.T                 #transpose the pandas dataframe
times_array = times_array_t.to_numpy()  #convert it to numpy array



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



#import m0s and their error and convert them to numpy 2D array
errorsmatrixpd = pd.read_csv("results/stdErrorMatrix.csv", header = None)
dataM0pd = pd.read_csv("results/dataPeakMatrix.csv", header = None)
errorsmatrix = errorsmatrixpd.to_numpy()
dataM0 = dataM0pd.to_numpy()



#extract the diagonal elements to the 1D array
N = times_array.size-1
errors = np.zeros(N+1)
m0s = np.zeros(N+1)
for i in range(N):
    errors[i] = (errorsmatrix[i,i])
    m0s[i] = (dataM0[i,i])
#append the last element once again so we can use the method plt.step()
errors[N] = (errorsmatrix[N-1,N-1])
m0s[N] = (dataM0[N-1,N-1])



#create ticks for x axis 
Xticks = np.linspace(times_array[0][0],times_array[0][-1],num = 10)
Xticks_normal = np.vectorize(lambda t:  datetime.utcfromtimestamp((t + 9) * 3600))(Xticks)   #convert unix time in hours to classic date and time



#plot lines that connect start of pause with the function value of the previous interval and end of pause with the function value of the next interval
for i in range(p_start_normal.size):
    plt.plot([p_start_normal[i], p_end_normal[i]],[m0s[i],m0s[i+1]], color = 'lightsteelblue', linewidth = 1.5)



#plot intervals with the function values and their error
for i in range(i_s_normal.size):
    plt.plot([i_s_normal[i],i_e_normal[i]], [m0s[i],m0s[i]], color = "red", linewidth = 2)
    plt.fill_between([i_s_normal[i],i_e_normal[i]], m0s[i]-errors[i], m0s[i]+errors[i], step= None, color='k', alpha=0.15)    #plot the errorsas grey rectangles   



#modify plot properties
plt.title("")
plt.xlabel('Date (and time)', fontsize=18)                  #set x-axis label
plt.ylabel('$E_{cms}$ [MeV]', fontsize=18)                  #set y-axis label
plt.xticks(Xticks_normal, rotation =35)                     #set x ticks we prepared earlier
plt.grid(True,"both","both",)
plt.show()      #show plot   




                                                                                                                       
