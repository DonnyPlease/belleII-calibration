import matplotlib.pyplot as plt
import csv;
import pandas as pd;
import numpy as np;
from matplotlib.ticker import ScalarFormatter


#import sizes and calculate loss of case when there is only one interval
test0, train0 = 0, 0;
sizes = pd.read_csv("results/sizes.txt", header = None)
for i in range(10):
    train_loss_matrix = pd.read_csv("results/trainLossMatrix"+str(i)+".csv", header = None).to_numpy()
    test_loss_matrix = pd.read_csv("results/testLossMatrix"+str(i)+".csv", header = None).to_numpy()
    train0 += train_loss_matrix[0][-2]/sizes[0][i]/10
    test0 += test_loss_matrix[0][-2]/sizes[1][i]/10


#import results and put the result of the calculation above in the beginning
data = pd.read_csv("results/result.csv")
data.loc[-1] = (0, train0, test0)
data.index = data.index + 1
data = data.sort_index()



N=8744092
pd.options.display.precision = 8



#plot results
plt.plot(data.n+1, data.train_error*N, label='Train loss')
plt.plot(data.n+1, data.test_error*N, label='Test loss')



#highlight the minimum of test loss
min_index = np.argmin(data.test_error) + 1
plt.scatter(min_index, data.test_error[min_index]*N, color = "red", label = "Minimum in $\it{k}$ = " + str(min_index))



#personalize the plot
plt.xlabel("Number of intervals $\it{k}$",fontsize = 16)
plt.ylabel("Value of loss", fontsize = 16)
plt.title("")
y_formatter = ScalarFormatter(useOffset=False)
ax = plt.axes()
ax.yaxis.set_major_formatter(y_formatter)
plt.xticks(np.linspace(0,50,11))
plt.grid(True,'both','both')
plt.legend()
plt.show()
