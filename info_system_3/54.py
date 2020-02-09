from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
import csv

temperature = []
with open("TokyoTemperatureSince2001.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        temperature.append(float(row[1]))

temperature = array(temperature)
n = len(temperature)
test_n = int(n*0.7)
month = arange(n)

model_number = 20

results = empty(model_number)

for m in range(1,model_number+1):
    x = empty((test_n-m,m))
    y = empty(test_n-m)
    for i in range(test_n-m):
        y[i] = temperature[m+1+i-1]
        for j in range(m):
            x[i][j] = temperature[m+i-j-1]
    
    alpha = la.solve(dot(x.T,x),dot(x.T,y))
    temp = zeros_like(temperature)
    for i in range(m):
        temp[i] = temperature[i]
    for i in range(m,n):
        for j in range(len(alpha)):
            temp[i] += alpha[j]*temperature[i-j-1]
    if m == 12:
        print(alpha)
        plt.scatter(month,temperature)
        plt.plot(month,temp,color = "C1")
        plt.scatter(month,temp,color = "C1",marker = "x")
        plt.xlabel("month from begin")
        plt.ylabel("temperature")
        plt.savefig("540.png",bbox_inches="tight", pad_inches=0.05)
        plt.show()
    

    results[m-1] = (n-test_n)*log(sum(((temp-temperature)**2)[test_n:])/(n-test_n)) + 2*m

plt.scatter(range(1,model_number+1),results)
plt.xlabel("model degree")
plt.ylabel("AIC")
plt.savefig("541.png",bbox_inches="tight", pad_inches=0.05)
plt.show()
print(results)
print(n)
print(test_n)