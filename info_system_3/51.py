from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
import random

#やりたいこと、まず適当にテストとtrainで分ける、
#次に全部でAICを試す。
year = []
price = []

with open("mazdas.txt") as f:
    lines = f.readlines()
    for line in lines:
        line = line.split()
        if len(line) == 2:
            year.append(line[0])
            price.append(line[1])

year.pop(0)
price.pop(0)

year = [float(each_year) for each_year in year]
price = [float(each_price) for each_price in price]

n = len(year)
year = array(year)
price = array(price)

m_max = 10

flag = False

def only_test_train():
    test_size = int(n*0.3)
    train_size = n - test_size
    test_times = 100000
    Q = zeros(m_max-1)
    for i in range(3,m_max):  
        m = i #何次の近似をするか選択
        temp_Q = 0 #Qに入れる前のぞそれぞれのQの計算場所
        for j in range(test_times): #randomで選ぶのでtest_times回試行
            random_index_list = list(range(n)) #testとtrainをそれぞれ選ぶ
            random_test_index_list = random.sample(random_index_list,test_size)
            train_year = []
            train_price = []
            test_year = []
            test_price = []
            for k in range(n):
                if k in random_test_index_list:
                    test_year.append(year[k])
                    test_price.append(price[k])
                else:
                    train_year.append(year[k])
                    train_price.append(price[k])
            test_year = array(test_year)
            test_price = array(test_price)
            train_year = array(train_year)
            train_price = array(train_price)
            
            X_m = empty((train_size,m+1))
            for k in range(m+1):
                X_m[:, k] = train_year**k
            y = train_price
            alp = la.solve(dot(X_m.T,X_m), dot(X_m.T,y))
            price_from_test_year = zeros_like(test_year)
            for k in range(len(alp)):
                price_from_test_year += alp[k]*test_year**k
            temp_Q += sum((test_price-price_from_test_year)**2)/test_size #まずtestの最小二乗和を入れる。

            if flag == False:
                plt.scatter(train_year,train_price,label = "train")
                plt.scatter(test_year,test_price,label = "test")
                x_v = arange(65,95,1e-2)
                y_v = zeros_like(x_v)
                for k in range(len(alp)):
                    y_v += alp[k]*x_v**k
                plt.plot(x_v,y_v)
                plt.xlabel("year")
                plt.ylabel("price")
                plt.ylim(-5000,45000)
                plt.legend(loc = "lower right")
                plt.savefig("51_ot_3.png",bbox_inches="tight", pad_inches=0.05)
                plt.show()
                return

        temp_Q /= test_times #test_timesで平均を取る。
        Q[i-1] = temp_Q #それぞれ入れていく。
    print(Q)

def only_AIC():
    AIC = zeros(m_max-1)
    for i in range(1,m_max):
        m = i
        X_m = empty((n,m+1))
        for j in range(m+1):
            X_m[:,j] = year**j
        y = price
        alp = la.solve(dot(X_m.T,X_m), dot(X_m.T,y))
        price_from_year = zeros_like(year)
        x_v = arange(65,95,1e-2)
        y_v = zeros_like(x_v)
        for k in range(len(alp)):
            y_v += alp[k]*x_v**k
            price_from_year += alp[k]*year**k
        if i == 8:
            plt.scatter(year,price)
            plt.plot(x_v,y_v)
            plt.ylim(-5000,45000)
            plt.xlabel("year")
            plt.ylabel("price")
            plt.ylim(-5000,45000)
            plt.savefig("51_oa_8.png",bbox_inches="tight", pad_inches=0.05)
            plt.show()

        Q = sum((price - price_from_year)**2)
        AIC[i-1] = n*log(Q/n) + 2*(len(alp))
    print(AIC)

only_AIC()