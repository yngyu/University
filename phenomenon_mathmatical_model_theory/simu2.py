import heapq
import queue

import numpy as np
import matplotlib.pyplot as plt

class MM1:
    mu = 1.0 # サービス確率
    Tstart = 1000000 # 記録開始タイムステップ
    Tend   = 10000000 # 記録終了タイムステップ

    come = True
    go = False

    experiment_rho = []
    experiment_Lambds = [0.05 * i for i in range(1, 22)]
    experiment_Ns = []
    experiment_Ws = []

    def __init__(self):
        self.experiment_rho.clear()
        self.experiment_Lambds = [0.05 * i for i in range(1, 22)]
        self.experiment_Ns.clear()
        self.experiment_Ws.clear()

    def simulate(self):
        self.__init__()
        for Lambda in self.experiment_Lambds:
            n = 0 # 滞在人数
            nSum = 0 # 滞在人数の累積
            wSum = 0 # 来た人数が待った累計
            allmembers = 0 # 来た人数の累計
            t = 0
            t_before = -1

            come_count = 0
            go_count = 0

            pq = [] # イベント行列
            heapq.heappush(pq, (0, self.come))
            come_count += 1

            que = queue.Queue() # 待ち行列に入っている人間

            while True:
                top = heapq.heappop(pq)
                t = top[0]
                flag = top[1]

                # Nsについて記録
                if t >= self.Tstart:
                    # 終わりの時間が来ていたら
                    if t > self.Tend:
                        nSum += n*(self.Tend - t_before)
                        break
                    else:
                        # before がまだスタート前だったら
                        if t_before < self.Tstart:
                            nSum += n*(t - self.Tstart)
                        else:
                            nSum += n*(t - t_before)

                if flag == self.come:
                    n += 1
                    come_count -= 1
                    que.put(t)
                    # Wsについて記録
                    if t >= self.Tstart and t < self.Tend:
                        allmembers += 1
                else:
                    n -= 1
                    go_count -= 1
                    come_time = que.get()
                    # Ws について記録
                    if t >= self.Tstart and t < self.Tend:
                        wSum += t - come_time

                t_before = t

                # 次のイベントを計算
                if come_count == 0:
                    heapq.heappush(pq, (t + self.calculate_next_come(Lambda), self.come))
                    come_count += 1
                if n > 0 and go_count == 0:
                    heapq.heappush(pq, (t + self.calculate_next_go(self.mu), self.go))
                    go_count += 1

            rho = Lambda/self.mu
            self.experiment_rho.append(rho)
            self.experiment_Ns.append(nSum/(self.Tend - self.Tstart))
            self.experiment_Ws.append(wSum/allmembers)
    
        return self.experiment_Lambds, self.experiment_rho, self.experiment_Ns, self.experiment_Ws
        
    def calculate_next_come(self, Lambda):
        return np.random.exponential(1.0/Lambda)
    
    def calculate_next_go(self, mu):
        return np.random.exponential(1.0/mu)

    def calculate_theoritical_result(self):
        Lambda = np.linspace(0, 0.99, 9900)
        rho = Lambda/self.mu
        Ns = rho/(1 - rho)
        return rho, Ns
    
    def plot(self, ax):
        experiment_Lambds, experiment_rho, experiment_Ns, experiment_Ws = self.simulate()
        theoritical_rho, theoritical_Ns = self.calculate_theoritical_result()
        ax.plot(theoritical_rho, theoritical_Ns, c = "Black")
        ax.scatter(experiment_rho, experiment_Ns, c = "Green", marker="s", label = "S.D. = 1")

class MG1(MM1):
    come_loc = 1.0
    come_scale = 0.0

    def calculate_next_go(self, Lambda):
        return np.random.normal(self.come_loc, self.come_scale)
    
    def calculate_theoritical_result(self):
        Lambda = np.linspace(0, 0.99, 9900)
        rho = Lambda/self.mu
        Cs = self.come_scale/self.come_loc
        Wq = (1 + Cs**2)/2*rho/(1 - rho)*self.come_loc # ポラチェック・ヒンチンの公式
        Ns = (Wq + self.come_loc)*Lambda
        return rho, Ns
    
    def plot(self, ax):
        experiment_Lambds, experiment_rho, experiment_Ns, experiment_Ws = self.simulate()
        theoritical_rho, theoritical_Ns = self.calculate_theoritical_result()
        ax.plot(theoritical_rho, theoritical_Ns, c = "Black")
        ax.scatter(experiment_rho, experiment_Ns, c = "Blue", marker="x", label = "S.D. = 0")

class MG3(MG1):
    come_scale = 3.0

    def calculate_next_go(self, Lambda):
        theta = self.come_scale**2/self.come_loc
        k = self.come_loc**2/self.come_scale**2

        return np.random.gamma(k, theta)
    
    def plot(self, ax):
        experiment_Lambds, experiment_rho, experiment_Ns, experiment_Ws = self.simulate()
        theoritical_rho, theoritical_Ns = self.calculate_theoritical_result()
        ax.plot(theoritical_rho, theoritical_Ns, c = "Black")
        ax.scatter(experiment_rho, experiment_Ns, c = "Red", marker="o", label = "S.D. = 3")

fig = plt.figure()
ax = fig.add_subplot(111)
plt.ylim(0, 6)

mm1 = MM1()
mm1.plot(ax)

mg1 = MG1()
mg1.plot(ax)

mg3 = MG3()
mg3.plot(ax)

plt.legend(loc = "upper left")
plt.xlabel(r"$\rho$")
plt.ylabel(r"$N_s$")
plt.savefig("result2_2.png", bbox_inches='tight')
plt.show()
