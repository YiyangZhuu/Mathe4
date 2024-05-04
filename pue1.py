# 1.Programmieruebung #
# 411927 Nikoly Panov; 429943 Yanliang Zhu; 442567 Shijie Xu; 447931 Yiyang Zhu

import numpy as np
import matplotlib as plt

# Const Parameters
pi = np.pi
N = 4       #Matrix dimension
M = 4
X = np.empty(N)     #Stuetzstellen initialization
Y = np.empty(M)
for j in range(N):
    X[j] = pi*(2*j+1)/(2*N)
for k in range(M):
    Y[k] = pi*(2*k+1)/(2*M)


def Funktion_auswerter(f): #Funktion auswerten an Stuetzstellen x_j und y_k
    F = np.zeros((N,M))
    for j in range(N):
        for k in range(M):
           F[j][k] = f(X[j],Y[k])
    return F

def DCT(F): #F ist die Matrix (F_jk)=f(x_j,y_k) ausgewertet an x_j und y_k
    D = np.empty((N,M))
    c = 4/(N*M)
    for j in range(N):
        for k in range(M):
            temp = 0
            for _j in range(N):
                for _k in range(M):
                    temp += F[j][k] * np.cos(_j*X[_j]) * np.cos(_k*Y[_k])
            D[j][k] = c * temp
    return D #D = (d_jk)

def TDCT(D,x,y):
    return A

# Aufgabe 1
# Testfunktion (a)
def f1(x,y):
    return np.cos(2*x) + np.cos(3*y)

# Testfunktion (b)
def f2(x,y):
    return (x-pi/2)*(x-pi/2) + (y-pi/2)*(y-pi/2)

def Test():
    return

# Aufgabe 2

# Aufgabe 3
def Zickzack():
    return


if __name__ == "__main__":
    F = Funktion_auswerter(f2)
    D = DCT(F)
    print(D)
