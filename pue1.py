# 1.Programmieruebung #
# 411927 Nikoly Panov; 429943 Yanliang Zhu; 442567 Shijie Xu; 447931 Yiyang Zhu

import numpy as np
import matplotlib.pyplot as plt

# Const Parameters
pi = np.pi
sigma = np.array([[10,15,25,37,51,66,82,100],        #Datenkomprimierungsfaktor
                 [15,19,28,39,52,67,83,101],
                 [25,28,35,45,58,72,88,105],
                 [51,52,58,66,76,89,103,119],
                 [66,67,72,79,89,101,114,130],
                 [82,83,88,94,103,114,127,142],
                 [100,101,105,111,119,130,142,156]])


def Funktion_auswerter(f,N,M): #Funktion auswerten an Stuetzstellen x_j und y_k
    X = np.empty(N)     #Stuetzstellen initialization
    Y = np.empty(M)
    for j in range(N):
        X[j] = pi*(2*j+1)/(2*N)
    for k in range(M):
        Y[k] = pi*(2*k+1)/(2*M)
    
    F = np.zeros((N,M))
    for j in range(N):
        for k in range(M):
           F[j][k] = f(X[j],Y[k])
    return F,X,Y

def DCT(F,X,Y,N,M): #F ist die Matrix (F_jk)=f(x_j,y_k) ausgewertet an x_j und y_k
    D = np.empty((N,M))
    c = 4/(N*M)
    for j in range(N):
        for k in range(M):
            temp = 0
            for _j in range(N):
                for _k in range(M):
                    temp += F[_j][_k] * np.cos(_j*X[_j]) * np.cos(_k*Y[_k])
            D[j][k] = c * temp
    return D #D = (d_jk)

def TDCT(D,x,y,N,M):     #x,y vector   
    #allgemein TDCT, nicht fuer JPEG geeignet
    #bei JPEG bitte noch D zuerst faktoralisieren
    P = len(x)
    Q = len(y)
    A = np.zeros((P,Q))
    for p in range(P):
        for q in range(Q):
            #_j = _k = 0
            A[p][q] += 0.25 * D[0][0]    # * cos(0*x) * cos(0*y) = 1
            #_j = 0
            for k in range(1,M):
                A[p][q] += 0.5 * D[0][k] * np.cos(k*y[q])    # * cos(0*x) = 1
            #_k = 0
            for j in range(1,N):
                A[p][q] += 0.5 * D[j][0] * np.cos(j*x[p])    # * cos(0*y) = 1 
            #sonst
            for j in range(1,N):
                for k in range(1,M):
                    A[p][q] += D[j][k] * np.cos(j*x[p]) * np.cos(k*y[q])
    return A

# Aufgabe 1
# Testfunktion (a)
def f1(x,y):
    return np.cos(2*x) + np.cos(3*y)

# Testfunktion (b)
def f2(x,y):
    return (x-pi/2)*(x-pi/2) + (y-pi/2)*(y-pi/2)

def Test_A1():
    N,M=4,4
    F1,X,Y = Funktion_auswerter(f1,N,M)
    D1 = DCT(F1,X,Y,N,M)
    A1 = TDCT(D1,X,Y,N,M)
    print(D1)
    print(A1)
    return

def Test_A2():
    Parameter_N = 40
    P=100
    Array_MAX=[]
    Array_N = np.linspace(1,Parameter_N,Parameter_N)  #1:40
    x = np.array([p*pi/P for p in range(1,P)])  #len=P-1
    y = np.array([q*pi/P for q in range(1,P)])  #len=P-1
    
    for N in range(1,Parameter_N+1):   
        F2,X,Y = Funktion_auswerter(f2,N,N)  #MatrixF2 N*N
        D2 = DCT(F2,X,Y,N,N)
        A2 = TDCT(D2,x,y,N,N)   #MatrixA (P-1)*(P-1)
        print(N)
        MAX = 0
        for p in range(1,P):
            for q in range(1,P):
                TEMP = abs(        f2(x[p-1],y[q-1]) - A2[p-1][q-1]        )
                if TEMP >MAX:
                    MAX = TEMP
        Array_MAX.append(MAX)
    
    plt.plot(Array_N, Array_MAX)
    plt.title('Max_Fehler(N)')
    plt.xlabel('N')
    plt.ylabel('Fehler')
    plt.show()



# Aufgabe 2
# TODO

# Aufgabe 3
# TODO 

def Zickzack():
    return


if __name__ == "__main__":
    Test_A1()
    Test_A2()