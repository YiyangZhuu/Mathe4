# 2.Programmieruebung #
# Testat am 18,06,2024  11:30-11:45  bei Tamme Claus
# explicit math. Form see skript

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from mpl_toolkits.mplot3d import Axes3D

# Parameters
n = 60
h = 1 / n
beta = 5 * np.pi / 6
epsilon_values = [1, 1e-2, 1e-4]

# Grid setup
x = np.linspace(0, 1, n + 1)
y = np.linspace(0, 1, n + 1)
X, Y = np.meshgrid(x, y)

# Function handles for konvektions-Diffusionsproblem
f_kd = lambda x, y: 1  
g_kd = lambda x, y: 0

# Function handles for Poisson
f_poisson = lambda x, y: -4  
g_poisson = lambda x, y: x**2 + y**2  


# Function definitions 
def create_A_b(f, g, epsilon, _upwind = False):
    _n = n-1
    m = _n*_n
    A = lil_matrix((m,m))
    b = np.zeros(m)

    # A #
    #apply zentrale Diffenzenquotient for x,y-derivative
    if(_upwind == False):
        temp1 = -epsilon/(h*h)
        temp2 = np.cos(beta)/(2*h)
        c1 = temp1-temp2
        c2 = temp1+temp2

        #put T in the matrix
        for i in range(_n):
            A[i*_n,i*_n] = -4*temp1
            A[i*_n,i*_n+1] = c2
            for j in range(1,_n-1):
                A[i*_n+j,i*_n+j-1] = c1
                A[i*_n+j,i*_n+j] = -4*temp1
                A[i*_n+j,i*_n+j+1] = c2
            A[i*_n+_n-1,i*_n+_n-2] = c1
            A[i*_n+_n-1,i*_n+_n-1] = -4*temp1
        #put I in the matrix
        for i in range(m-_n):
            A[i,i+_n] = c2
            A[i+_n,i] = c1

    #apply upwind for x,y-derivative
    else:
        temp1 = -epsilon/(h*h)
        temp2 = np.cos(beta)/h
        c1 = temp1-temp2
        c2 = temp1+temp2
        if(np.cos(beta) >= 0):
            #put T in the matrix
            for i in range(_n):
                A[i*_n,i*_n] = -4*temp1 + 2*temp2
                A[i*_n,i*_n+1] = temp1
                for j in range(1,_n-1):
                    A[i*_n+j,i*_n+j-1] = c1
                    A[i*_n+j,i*_n+j] = -4*temp1 + 2*temp2
                    A[i*_n+j,i*_n+j+1] = temp1
                A[i*_n+_n-1,i*_n+_n-2] = c1
                A[i*_n+_n-1,i*_n+_n-1] = -4*temp1 + 2*temp2
            #put I in the matrix
            for i in range(m-_n):
                A[i,i+_n] = temp1
                A[i+_n,i] = c1
        else:
            #put T in the matrix 
            for i in range(_n):
                A[i*_n,i*_n] = -4*temp1 - 2*temp2
                A[i*_n,i*_n+1] = c2
                for j in range(1,_n-1):
                    A[i*_n+j,i*_n+j-1] = temp1
                    A[i*_n+j,i*_n+j] = -4*temp1 - 2*temp2
                    A[i*_n+j,i*_n+j+1] = c2
                A[i*_n+_n-1,i*_n+_n-2] = temp1
                A[i*_n+_n-1,i*_n+_n-1] = -4*temp1 - 2*temp2
            #put I in the matrix
            for i in range(m-_n):
                A[i,i+_n] = c2
                A[i+_n,i] = temp1

    # b #
    c = 1/(h*h)
    # b1 #
    b[0] = f(h,h) + c * (g(h,0)+g(0,h))
    for i in range(1,_n-1):
        b[i] = f((i+1)*h,h) + c * g((i+1)*h,0)
    b[_n-1] = f(1-h,h) + c * (g(1-h,0)+g(1,h))
    # bj #
    for j in range(1,_n-1):
        b[j*_n] = f(h,j*h) + c * g(0,j*h)
        for i in range(1,_n-1):
            b[j*_n+i] = f((i+1)*h,j*h)
        b[j*_n+_n-1] = f(1-h,j*h) + c * g(1,j*h)
    # bn-1 #
    b[m-1-_n] = f(h,1-h) + c * (g(h,1)+g(0,1-h))
    for i in range(1,_n-1):
        b[m-1-_n+i] = f((i+1)*h,1-h) + c * g((i+1)*h,1)
    b[m-1] = f(1-h,1-h) + c * (g(1-h,1)+g(1,1-h))
        
    return A, b

def Aufgabe_a():
    for epsilon in epsilon_values:
        A, b = create_A_b(f_kd, g_kd, epsilon, _upwind = True)
        A = A.tocsr()

        u = spsolve(A, b)
        U = u.reshape((n-1,n-1))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X[1:-1, 1:-1], Y[1:-1, 1:-1], U, cmap='viridis')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('U')
        ax.set_title(f'Convection-Diffusion Solution, ε = {epsilon}')
        plt.show()
    return 0

def Aufgabe_b():
    return

def Aufgabe_c():
    return

if __name__=="__main__":
    Aufgabe_a()