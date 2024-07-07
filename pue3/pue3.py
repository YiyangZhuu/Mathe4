# 3. Programmieruebung 
# Testat am 16,07,2024  15:30-15:45  bei Tamme Claus

import numpy as np
import matplotlib.pyplot as plt

#parameters
N=50

# Aufgabe 1
# Jacobi 
def JC(A, b, delta=1e-4, max_iteration=1000):
    m = len(b)
    u = np.zeros(m)
    u_new = np.zeros(m)
    steps = 0
    error = 1
    L = np.tril(A, k=-1)
    R = np.triu(A, k=1)
    D = np.diag(A)
    D = np.diag(D)
    D_inv = np.linalg.inv(D)
    B = -np.dot(D_inv , L + R)
    y = np.dot(D_inv, b)
    while(error >= delta and steps <= max_iteration):
        u = u_new
        u_new = np.dot(B, u) + y
        error = np.linalg.norm(b - np.dot(A,u),ord=np.inf)
        steps += 1
        print("JC: " + str(steps))
    return u, steps

# Gauss-Seidel
def GS(A,b, delta=1e-4, max_iteration=1000):
    m = len(b)
    u = np.zeros(m)
    u_new = np.zeros(m)
    steps = 0
    error = 1
    L = np.tril(A, k=-1)
    R = np.triu(A, k=1)
    D = np.diag(A)
    D = np.diag(D)
    D_L = D+L
    D_L_inv = np.linalg.inv(D_L)
    B = -np.dot(D_L_inv , R)
    y = np.dot(D_L_inv, b)
    while(error >= delta and steps <= max_iteration):
        u = u_new
        u_new = np.dot(B, u) + y
        error = np.linalg.norm(b - np.dot(A,u),ord=np.inf)
        steps += 1
        print("GS: " + str(steps))
    return u, steps

# visualisation
def plot(u, n=N, title: str = "title"):
    u = np.array(u)
    U = u.reshape((n+1), (n+1))
    x = np.linspace(0,1,n+1)
    y = np.linspace(0,1,n+1)
    X, Y = np.meshgrid(x, y)
   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, U, cmap='viridis')

    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig(title)
    plt.show()
    return

f_poi = lambda x, y: -4
g_poi = lambda x, y: x**2 + y**2
def create_A_b_Poisson(f, g, epsilon=1, n=N):
    h = 1/n
    A=np.zeros(((n+1)**2,(n+1)**2))
    b=np.zeros((n+1)**2)
    for j in range(0,n+1):
            for i in range(0,n+1):
                  if i==0 or j==0 or i==n or j==n:
                      A[(n+1)*j + i][(n+1)*j + i] =1
                      b[(n+1)*j + i] = g(i*h, j*h)
                  else:
                      A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h) 
                      A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h) 
                      A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h)
                      A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h) 
                      A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h) 
                      b[(n+1)*j + i]                                = f(i*h, j*h)
    return A,b

def Aufgabe_a():
    A,b = create_A_b_Poisson(f_poi, g_poi)
  
    u, step = JC(A, b)
    plot(u, N, "Aufgabe_a_JC")

    u, step = GS(A, b)
    plot(u, N, "Aufgabe_a_GS")

def Aufgabe_b():
    #TODO
    pass

if __name__=="__main__":
    Aufgabe_a()