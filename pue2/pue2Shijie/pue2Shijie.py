# 2.Programmieruebung #
# Testat am 18,06,2024  11:30-11:45  bei Tamme Claus
# explicit math. Form see skript

import numpy as np
import matplotlib.pyplot as plt

# Parameters
n = 60
h = 1 / n
beta = 5 * np.pi / 6
epsilon_values = [1, 1e-2, 1e-4]
cos = np.cos
sin = np.sin

# Function handles for konvektions-Diffusionsproblem
f_kd = lambda x, y: 1  
g_kd = lambda x, y: 0

def create_A_b(f, g, epsilon, _upwind = False):
    A=np.zeros(((n+1)**2,(n+1)**2))
    b=np.zeros((n+1)**2)
    if _upwind == False:     
        for j in range(0,n+1):
            for i in range(0,n+1):
                  if i==0 or j==0 or i==n or j==n:
                      A[(n+1)*j + i][(n+1)*j + i] =1
                      b[(n+1)*j + i] = g(i*h, j*h)
                  else:
                      A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h) - cos(beta)/(2*h)
                      A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h) + cos(beta)/(2*h)
                      A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h)
                      A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h) - sin(beta)/(2*h)
                      A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h) + sin(beta)/(2*h)
                      b[(n+1)*j + i]                                = f(i*h, j*h)
    if _upwind == True:
         for j in range(0,n+1):
            for i in range(0,n+1):
                  if i==0 or j==0 or i==n or j==n:
                      A[(n+1)*j + i][(n+1)*j + i] =1
                      b[(n+1)*j + i] = g(i*h, j*h)
                  else:
                      cos_h = cos(beta)/h
                      sin_h  = sin(beta)/h
                      if cos(beta)>=0:
                           if sin(beta)>=0:
                              A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h)    - cos_h  
                              A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h) 
                              A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h) +cos_h + sin_h
                              A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h)                   - sin_h
                              A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h) 
                              b[(n+1)*j + i]                                = f(i*h, j*h)
                           else:
                              A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h)    - cos_h
                              A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h)                  
                              A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h) +cos_h  -sin_h
                              A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h) 
                              A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h)                   +sin_h
                              b[(n+1)*j + i]                                = f(i*h, j*h)
                      else:
                             if sin(beta)>=0:
                              A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h) 
                              A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h)    +cos_h  
                              A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h)  -cos_h  +sin_h
                              A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h)                    -sin_h
                              A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h) 
                              b[(n+1)*j + i]                                = f(i*h, j*h)
                             else:
                              A[(n+1)*j + i][(n+1)*j + (i-1)]   = -epsilon/(h*h) 
                              A[(n+1)*j + i][(n+1)*j + (i+1)]  = -epsilon/(h*h)    +cos_h 
                              A[(n+1)*j + i][(n+1)*j + i]          = 4*epsilon/(h*h)  -cos_h   -sin_h
                              A[(n+1)*j + i][(n+1)*(j-1)+ i]    = -epsilon/(h*h) 
                              A[(n+1)*j + i][(n+1)*(j+1)+ i]   = -epsilon/(h*h)                   +sin_h
                              b[(n+1)*j + i]                                = f(i*h, j*h)
    return A,b

def Aufgabe1():
    A,b = create_A_b(f_kd, g_kd, epsilon_values[2], _upwind = True)
    u = np.linalg.solve(A, b)
    U = u.reshape((n+1), (n+1))
    x = np.linspace(0,1,n+1)
    y = np.linspace(0,1,n+1)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, U, cmap='viridis')

    ax.set_title('3D Surface Plot')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig('A1_epsilon=10_-4', dpi=300, bbox_inches='tight')
    plt.show()

def Aufgabe2():
    A,b = create_A_b(f_kd, g_kd, epsilon_values[2], _upwind = False)
    u = np.linalg.solve(A, b)
    U = u.reshape((n+1), (n+1))
    x = np.linspace(0,1,n+1)
    y = np.linspace(0,1,n+1)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, U, cmap='viridis')

    ax.set_title('eps=1 _upwind=False')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig('A2_epsilon=10-4', dpi=300, bbox_inches='tight')
    plt.show()
   
Aufgabe2()