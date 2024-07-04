# 3. Programmieruebung 
# Testat am 16,07,2024  15:30-15:45  bei Tamme Claus

import numpy as np
import matplotlib as plt

# Aufgabe 1
# Jacobi 
def JC(A, b, delta=1e-4, max_iteration=1000):
    m = len(b)
    steps = 0
    u = np.zeros(m)
    u_new = np.zeros(m)
    error = 1
    while(error >= delta and steps <= max_iteration):
        for i in range(m):
            temp = 0
            for j in range(m):
                if (j != i):
                    temp += A[i][j] * u[j] 
            u_new[i] = (b[i] - temp) / A[i][i]
        u = u_new
        error = np.linalg.norm(b - np.dot(A,u),ord=np.inf)
        steps += 1
    return u, steps

# Gauss-Seidel
def GS(A,b, delta=1e-4, max_iteration=1000):
    m = len(b)
    steps = 0
    u = np.zeros(m)
    u_new = np.zeros(m)
    error = 1
    while(error >= delta and steps <= max_iteration):
        for i in range(m):
            temp1 = 0
            temp2 = 0
            for j in range(i):
                temp1 += A[i][j] * u_new[j]
            for j in range(i+1,m):
                temp2 += A[i][j] * u[j]
            u_new[i] = (b[i] - temp1 - temp2) / A[i][i]
        u = u_new
        error = np.linalg.norm(b - np.dot(A,u),ord=np.inf)
        steps += 1
    return u, steps

# visualisation
def plot(u, n):
    x = np.linspace(0, 1, n+1)
    y = np.linspace(0, 1, n+1)
    X, Y = np.meshgrid(x, y)
    U = u.reshape((n-1,n-1))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X[1:-1, 1:-1], Y[1:-1, 1:-1], U, cmap='viridis')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('U')
    plt.show()
    return

# test
def testJC():
    A = np.array([[2,1],[5,7]])
    b = np.array([11,13])
    x1,s1 = JC(A,b)
    x2 = np.linalg.solve(A,b)
    print(x1,s1)
    print(x2)
    return

def testGS():
    A = np.array([[16,3],[7,-11]])
    b = np.array([11,13])
    x1,s1 = GS(A,b)
    x2 = np.linalg.solve(A,b)
    print(x1,s1)
    print(x2)
    return


# TODO !!!
# create A and b for Posisson-Problem
def create_A(n):
    # TODO
    return A

def create_b(n):
    # TODO
    return b

# Aufgabe b
# TODO

if __name__=="__main__":
    # do something
    #print("hello world")
    testJC()
    testGS()