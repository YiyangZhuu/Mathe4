# 2.Programmieruebung #
# Testat am 18,06,2024  11:30-11:45  bei Tamme Claus

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# Parameters
n = 60
h = 1 / n
beta = 5 * np.pi / 6
epsilon_values = [1, 1e-2, 1e-4]

# Grid setup
x = np.linspace(0, 1, n + 1)
y = np.linspace(0, 1, n + 1)
X, Y = np.meshgrid(x, y)

# Function handles for boundary conditions
f_cd = 1  
g_cd = 0  

# Parameters for Aufgabe c)
f_poisson = -4  
g_poisson = lambda x, y: x**2 + y**2  

# Function definitions (TODO)
def create_A_b(n, h, f):
    return A, b

def apply_boundary_conditions(A, b, n, g):
    return A, b


def Aufgabe_a():
    return 0

def Aufgabe_b():
    return

def Aufgabe_c():
    return