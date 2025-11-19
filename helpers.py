import math
import numpy as np
def vitesse_developpée(x,y, params):
    vy = 0
    vx = 3/2*params.U*(1 - (2*y/params.H)**2)
    return vx, vy


def vitesse_developpement(x,y, params):
    vy = 0
    delta_x = math.sqrt(params.mu*x/params.rho/params.U)
    vx = params.U*(1 - math.exp(-y/delta_x))
    return vx, vy

def vitesse(x,y, params, Ldev = 0.05):
    if x < Ldev:
        return vitesse_developpement(x,y, params)
    else:
        return vitesse_developpée(x,y, params)
    
def position(X : tuple,Y : tuple,nx : int, ny : int):
    """ Fonction générant deux matrices de discrétisation de l'espace

    Entrées:
        - X : Bornes du domaine en x, X = [x_min, x_max]
        - Y : Bornes du domaine en y, Y = [y_min, y_max]
        - nx : Discrétisation de l'espace en x (nombre de points)
        - ny : Discrétisation de l'espace en y (nombre de points)

    Sorties (dans l'ordre énuméré ci-bas):
        - x : Matrice (array) de dimension (ny x nx) qui contient la position en x
        - y : Matrice (array) de dimension (ny x nx) qui contient la position en y
            * Exemple d'une matrice position :
            * Si X = [-1, 1] et Y = [0, 1]
            * Avec nx = 3 et ny = 3
                x = [-1    0    1]
                    [-1    0    1]
                    [-1    0    1]

                y = [1    1    1  ]
                    [0.5  0.5  0.5]
                    [0    0    0  ]
    """


    x_values = np.linspace(X[0], X[1], nx)
    y_values = np.linspace(Y[0], Y[1], ny)  

    matrice = np.zeros((ny, nx))
    for i in range(ny):
        matrice[i, :] = x_values

    matrice_y = np.zeros((ny, nx))
    for j in range(nx):
        matrice_y[:, j] = y_values

    return matrice, matrice_y

def idx(i, j, ny):
    return i + j * ny