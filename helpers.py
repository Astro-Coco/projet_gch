import math
import numpy as np
from data_structures import Params, Results

def vitesse_developpée(x,y, params):
    vy = 0
    vx = 3/2*params.U_in*(1 - (2*y/params.H)**2)
    return vx, vy


def vitesse_developpement(x,y, params):
    vy = 0

    Re = params.rho*params.U_in*params.H/params.mu
    Ldev = params.Ldev*Re*params.H

    def fx() : return (1 - math.exp(-x/Ldev))

    vx = params.U_in*(1- fx()) + fx()*vitesse_developpée(x,y, params)[0]
    return vx, vy

def vitesse(x,y, params):
    #Plus petite que la longueur de développement, alors developpemt
    #Si plus grand, alors développée?
    #Si c'est pas la logique, juste mettre developpemt partourt 
    #et l'équation tend ver U_dev

    return vitesse_developpement(x,y, params)

    
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

    matrice_x = np.zeros((ny, nx))
    for i in range(ny):
        matrice_x[i, :] = x_values

    matrice_y = np.zeros((ny, nx))
    for j in range(nx):
        matrice_y[:, j] = y_values

    return matrice_x, matrice_y

def idx(i, j, ny):
    return i + j * ny

def vecteur_en_matrice(T_vec, nx, ny):
    T_mat = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            k = idx(i, j, ny)  # même mapping que dans mdf_assemblage
            T_mat[i, j] = T_vec[k]
    return T_mat
