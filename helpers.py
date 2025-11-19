import math
import numpy as np

def vitesse_developpée(x,y, params):
    vy = 0
    vx = 3/2*params.U_in*(1 - (2*y/params.H)**2)
    return vx, vy


def vitesse_developpement(x,y, params):
    vy = 0

    
    def fx() : return (1 - math.exp(-x/params.Ldev))

    vx = params.U_in*(1- fx()) + fx()*vitesse_developpée(x,y, params)[0]
    return vx, vy

def vitesse(x,y, params):
    #Plus petite que la longueur de développement, alors developpemt
    #Si plus grand, alors développée?
    #Si c'est pas la logique, juste mettre developpemt partourt 
    #et l'équation tend ver U_dev
    if x < params.Ldev:
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



def mdf_assemblage(X : tuple, Y : tuple, nx : int, ny : int, params : Params):
    """ Fonction assemblant la matrice A et le vecteur b

    Entrées:
        - X : Bornes du domaine en x, X = [x_min, x_max]
        - Y : Bornes du domaine en y, Y = [y_min, y_max]
        - nx : Discrétisation de l'espace en x (nombre de points)
        - ny : Discrétisation de l'espace en y (nombre de points)
        - params : Structure de données contenant les paramètres physiques du problème

    Sorties (dans l'ordre énuméré ci-bas):
        - A : Matrice (array)
        - b : Vecteur (array)
    """
    

    x_mat, y_mat = position(X, Y, nx, ny)
    

    dx = (X[1] - X[0]) / (nx - 1)
    dy = (Y[1] - Y[0]) / (ny - 1)

    #Simplification des paramètres
    const = params.rho * params.cp /2/params.k/dx

    N = nx * ny
    A = np.zeros((N, N))
    B = np.zeros(N)

    # Dirichlet 
    T_edge = params.T_w
    T_entrée = params.T_in


    for i in range(ny):
        for j in range(nx):
            k = idx(i, j, ny)

            #Détection des bords
            gauche_lim   = (j == 0)
            droite_lim  = (j == nx - 1)
            haut_lim    = (i == 0)         
            bas_lim = (i == ny - 1)

            # Supposer que les conditions T_w appliqué en haut et bas sont prioritaires?
            
            if gauche_lim:
                A[k, k] = 1.0
                B[k] = T_entrée
            elif droite_lim:
                A[k, k] = 3/2/dx
                k_gauche = idx(i, j - 1, ny)
                k_gauche2 = idx(i, j - 2, ny)
                A[k, k_gauche] = -2/dx
                A[k, k_gauche2] = 1/2/dx
                B[k] = 0
            elif haut_lim or bas_lim:
                    A[k, k] = 1.0
                    B[k] = T_edge

            else:
                # ---- Interior / Neumann edges ----
                u = vitesse(x_mat[i, j], y_mat[i, j], params)[0]
                #évaluation de v même si on ne l'utilise pas, vitesse nulle

                #possible utilisation du programme dans le futur avec v non nul
                v = vitesse(x_mat[i, j], y_mat[i, j], params)[1]

                #Évaluation des indices k pour les voisins
                k_gauche = idx(i, j - 1, ny)
                k_droite = idx(i, j + 1, ny)
                k_haut = idx(i - 1, j, ny)
                k_bas = idx(i + 1, j, ny)

                # Assemblage de la matrice A et du vecteur b
                A[k, k] = 2/dx/dx + 2/dy/dy
                A[k, k_gauche] = -const*u - 1/dx/dx
                A[k, k_droite] = const*u - 1/dx/dx
                A[k, k_bas] = -1/dy/dy
                A[k, k_haut] = -1/dy/dy
                B[k] = 0
    return A, B
