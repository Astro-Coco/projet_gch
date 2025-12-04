import numpy as np
from scipy.linalg import solve
from data_structures import *
from helpers import *


def T_analytique(x, y):
    """
    Solution analytique utilisée pour la vérification de convergence.

    Parameters
    ----------
    x : float 
        Coordonnée(s) en x.
    y : float 
        Coordonnée(s) en y.

    Returns
    -------
    np.ndarray
        Valeur(s) de la température analytique aux points (x, y).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    return np.exp(10 * x) * np.cos(10 * y) + 10


def mdf_assemblage_verification(
    X: tuple[float, float],
    Y: tuple[float, float],
    nx: int,
    ny: int,
    params: Params,
    T_analytique
) -> tuple[np.ndarray, np.ndarray]:
    """
    Assemble A et B pour la résolution de Laplace avec conditions analytique (Dirichlet).
    Fonction utilisée pour la vérification de convergence, basée sur mdf_assemblage définie dans helpers

    Parameters
    ----------
    X : tuple[float, float]
        Bornes du domaine en x, X = (x_min, x_max).
    Y : tuple[float, float]
        Bornes du domaine en y, Y = (y_min, y_max).
    nx : int
        Nombre de points de discrétisation en x.
    ny : int
        Nombre de points de discrétisation en y.
    params : Params
        Paramètres physiques.
    T_analytique : callable
        Fonction T(x, y) donnant la solution analytique sur le bord.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (A, B) système linéaire pour le problème de Laplace 2D.
    """
    # Assemblage pour le retrait de l'advection
    x_mat, y_mat = position(X, Y, nx, ny)
    dx = (X[1] - X[0]) / (nx - 1)
    dy = (Y[1] - Y[0]) / (ny - 1)

    N = nx * ny
    A = np.zeros((N, N))
    B = np.zeros(N)

    for j in range(nx):
        for i in range(ny):
            k = idx(i, j, ny)
            x_coord = x_mat[i, j]
            y_coord = y_mat[i, j]

            gauche = (j == 0)
            droite = (j == nx - 1)
            haut = (i == 0)
            bas = (i == ny - 1)

            if gauche or droite or haut or bas:
                A[k, k] = 1.0
                B[k] = T_analytique(x_coord, y_coord)
            else:
                k_g = idx(i, j - 1, ny)
                k_d = idx(i, j + 1, ny)
                k_h = idx(i - 1, j, ny)
                k_b = idx(i + 1, j, ny)

                A[k, k] = 2 / dx**2 + 2 / dy**2
                A[k, k_g] = -1 / dx**2
                A[k, k_d] = -1 / dx**2
                A[k, k_h] = -1 / dy**2
                A[k, k_b] = -1 / dy**2
                B[k] = 0.0

    return A, B


def laplace_2d(
    Nx: int,
    Ny: int,
    L: float,
    H: float,
    params: Params,
    T_analytique
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Résout le problème de Laplace 2D pour un maillage donné et compare à T_analytique.

    Parameters
    ----------
    Nx : int
        Nombre de points en x.
    Ny : int
        Nombre de points en y.
    L : float
        Longueur du domaine en x.
    H : float
        Hauteur du domaine en y.
    params : Params
        Paramètres physiques.
    T_analytique : callable
        Solution analytique T(x, y) utilisée pour les conditions et la comparaison.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, float, float]
        (T_num, T_an, dx, dy) où :
        - T_num : solution numérique (Ny, Nx)
        - T_an  : solution analytique évaluée sur le maillage (Ny, Nx)
        - dx    : pas en x
        - dy    : pas en y
    """
    X_coords = (0, L)
    Y_coords = (0, H)

    A, B = mdf_assemblage_verification(X_coords, Y_coords, Nx, Ny, params, T_analytique)
    T_vec = solve(A, B)
    T_num = vecteur_en_matrice(T_vec, Nx, Ny)

    x_mat, y_mat = position(X_coords, Y_coords, Nx, Ny)
    T_an = T_analytique(x_mat, y_mat)

    dx = L / (Nx - 1)
    dy = H / (Ny - 1)
    return T_num, T_an, dx, dy


def erreur_L2(T_num: np.ndarray, T_an: np.ndarray) -> float:
    """
    Calcule l'erreur L2 moyenne entre la solution numérique et analytique.

    Parameters
    ----------
    T_num : np.ndarray
        Solution numérique 2D.
    T_an : np.ndarray
        Solution analytique 2D.

    Returns
    -------
    float
        Erreur L2 moyenne.
    """
    err = T_num[1:-1, 1:-1] - T_an[1:-1, 1:-1]
    return np.sqrt(np.mean(err**2))


def ordre(
    erreur_grossier: float,
    erreur_fin: float,
    facteur_raffinement: float
) -> float:
    """
    Estime l'ordre de convergence p à partir de deux erreurs et d'un facteur de raffinement.

    Parameters
    ----------
    erreur_grossier : float
        Erreur sur le maillage grossier.
    erreur_fin : float
        Erreur sur le maillage raffiné.
    facteur_raffinement : float
        Facteur de raffinement spatial (ex. 2 pour Nx_fin = 2 * Nx_grossier).

    Returns
    -------
    float
        Ordre de convergence estimé p. np.nan si une des erreurs vaut 0.
    """
    if erreur_fin == 0 or erreur_grossier == 0:
        return np.nan
    p = np.log(erreur_grossier / erreur_fin) / np.log(facteur_raffinement)
    return p


# Hauteur doit être plus grande pour avoir des erreurs assez grandes
H_val = 1.0
L_val = 1.0

params = Params(
    H=H_val,
    mu=0.001,
    rho=1000,
    T_in=298,
    T_w=373,
    cp=4186,
    k=0.6,
    L=1.0,
    U_in=1.0,
    Ldev=0.05,
    n=1
)

facteur_raffinement = 2
initial_Nx = 40
initial_Ny = 40

grids = [
    (initial_Nx, initial_Ny),
    (initial_Nx * facteur_raffinement, initial_Ny * facteur_raffinement),
    (initial_Nx * facteur_raffinement**2, initial_Ny * facteur_raffinement**2)
]

err_L2 = []

for k_grid, (Nx, Ny) in enumerate(grids):
    T_num, T_an, dx, dy = laplace_2d(Nx, Ny, L_val, H_val, params, T_analytique)
    E = erreur_L2(T_num, T_an)
    err_L2.append(E)

if len(err_L2) >= 2:
    p_32 = ordre(err_L2[0], err_L2[1], facteur_raffinement)
    p_21 = ordre(err_L2[1], err_L2[2], facteur_raffinement)
    print(f"Ordre entre M3 et M2 : p = {p_32:.4f}")
    print(f"Ordre entre M2 et M1 : p = {p_21:.4f}")
    print(f"Ordre observé final : p ≈ {p_21:.4f}")
