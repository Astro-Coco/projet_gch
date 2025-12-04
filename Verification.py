import numpy as np
from scipy.linalg import solve
from data_structures import *
from helpers import *

def T_analytique(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    return np.exp(10 * x) * np.cos(10 * y) + 10

def mdf_assemblage_verification(X: tuple, Y: tuple, nx: int, ny: int, params: Params, T_analytique):
    
    #Assemblage pour le retrait de l'advection
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


def laplace_2d(Nx, Ny, L, H, params, T_analytique):
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


def erreur_L2(T_num, T_an):
    err = T_num[1:-1, 1:-1] - T_an[1:-1, 1:-1]
    return np.sqrt(np.mean(err**2))


def ordre(erreur_grossier, erreur_fin, facteur_raffinement):
    if erreur_fin == 0 or erreur_grossier == 0:
        return np.nan
    p = np.log(erreur_grossier / erreur_fin) / np.log(facteur_raffinement)
    return p

#Hauteur doit être plus grande pour avoir des erreurs assez grandes
H_val = 1.0
L_val = 1.0

params = Params( H = H_val,
    mu = 0.001,
    rho = 1000,
    T_in = 298,
    T_w = 373,
    cp = 4186,
    k = 0.6,
    L = 1.0,
    U_in = 1.0,
    Ldev = 0.05,
    n=1)

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
