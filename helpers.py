import math
import numpy as np
from data_structures import *


def vitesse_developpée(x: float, y: float, params) -> tuple[float, float]:
    """
    Calcule le champ de vitesse développé à partir d'un profil analytique.

    Parameters
    ----------
    x : float
        Position horizontale (non utilisée dans ce modèle, incluse pour compatibilité).
    y : float
        Position verticale où la vitesse est évaluée.
    params :
        Objet contenant au minimum les attributs `U_in` (vitesse de référence)
        et `H` (hauteur caractéristique du canal).

    Returns
    -------
    tuple[float, float]
        (vx, vy) où :
        - vx : composante horizontale de la vitesse au point (x, y)
        - vy : composante verticale de la vitesse (nulle dans ce modèle)
    """
    vy = 0
    vx = 3 / 2 * params.U_in * (1 - (2 * y / params.H) ** 2)
    return vx, vy



def vitesse_developpement(x: float, y: float, params) -> tuple[float, float]:
    """
    Calcule le champ de vitesse dans la zone de développement de la couche limite.

    La vitesse est interpolée entre la vitesse uniforme d'entrée et le profil
    développé en fonction d'une longueur de développement caractéristique.

    Parameters
    ----------
    x : float
        Position horizontale dans le canal.
    y : float
        Position verticale où la vitesse est évaluée.
    params :
        Objet contenant au minimum les attributs :
        - U_in : vitesse de référence à l'entrée
        - Ldev : facteur de longueur de développement (sans dimension)

    Returns
    -------
    tuple[float, float]
        (vx, vy) où :
        - vx : composante horizontale de la vitesse au point (x, y)
        - vy : composante verticale de la vitesse (nulle dans ce modèle)
    """
    vy = 0

    def fx() -> float:
        return 1 - math.exp(-x / params.Ldev)

    vx = params.U_in * (1 - fx()) + fx() * vitesse_developpée(x, y, params)[0]
    return vx, vy



def vitesse(x: float, y: float, params) -> tuple[float, float]:
    """
    Calcule le champ de vitesse au point (x, y).

    Actuellement, cette fonction utilise directement le modèle de
    développement de la couche limite pour obtenir le profil de vitesse.

    Parameters
    ----------
    x : float
        Position horizontale dans le canal.
    y : float
        Position verticale où la vitesse est évaluée.
    params :
        Objet contenant les paramètres nécessaires au calcul de
        `vitesse_developpement`.

    Returns
    -------
    tuple[float, float]
        (vx, vy) où :
        - vx : composante horizontale de la vitesse au point (x, y)
        - vy : composante verticale de la vitesse au point (x, y)
    """
    # Utiliser le profil en développement partout (converge vers le profil développé)
    return vitesse_developpement(x, y, params)



def position(X: tuple, Y: tuple, nx: int, ny: int):
    """Fonction générant deux matrices de discrétisation de l'espace

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


def vecteur_en_matrice(T_vec: np.ndarray, nx: int, ny: int) -> np.ndarray:
    """
    Convertit un vecteur de température en matrice 2D (ny x nx).

    Le mapping utilise la même convention d'indexation que `idx(i, j, ny)` et
    que dans `mdf_assemblage`.

    Parameters
    ----------
    T_vec : np.ndarray
        Vecteur 1D contenant les valeurs (taille nx*ny).
    nx : int
        Nombre de points en x (colonnes).
    ny : int
        Nombre de points en y (lignes).

    Returns
    -------
    np.ndarray
        Matrice 2D de taille (ny, nx) contenant les valeurs de T_vec réarrangées.
    """
    T_mat = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            k = idx(i, j, ny)  # même mapping que dans mdf_assemblage
            T_mat[i, j] = T_vec[k]
    return T_mat




def mdf_assemblage(
    X: tuple[float, float],
    Y: tuple[float, float],
    nx: int,
    ny: int,
    params: Params
) -> tuple[np.ndarray, np.ndarray]:
    """
    Assemble la matrice A et le vecteur B pour le schéma MDF stationnaire.

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
        Structure contenant les paramètres physiques (rho, cp, k, mu, U_in, H, T_w, T_in, Ldev, etc.).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (A, B) où :
        - A : matrice du système linéaire de taille (N, N) avec N = nx * ny
        - B : vecteur colonne du système linéaire de taille (N,)
    """

    x_mat, y_mat = position(X, Y, nx, ny)

    dx = (X[1] - X[0]) / (nx - 1)
    dy = (Y[1] - Y[0]) / (ny - 1)

    # Simplification des paramètres
    const = params.rho * params.cp / 2 / params.k / dx

    N = nx * ny
    # Initialisation des matrices et vecteurs
    A = np.zeros((N, N))
    B = np.zeros(N)

    # Conditions de Dirichlet
    T_edge = params.T_w
    T_entrée = params.T_in

    for i in range(ny):
        for j in range(nx):
            k = idx(i, j, ny)

            # Détection des bords
            gauche_lim = (j == 0)
            droite_lim = (j == nx - 1)
            haut_lim = (i == 0)
            bas_lim = (i == ny - 1)

            # On suppose que les conditions appliquées aux extrémités sont prioritaires
            if gauche_lim:
                A[k, k] = 1.0
                B[k] = T_entrée
            elif droite_lim:
                A[k, k] = 3 / 2 / dx
                k_gauche = idx(i, j - 1, ny)
                k_gauche2 = idx(i, j - 2, ny)
                A[k, k_gauche] = -2 / dx
                A[k, k_gauche2] = 1 / 2 / dx
                B[k] = 0
            elif haut_lim or bas_lim:
                A[k, k] = 1.0
                B[k] = T_edge

            else:
                # ---- Intérieur / Bords Neumann ----
                u = vitesse(x_mat[i, j], y_mat[i, j], params)[0]
                # v évaluée pour compatibilité future (vitesse actuelle supposée nulle)
                v = vitesse(x_mat[i, j], y_mat[i, j], params)[1]

                # Évaluation des indices k pour les voisins
                k_gauche = idx(i, j - 1, ny)
                k_droite = idx(i, j + 1, ny)
                k_haut = idx(i - 1, j, ny)
                k_bas = idx(i + 1, j, ny)

                # Assemblage de la matrice A et du vecteur B
                A[k, k] = 2 / dx / dx + 2 / dy / dy
                A[k, k_gauche] = -const * u - 1 / dx / dx
                A[k, k_droite] = const * u - 1 / dx / dx
                A[k, k_bas] = -1 / dy / dy
                A[k, k_haut] = -1 / dy / dy
                B[k] = 0

    return A, B




# SECTION COUCHES LIMITES
def couches_limites(
    resultats,
    Params
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcule les épaisseurs de couches limites thermique et de vitesse.

    La couche limite est définie comme la distance y à partir de la paroi où :
    - T(y) atteint 99 % de la différence entre la paroi et le centre (thermique)
    - U(y) atteint 99 % de la vitesse au centre (vitesse)

    Parameters
    ----------
    resultats :
        Objet contenant au minimum :
        - T_mat : np.ndarray, champ de température (ny, nx)
        - U_mat : np.ndarray, champ de vitesse axiale (ny, nx)
        - y_mat : np.ndarray, coordonnées verticales (ny, nx)
    Params :
        Objet contenant au minimum :
        - T_w : température de paroi

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        (y_limite_T, y_limite_U, y_limite_T_symm, y_limite_U_symm) où :
        - y_limite_T      : épaisseur de couche limite thermique (y ≥ 0) pour chaque x
        - y_limite_U      : épaisseur de couche limite de vitesse (y ≥ 0) pour chaque x
        - y_limite_T_symm : valeur symétrique de la couche limite thermique (y ≤ 0)
        - y_limite_U_symm : valeur symétrique de la couche limite de vitesse (y ≤ 0)
    """
    T_mat = resultats.T_mat
    U_mat = resultats.U_mat
    y_mat = resultats.y_mat
    T_w = Params.T_w

    ny, nx = T_mat.shape

    y_limite_T = np.zeros(nx)
    y_limite_U = np.zeros(nx)
    # Utilisation de la symétrie pour obtenir la couche limite de l'autre côté
    y_limite_T_symm = np.zeros(nx)
    y_limite_U_symm = np.zeros(nx)

    for j in range(nx):
        T_col = T_mat[:, j]
        U_col = U_mat[:, j]
        y_col = y_mat[:, j]

        # Valeurs au centre du canal
        if ny % 2 == 0: 
            i_centre_haut = ny // 2
            i_centre_bas = ny // 2 + 1
            T_centre = (T_col[i_centre_haut] + T_col[i_centre_bas]) / 2
            U_centre = (U_col[i_centre_haut] + U_col[i_centre_bas]) / 2
        else:
            i_centre = ny // 2
            T_centre = T_col[i_centre]
            U_centre = U_col[i_centre]

        T_limite = T_w + 0.99 * (T_centre - T_w)
        U_limite = 0.99 * U_centre

        # Première position où T atteint le seuil et sa position symétrique
        idx_T = np.where(T_col <= T_limite)[0]
        if len(idx_T) == 0:
            y_limite_T[j] = np.nan
            y_limite_T_symm[j] = ny - 1
        else:
            y_limite_T[j] = y_col[idx_T[0]]
            y_limite_T_symm[j] = y_col[ny - idx_T[0] - 1]

        # Première position où U atteint le seuil et sa position symétrique
        idx_U = np.where(U_col >= U_limite)[0]
        if len(idx_U) == 0:
            y_limite_U[j] = np.nan
            y_limite_U_symm[j] = ny - 1
        else:
            y_limite_U[j] = y_col[idx_U[0]]
            y_limite_U_symm[j] = y_col[ny - idx_U[0] - 1]

    return y_limite_T, y_limite_U, y_limite_T_symm, y_limite_U_symm



def filtre_hysteresis(y: np.ndarray, epsilon: float) -> np.ndarray:
    """
    Applique un filtre d'hystérésis 1D pour lisser les variations brusques.

    Pour chaque point, si la différence avec le point précédent dépasse
    le seuil epsilon, la valeur est remplacée par la valeur précédente.
    Les valeurs NaN sont ignorées (pas de filtrage à ces indices).

    Parameters
    ----------
    y : np.ndarray
        Tableau 1D des valeurs à filtrer.
    epsilon : float
        Seuil maximal de variation accepté entre deux points consécutifs.

    Returns
    -------
    np.ndarray
        Tableau filtré de même taille que y.
    """
    y_f = y.copy()
    for i in range(1, len(y_f)):
        if np.isnan(y_f[i]) or np.isnan(y_f[i - 1]):
            continue
        # Seuil de variation acceptée
        if abs(y_f[i] - y_f[i - 1]) > epsilon:
            y_f[i] = y_f[i - 1]   # stabilisation
    return y_f
