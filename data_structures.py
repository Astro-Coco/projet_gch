import numpy as np
class Params:
    """
    Contient les paramètres physiques et géométriques du problème.

    Parameters
    ----------
    H : float
        Demi-hauteur ou hauteur caractéristique du canal.
    mu : float
        Viscosité dynamique du fluide.
    rho : float
        Masse volumique du fluide.
    T_in : float
        Température d'entrée du fluide.
    T_w : float
        Température de paroi.
    cp : float
        Capacité calorifique massique à pression constante.
    k : float
        Conductivité thermique.
    L : float
        Longueur du domaine en x.
    U_in : float
        Vitesse d'entrée moyenne (ou de référence).
    Ldev : float
        Facteur sans dimension pour la longueur de développement.
    n : int
        Paramètre numérique additionnel (ex. nombre de points, puissance, etc.).
    """
    def __init__(
        self,
        H: float,
        mu: float,
        rho: float,
        T_in: float,
        T_w: float,
        cp: float,
        k: float,
        L: float,
        U_in: float,
        Ldev: float,
        n: int
    ) -> None:
        self.H = H
        self.mu = mu
        self.rho = rho
        self.T_in = T_in
        self.T_w = T_w
        self.cp = cp
        self.k = k
        self.L = L
        self.U_in = U_in

        Re = self.rho * self.U_in * self.H / self.mu
        self.Ldev = Ldev * Re * self.H
        self.n = n



class Results:
    """
    Regroupe les champs de résultats numériques et les épaisseurs de couches limites.

    Parameters
    ----------
    x : np.ndarray
        Matrice (ny, nx) des positions en x.
    y : np.ndarray
        Matrice (ny, nx) des positions en y.
    T : np.ndarray
        Matrice (ny, nx) des températures.
    U : np.ndarray
        Matrice (ny, nx) des vitesses axiales.
    y_limite_T : np.ndarray, optional
        Épaisseur de couche limite thermique (y ≥ 0) pour chaque x.
    y_limite_U : np.ndarray, optional
        Épaisseur de couche limite de vitesse (y ≥ 0) pour chaque x.
    y_limite_T_symm : np.ndarray, optional
        Symétrique de la couche limite thermique (y ≤ 0).
    y_limite_U_symm : np.ndarray, optional
        Symétrique de la couche limite de vitesse (y ≤ 0).
    """
    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
        T: np.ndarray,
        U: np.ndarray,
        y_limite_T: np.ndarray | None = None,
        y_limite_U: np.ndarray | None = None,
        y_limite_T_symm: np.ndarray | None = None,
        y_limite_U_symm: np.ndarray | None = None
    ) -> None:
        # L'initialisation par défaut à None permet de calculer ces valeurs plus tard
        self.x_mat = x
        self.y_mat = y
        self.T_mat = T
        self.U_mat = U
        self.y_limite_T = y_limite_T
        self.y_limite_U = y_limite_U
        self.y_limite_T_symm = y_limite_T_symm
        self.y_limite_U_symm = y_limite_U_symm

