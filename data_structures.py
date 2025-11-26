class Params:
    def __init__(self, H, mu, rho, T_in, T_w, cp, k, L, U_in, Ldev, n):
        self.H = H
        self.mu = mu
        self.rho = rho
        self.T_in = T_in
        self.T_w = T_w
        self.cp = cp
        self.k = k
        self.L = L
        self.U_in = U_in

        Re = self.rho*self.U_in*self.H/self.mu
        self.Ldev = Ldev*Re*self.H
        self.n= n


class Results:
    def __init__(self, x, y, T, U,y_limite_T=None,y_limite_U=None):
        self.x_mat = x
        self.y_mat = y
        self.T_mat = T
        self.U_mat = U
        self.y_limite_T=y_limite_T
        self.y_limite_U=y_limite_U

