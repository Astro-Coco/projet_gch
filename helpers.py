import math
def vitesse_developpée(x,y, params):
    vy = 0
    vx = 3/2*params.U*(1 - (2*y/params.H)**2)
    return vx, vy


def vitesse_developpement(x,y, params):
    vy = 0
    delta_x = math.sqrt(params.mu*x/params.rho/params.U)
    vx = params.U*(1 - math.exp(-y/delta_x))
    return vx, vy