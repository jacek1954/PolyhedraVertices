import math

epsilon = 0.00000000001
pi = math.acos(-1)
p = 1
q = 1
r = 1


def zero(t):
    x = ang1(t)
    y = ang2(t)
    z = ang3(t)
    v = 2*pi-3*t-x-y-z
    return v


def ang1(vl):
    vl = 1 - 4*(1 - math.cos(vl))*(math.cos(pi/p))**2
    vl = math.acos(vl)
    return vl


def angdiff1(vl):
    nom = 4*((math.cos(pi/p))**2)*math.sin(vl)
    denom = math.sqrt((1 - (1 - 4*(1 - math.cos(vl))*(math.cos(pi/p))**2)**2))
    return nom/denom


def ang2(vl):
    vl = 1 - 4*(1 - math.cos(vl))*(math.cos(pi/q))**2
    vl = math.acos(vl)
    return vl


def angdiff2(vl):
    nom = 4*((math.cos(pi/q))**2)*math.sin(vl)
    denom = math.sqrt((1 - (1 - 4*(1 - math.cos(vl))*(math.cos(pi/q))**2)**2))
    return nom/denom


def ang3(vl):
    vl = 1 - 4*(1 - math.cos(vl))*(math.cos(pi/r))**2
    vl = math.acos(vl)
    return vl


def angdiff3(vl):
    nom = 4*((math.cos(pi/r))**2)*math.sin(vl)
    denom = math.sqrt((1 - (1 - 4*(1 - math.cos(vl))*(math.cos(pi/r))**2)**2))
    return nom/denom


def newton(first):
    nxt = first - (zero(first)/zerodiff(first))
    rest = zero(nxt)
    convergence = 0
    while math.fabs(rest) > epsilon:
        nxt = nxt - (zero(nxt)/zerodiff(nxt))
        rest = zero(nxt)
        convergence = convergence + 1
        if convergence > 100:
            print('too many steps')
            break
    print('number of steps = ', convergence)
    return nxt


def zerodiff(t):
    diff = -3 - angdiff1(t) - angdiff2(t) - angdiff3(t)
    return diff


p = 3
q = 3
r = 5/2
solution = newton(pi/3)
print('snub angle = ', solution)
print('p angle = ', ang1(solution))
print('q angle = ', ang2(solution))
print('r angle = ', ang3(solution))
ro = 1/(2*math.sin(solution/2))
print('base of pyramid radius = ', ro)
R = 1/(2*math.sqrt(1 - ro**2))
print('circumradius = ', R)
