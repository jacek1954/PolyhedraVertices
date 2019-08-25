import math


def valuelinfun(mat, vec):
    ret = []
    for i in range(len(vec)):
        coord = 0
        for j in range(len(vec)):
            coord = coord + mat[i][j] * vec[j]
        ret.append(coord)
    return ret


def rotaroundaxis(p, ax, sn, cs):
    if ax == 'x':
        m = [[1, 0, 0], [0, cs, -sn], [0, sn, cs]]
    elif ax == 'y':
        m = [[cs, 0, sn], [0, 1, 0], [-sn, 0, cs]]
    elif ax == 'z':
        m = [[cs, -sn, 0], [sn, cs, 0], [0, 0, 1]]
    ret = valuelinfun(m, p)
    return ret


def translate(p, vect):
    q = [p[0]-vect[0], p[1]-vect[1], p[2]-vect[2]]
    return q


def movetoxyplane(pointset):
    # pointset jest układem punktów w przestrzeni (co najmniej trzy punkty)
    # funkcja przeprowadza ten układ na płaszzyznę xy
    # pierwszy punkt układu otrzymuje współrzędne (0, 0)
    # drugi punkt trafia na oś x
    trvect = pointset[0]
    for i in range(1, len(pointset)):
        pointset[i] = translate(pointset[i], trvect)
    pointset[0] = [0, 0, 0]
    v = pointset[1]
    s = v[1]
    c = v[2]
    r = math.sqrt(s**2+c**2)
    sn = s/r
    cs = c/r
    for i in range(len(pointset) - 1):
        pointset[i+1] = rotaroundaxis(pointset[i+1], 'x', sn, cs)
    v = pointset[1]
    s = v[2]
    c = v[0]
    r = math.sqrt(s**2+c**2)
    sn = s/r
    cs = c/r
    for i in range(len(pointset) - 1):
        pointset[i+1] = rotaroundaxis(pointset[i+1], 'y', sn, cs)
    v = pointset[2]
    s = -v[2]
    c = v[1]
    r = math.sqrt(s**2+c**2)
    sn = s/r
    cs = c/r
    for i in range(len(pointset) - 2):
        pointset[i+2] = rotaroundaxis(pointset[i+2], 'x', sn, cs)
    return [pointset]


while 0:
    u = [[1, 2, 3], [3, 6, 0], [-2, 3, 1]]

    p = [2, 4, -3]
    # obrót wokół x tak, aby punkt znalazł się na xz (x' = x, y = 0)
    v = p
    s = v[1]
    c = v[2]
    r = math.sqrt(s ** 2 + c ** 2)
    sn = s / r
    cs = c / r
    p = rotaroundaxis(p, 'x', sn, cs)
    print(p)

    # obrót wokół y tak, aby punkt znalazł się na xy (y' = y, z = 0)
    v = p
    s = v[2]
    c = v[0]
    r = math.sqrt(s ** 2 + c ** 2)
    sn = s / r
    cs = c / r
    p = rotaroundaxis(p, 'y', sn, cs)
    print(p)

    # obrót wokół x tak, aby punkt znalazł się na xy (x' = x, z = 0)
    v = p
    s = -v[2]
    c = v[1]
    r = math.sqrt(s ** 2 + c ** 2)
    sn = s / r
    cs = c / r
    p = rotaroundaxis(p, 'x', sn, cs)
    print(p)

    print(u)
    w = movetoxyplane(u)
    print(w)
