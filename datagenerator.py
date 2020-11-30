import math
import numpy as np
import movetoplane
import sbrouts
import coordinates
import danesided

zero = [0, 0, 0]
fi = (math.sqrt(5) + 1)/2  # golden number
# epsilon = 0.00000001
epsilon = 0.0001
pi = math.acos(-1)
symbol = [[3, 1], [3, 1], [5, 1], [3, 1], [5, 3], [3, 1], 1]  # sided
# symbol = [[3, 1], [3, 1], [3, 1], [3, 1], [5, 2], [3, 1], 1]  # seside
# symbol = [[5, 1], [3, 1], [2, 1], [3, 1], [5, 2], [3, 1], 1]  # siddid
# symbol = [[3, 1], [3, 1], [2, 1], [3, 1], [5, 3], [3, 1], 1]  # gisid
# r_sided = 1.12689791279993934351  # circumscribed radius of sided
# r_seside = 1.4581903307387025510  # circumscribed radius of seside
# r_siddid = 1.2744398820380217900
# r = r_sided
# r = r_seside
# r = r_siddid
r = sbrouts.CircumradiusOfSnubPoly(symbol)
rv = math.sqrt(1-(1/(4*(r**2))))  # circumscribed radius of the vertex figure
print('rv = ', rv)

h = (1/(2*r))  # height of the vertex pyramid
northpole = [0, 0, r]  # default vertex of the pyramid
greenwich = [rv, 0, r - h]  # default first point of the vertex figure
# sided = [1, 1, fi, 1, 1-fi, 1]
# seside = [1, 1, 1, 1, fi-1, 1]
if r > 1:  # difference between latitudes of northpole and greenwich
    stepangle = math.asin(rv / r)
else:
    stepangle = pi - math.asin(rv / r)
print('stepangle = ', stepangle)
triangles = []
squares = []
pentagons = []
hexagons = []
octagons = []
decagons = []
pentagrams = []
octagrams = []
decagrams = []
faces = [triangles, squares, pentagons, hexagons, octagons, decagons, pentagrams, octagrams, decagrams]


def FaceGenderNumber(symbol):  # converts face description n/d to number in faces list
    n = symbol[0]  # nominator
    d = symbol[1]  # denominator
    if n == 2:  # digon
        return -1
    if n == 3:  # triangle
        return 0
    elif n == 4:  # square
        return 1
    elif n == 5:  # pentagon or pentagram
        if d == 1 or d == 4:
            return 2
        else:
            return 6
    elif n == 6:  # hexagon
        return 3
    elif n == 8:  # octagon or octagram
        if d == 1 or d == 7:
            return 4
        else:
            return 7
    elif n == 10:  # decagon or decagram
        if d == 1 or d == 9:
            return 5
        else:
            return 8


def CheckFace(face, gender):  # checks if the face has already been found
    # face - given face, gender - face gender
    # print('CheckFace')
    # print('face = ', face)
    # print('gender = ', gender)
    for i in range(len(faces[gender])):
        checkf = 0
        for j in range(len(face)):
            # print('lenfacesgen  = ', len(faces[gender]))
            for k in range(len(faces[gender][0])):
                if sbrouts.dist(face[j], faces[gender][i][k]) < epsilon:
                    checkf = checkf + 1
                    break
        if checkf == len(face):
            # print('face exists')
            return 1  # face exists
    # print('face not exists')
    return 0  # face not exists


def FaceRefill(fragment, symbol):  # completes missing points of the face
    center = CircumCenterOfTriangle(fragment)
    first = fragment[0]
    second = fragment[1]
    last = fragment[2]
    nominator = symbol[0]
    denominator = symbol[1]
    angle = ((2*pi)/nominator)*denominator
    # print('angle =', (angle*180)/pi)
    current = second
    ret = [first, second]
    for i in range(nominator - 3):
        following = RotAroundAxis(center, current, angle)
        ret.append(following)
        current = following
    ret.append(last)
    return ret


def NewFaces(pyramid, config):  # adds new faces designated by a given pyramid
    peak = pyramid[0]
    # del pyramid[0]
    # base = pyramid
    base = []
    for k in range(1, len(pyramid)):
        base.append(pyramid[k])
    mod = len(base)
    countf = 0
    for i in range(mod):
        fragm = [peak, base[i], base[(i+1) % mod]]
        gen = FaceGenderNumber(config[i])
        if gen == -1:
            # print('i = ', i)
            continue
        face = FaceRefill(fragm, config[i])
        if CheckFace(face, gen) == 0:
            faces[gen].append(face)
            countf = countf+1
            # if gen == 2:
            #     print('faceslen = ', len(faces[gen]))
            #     for m in range(len(faces[gen])):
            #         print(faces[gen][m])
    return countf


# def CentralAngle(radius, p, q):  # compute angular distance between points
#     d = sbrouts.dist(p, q)
#     alfa = 2*math.asin(d/(2*radius))
#     ret = []
#     ret.append(math.sin(alfa))
#     ret.append(math.cos(alfa))
#     print('ret= ', ret)
#     return ret


def AxisDirection(p, q):  # vector product
    a = p[1]*q[2]-p[2]*q[1]
    b = -(p[0]*q[2]-p[2]*q[0])
    c = p[0]*q[1]-p[1]*q[0]
    v = [a, b, c]
    return v


def RotAroundAxis(axis, point, angle):  # rotation around any axis
    module = sbrouts.dist(axis, [0, 0, 0])
    x = axis[0]/module
    y = axis[1]/module
    z = axis[2]/module
    s = math.sin(angle)
    c = math.cos(angle)
    cc = 1-c
    w1 = [(x*x*cc)+c, (x*y*cc)-(z*s), (x*z*cc)+(y*s)]
    w2 = [(x*y*cc)+(z*s), c+(y*y*cc), (y*z*cc)-(x*s)]
    w3 = [(z*x*cc)-(y*s), (z*y*cc)+(x*s), c+(z*z*cc)]
    matrix = [w1, w2, w3]
    image = movetoplane.ValueLinFun(matrix, point)
    return image


def CircumCenterOfTriangle(triangle):  # center of circumcircle of the triangle
    a = triangle[1]
    b = triangle[2]
    c = triangle[0]
    u = []
    v = []
    ret = []
    d = []
    for i in range(3):
        u.append(a[i] - c[i])
        v.append(b[i] - c[i])
    # print(u, v, zero)
    mu = sbrouts.dist(u, zero)**2
    mv = sbrouts.dist(v, zero)**2
    uv = AxisDirection(u, v)
    # print(uv)
    muv = sbrouts.dist(uv, zero)**2
    # print(mu, mv, muv)
    for i in range(3):
        d.append(mu*v[i] - mv*u[i])
    # print(d)
    e = AxisDirection(d, uv)
    for i in range(3):
        ret.append((e[i]/(2*muv)) + c[i])
    return ret


def CircumCenterOfTriangle2(triangle):
    t = [0, 0, 0, 0]
    unit = [[1], [1], [1]]
    jnt = np.concatenate((unit, triangle), axis=1)
    for i in range(4):
        s = np.delete(jnt, i, 1)
        t[i] = sbrouts.determinant(s, 3)
    t[2] = -t[2]
    t = np.reshape(t, 4)
    # print('t = ', t)
    # t contains the coefficients of the triangle plane equation
    a = triangle[0]
    b = triangle[1]
    c = triangle[2]
    u = [0, 0, 0, 0]
    for i in range(1, 4):
        u[i] = (b[i-1] - a[i-1])/2
    u[0] = (b[0]**2-a[0]**2 + b[1]**2-a[1]**2 + b[2]**2-a[2]**2)/4
    # print('u = ', u)
    v = [0, 0, 0, 0]
    for i in range(1, 4):
        v[i] = (c[i-1] - a[i-1])/2
    v[0] = (c[0]**2-a[0]**2 + c[1]**2-a[1]**2 + c[2]**2-a[2]**2)/4
    # print('v = ', v)
    # u and v contain the coefficients of equations ad.ds=0 and ae.es=0, where d and e are midpoints of ab and ac
    # and s is the circumcenter of the triangle
    m = [u, v, t]  # m is the matrix of equation system
    n = np.reshape(m, (3, 4))
    # print('n = ', n)
    rside = np.delete(n, (1, 2, 3), 1)
    lside = np.delete(n, 0, 1)
    # print('lside = ', lside)
    # print('rside = ', rside)
    ret = np.zeros((3, 1))
    ret = np.linalg.solve(lside, rside)
    ret = np.reshape(ret, 3)
    return ret


def diagonal(config):  # creates a list of the sides of the pyramid's base
    # (p.3.q.3.r.3) is the vertex configuration of snub polyhedron
    p = config[0][0]/config[0][1]
    q = config[2][0]/config[2][1]
    r = config[4][0]/config[4][1]
    conf = []
    conf.append(2*math.cos(pi/p))
    conf.append(1)
    conf.append(2*math.cos(pi/q))
    conf.append(1)
    conf.append(2*math.cos(pi/r))
    conf.append(1)
    return conf


def NewPyramid(oldpyramid, direction):  # creates pyramids for any point of given pyramid
    oldpeak = oldpyramid[0]
    newpeak = oldpyramid[direction]
    axis1 = AxisDirection(oldpeak, newpeak)
    temppyramid = [newpeak]
    for i in range(1, len(oldpyramid)):
        image = RotAroundAxis(axis1, oldpyramid[i], stepangle)
        temppyramid.append(image)
    # print('temppyr = ', temppyramid)
    newpyramid = [newpeak]
    for j in range(1, len(temppyramid)):
        image = RotAroundAxis(newpeak, temppyramid[j], angles2[direction])
        newpyramid.append(image)
    return newpyramid


def CheckPoint(point):
    for i in range(len(vertices)):
        if math.fabs(sbrouts.dist(point, vertices[i])) <= epsilon:
            return 1
    return 0


def CheckDistance(pyramid, config):
    checklist = []
    for i in range(1, len(pyramid)-1):
        a = sbrouts.dist(pyramid[i], pyramid[i+1])
        b = math.fabs(config[i-1])
        if math.fabs(a-b) <= epsilon:
            checklist.append(1)
        else:
            checklist.append(0)
            # print('a = ', a, 'b = ', b)
    a = sbrouts.dist(pyramid[len(pyramid)-1], pyramid[1])
    b = math.fabs(config[len(config)-1])
    if math.fabs(a - b) <= sbrouts.epsilon:
        checklist.append(1)
    else:
        checklist.append(0)
    return checklist

diags = diagonal(symbol)
# print('diags = ', diags)
angles = []  # successive center angles between the points on the base of pyramid
for i in range(6):
    angles.append(2 * math.asin(diags[i] / (2 * rv)))
angles2 = [0]  # calculates the angles of rotation setting the pyramid to a given meridian
angles2.append(pi - angles[0])
angles2.append(pi + angles[0])
angles2.append(pi - angles[2])
angles2.append(pi + angles[2])
angles2.append(pi - angles[4])
angles2.append(pi + angles[4])
print('angles = ', angles)
print('angles2  =', angles2)
# northpole = danesided.V[0]
# greenwich = danesided.V[48]
origin = []
origin.append(northpole)
origin.append(greenwich)
point = greenwich
for i in range(len(angles) - 1):
    image = RotAroundAxis(northpole, point, angles[i])
    origin.append(image)
    point = image

vertices = [origin[0]]
pyramids = [origin]
facecount = NewFaces(origin, symbol)
# print('origin = ', origin)
# print('check distance', CheckDistance(origin, diags))
# print('last = ', sbrouts.dist(origin[5], origin[1]))
countnew = 0
allpyramids = [origin]
while pyramids != []:
    temppyramids = []
    for i in range(len(pyramids)):
        for j in range(1, len(pyramids[i])):
            check = CheckPoint(pyramids[i][j])
            # print('check = ', check)
            if check == 0:
                countnew = countnew + 1
                if countnew > 150:
                    break
                image = NewPyramid(pyramids[i], j)
                # print('image = ', image)
                # print('check distance', CheckDistance(image, diags))
                temppyramids.append(image)
                allpyramids.append(image)
                vertices.append(pyramids[i][j])
                newfaces = NewFaces(image, symbol)
                # print('newfaces = ', newfaces)
                facecount = facecount + newfaces
    pyramids = temppyramids
print('number of vertices = ', len(vertices))
print('number of faces  =', facecount)
for i in range(len(faces)):
    print(len(faces[i]))
print('number of pyramids', len(allpyramids))
