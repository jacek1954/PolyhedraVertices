import numpy as np
import sbrouts
import movetoplane
import math

def BowersName():
    # name = 'siddid'
    name = 'sided'
    # name = 'qrid'
    # name = 'gaquatid'
    # name = 'girsid'
    # name = 'gisid'
    return name


def FaceId():
    gen = 0
    num = 2
    return [gen, num]


zero = [0, 0, 0]
fi = (math.sqrt(5) + 1)/2  # golden number
epsilon = 0.00001
pi = math.acos(-1)

if BowersName() == 'sided':
    symbol = [[3, 1], [3, 1], [5, 1], [3, 1], [5, 3], [3, 1], 1]  # sided
elif BowersName() == 'gosid':
    pass
elif BowersName() == 'gaquatid':
    pass
elif BowersName() == 'did':
    pass
elif BowersName() == 'seside':
    symbol = [[3, 1], [3, 1], [3, 1], [3, 1], [5, 2], [3, 1], 1]  # seside
elif BowersName() == 'siddid':
    symbol = [[5, 1], [3, 1], [2, 1], [3, 1], [5, 2], [3, 1], 1]  # siddid
elif BowersName() == 'qrid':
    pass
elif BowersName() == 'girsid':
    pass
elif BowersName() == 'gisid':
    symbol = [[3, 1], [3, 1], [2, 1], [3, 1], [5, 3], [3, 1], 1]  # gisid
else:
    faces = []
r = sbrouts.CircumradiusOfSnubPoly(symbol)
rv = math.sqrt(1-(1/(4*(r**2))))  # circumscribed radius of the vertex figure
h = (1/(2*r))  # height of the vertex pyramid
northpole = [0, 0, r]  # default vertex of the pyramid
greenwich = [rv, 0, r - h]  # default first point of the vertex figure
if r > 1:  # difference between latitudes of northpole and greenwich
    stepangle = math.asin(rv / r)
else:
    stepangle = pi - math.asin(rv / r)
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
    for i in range(len(faces[gender])):
        checkf = 0
        for j in range(len(face)):
            for k in range(len(faces[gender][0])):
                if sbrouts.dist(face[j], faces[gender][i][k]) < epsilon:
                    checkf = checkf + 1
                    break
        if checkf == len(face):
            return 1  # face exists
    return 0  # face not exists


def FaceRefill(fragment, symbol):  # completes missing points of the face
    center = CircumCenterOfTriangle(fragment)
    first = fragment[0]
    second = fragment[1]
    last = fragment[2]
    nominator = symbol[0]
    denominator = symbol[1]
    angle = ((2*pi)/nominator)*denominator
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
    base = []
    for k in range(1, len(pyramid)):
        base.append(pyramid[k])
    mod = len(base)
    countf = 0
    for i in range(mod):
        fragm = [peak, base[i], base[(i+1) % mod]]
        gen = FaceGenderNumber(config[i])
        if gen == -1:
            continue
        face = FaceRefill(fragm, config[i])
        if CheckFace(face, gen) == 0:
            faces[gen].append(face)
            countf = countf+1
    return countf


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
    mu = sbrouts.dist(u, zero)**2
    mv = sbrouts.dist(v, zero)**2
    uv = AxisDirection(u, v)
    muv = sbrouts.dist(uv, zero)**2
    for i in range(3):
        d.append(mu*v[i] - mv*u[i])
    e = AxisDirection(d, uv)
    for i in range(3):
        ret.append((e[i]/(2*muv)) + c[i])
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


diags = diagonal(symbol)
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
countnew = 0
allpyramids = [origin]
while pyramids != []:
    temppyramids = []
    for i in range(len(pyramids)):
        for j in range(1, len(pyramids[i])):
            check = CheckPoint(pyramids[i][j])
            if check == 0:
                countnew = countnew + 1
                if countnew > 150:
                    break
                image = NewPyramid(pyramids[i], j)
                temppyramids.append(image)
                allpyramids.append(image)
                vertices.append(pyramids[i][j])
                newfaces = NewFaces(image, symbol)
                facecount = facecount + newfaces
    pyramids = temppyramids
print('number of vertices = ', len(vertices))
print('number of faces  =', facecount)
for i in range(len(faces)):
    print(len(faces[i]))
print('number of pyramids', len(allpyramids))

faces[1] = faces[2]
faces[2] = faces[6]


# tworzenie zbioru punktów do rysowania siatki - wspórzędne wierzchołków ściany i odcinków z FaceMap
rct = sbrouts.FaceMap(FaceId()[0], FaceId()[1], faces, BowersName())
pctcoord = []
for i in range(len(faces[FaceId()[0]][FaceId()[1]])):  # współrzędne wierzchołków ściany
    vrt = faces[FaceId()[0]][FaceId()[1]][i]
    vrt = np.reshape(vrt, 3)
    pctcoord.append(vrt)
for i in range((len(rct))):  # współrzędne odcinków mapy
    pctcoord.append(rct[i][0])
    pctcoord.append(rct[i][1])
cop = sbrouts.CoplanarCheck(FaceId(), faces)
# print('cop = ', cop)
if cop != []:  # współrzędne ściany współpłaszczyznowej (o ile istnieje)
    for i in range(len(cop)):
        pctcoord.append(cop[i])

# przeniesienie pctcoord na płaszzczyznę xy - jako convcoord
convcoord = movetoplane.MoveToXYPlane(pctcoord)

# zaokrąglenie danych do 5 miejsc po przecinku
for i in range(len(convcoord[0])):
    for j in range(3):
        convcoord[0][i][j] = round(convcoord[0][i][j], 5)

# przejście do 2 wymiaru - wyrzucenie niepotrzebnej trzeciej współrzędnej (zerowej)
for i in range(len(convcoord[0])):
    convcoord[0][i] = np.delete(convcoord[0][i], 2, 0)

# tworzenie pliku ze współrzędnymi

templatedata = open('danetest', 'w')
recf = str(len(faces[FaceId()[0]][0]))
templatedata.write(recf)  # zapisujemu na początku pliku liczbę wierzchołków ściany
for i in range(len(convcoord[0])):
    srecf = str(convcoord[0][i])
    templatedata.write(srecf)
templatedata.write('\n')
templatedata.close()

