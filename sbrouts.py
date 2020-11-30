import numpy as np
import math
import althead

epsilon = 0.0000001
fi = (math.sqrt(5) + 1)/2
pi = math.acos(-1)


def determinant(matrix, dimension):
    # oblicza wyznacznik macierzy
    if dimension == 1:
        return matrix[0]
    without1row = np.delete(matrix, 0, 0)
    ret = 0
    for i in range(dimension):
        minormac = np.delete(without1row, i, 1)
        ret = ret + ((-1)**i)*matrix[0][i]*determinant(minormac, dimension - 1)
    return ret


def IntersectionPointOfPlanes(pst):
    # oblicza wspólny punkt trzech płaszczyzn
    # elementy pst są współrzędnymi punktów leżących na trzech danych płaszczyznach
    un = [[1], [1], [1]]
    p = np.delete(pst, (2, 1), 0)
    r = np.reshape(p, (3, 3))
    msk = [[1, 2], [0, 2], [0, 1]]  # msk pozwala wybrać podzbiory 2 elementowe ze zbioru 3 elementowego
    p = np.delete(pst, msk[0], 0)
    # tablica t zawiera 4 współczynniki trzech równań danych płaszczyzn
    t = np.zeros((3, 4))
    for i in range(3):
        p = np.delete(pst, msk[i], 0)
        r = np.reshape(p, (3, 3))
        jnt = np.concatenate((un, r), axis=1)
        for j in range(4):
            s = np.delete(jnt, j, 1)
            t[i][j] = -determinant(s, 3)
    for i in range(3):
        t[i][2] = -t[i][2]
    # t jest obliczone, tworzymy układ równań i rozwiązujemy
    RightSide = np.delete(t, (1, 2, 3), 1)  # RightSide jest wektorem po prawej stronie układu równań
    LeftSide = np.delete(t, 0, 1)  # LeftSide jest macierzą lewej strony ukłądu równań
    w = determinant(LeftSide, 3)
    if math.fabs(w) < epsilon:
        raise np.linalg.LinAlgError
    solution = np.zeros((3, 1))
    solution = np.linalg.solve(LeftSide, RightSide)  # linalg rozwiązuje układ równań
    r = np.reshape(solution, 3)
    return r


def dist(p, q):
    # oblicza odległoość punktów p i q
    # print('p = ', p, 'q = ', q)
    w = 0
    for i in range(3):
        w = w + (p[i] - q[i])**2
    return math.sqrt(w)


def EqualityOfPoints(p, q):
    s = dist(p, q)
    if s <= epsilon:
        return 1
    else:
        return 0


def SumOfTwoVectors(v1, v2):
    for i in range(3):
        ret = [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]
    return ret


def MultiplcationOfVectorAndScalar(v, s):
    ret = [s*v[0], s*v[1], s*v[2]]
    return ret


def DotProduct(p, q):
    ret = 0
    for i in range(3):
        ret = ret + p[i]*q[i]
    return ret


def CrossProduct(p, q):  # vector product
    a = p[1]*q[2]-p[2]*q[1]
    b = -(p[0]*q[2]-p[2]*q[0])
    c = p[0]*q[1]-p[1]*q[0]
    v = [a, b, c]
    return v


def CoefficientsOfTheEquationOfThePlane(plane):
    # calculates a, b, c, d in equation ax+by+cz+d=0
    # plane - list of three points set out the plane
    a1 = plane[0][1]*(plane[1][2] - plane[2][2])
    a2 = plane[1][1]*(plane[2][2] - plane[0][2])
    a3 = plane[2][1]*(plane[0][2] - plane[1][2])
    a = a1 + a2 + a3
    b1 = plane[0][2]*(plane[1][0] - plane[2][0])
    b2 = plane[1][2]*(plane[2][0] - plane[0][0])
    b3 = plane[2][2]*(plane[0][0] - plane[1][0])
    b = b1 + b2 + b3
    c1 = plane[0][0]*(plane[1][1] - plane[2][1])
    c2 = plane[1][0]*(plane[2][1] - plane[0][1])
    c3 = plane[2][0]*(plane[0][1] - plane[1][1])
    c = c1 + c2 + c3
    d1 = plane[0][0]*(plane[1][1]*plane[2][2]-plane[2][1]*plane[1][2])
    d2 = plane[1][0]*(plane[2][1]*plane[0][2]-plane[0][1]*plane[2][2])
    d3 = plane[2][0]*(plane[0][1]*plane[1][2]-plane[1][1]*plane[0][2])
    d = d1 + d2 + d3
    d = -d
    return [a, b, c, d]


def IntersetionOfThreePlanes(plane, qlane, rlane):
    x = CoefficientsOfTheEquationOfThePlane(plane)
    y = CoefficientsOfTheEquationOfThePlane(qlane)
    z = CoefficientsOfTheEquationOfThePlane(rlane)
    n1 = [x[0], x[1], x[2]]
    d1 = x[3]
    n2 = [y[0], y[1], y[2]]
    d2 = y[3]
    n3 = [z[0], z[1], z[2]]
    d3 = z[3]
    denominator = DotProduct(n1, CrossProduct(n2, n3))
    if math.fabs(denominator) <= epsilon:
        return []
    e1 = MultiplcationOfVectorAndScalar(CrossProduct(n2, n3), d1)
    e2 = MultiplcationOfVectorAndScalar(CrossProduct(n1, n3), d2)
    e3 = MultiplcationOfVectorAndScalar(CrossProduct(n2, n1), d3)
    e1 = MultiplcationOfVectorAndScalar(e1, 1 / denominator)
    e2 = MultiplcationOfVectorAndScalar(e2, 1 / denominator)
    e3 = MultiplcationOfVectorAndScalar(e3, 1 / denominator)
    ret = SumOfTwoVectors(e1, e2)
    ret = SumOfTwoVectors(ret, e3)
    return ret


def WhichSideOfPlane(plane, point):
    # plane - list of three points set out the plane
    # point - the point whose position we are examining in relation to the plane
    coeff = CoefficientsOfTheEquationOfThePlane(plane)
    s = coeff[0] * point[0] + coeff[1] * point[1] + coeff[2] * point[2] + coeff[3]
    if math.fabs(s) <= epsilon:
        return 0  # point lies on the plane
    elif s > 0:
        return 1  # point lies on the same side as the point (coeff[0], coeff[1], coeff[2])
    else:
        return -1  # point lies on the opposite side to the point (coeff[0], coeff[1], coeff[2])


def VertexMaps(pyramids, plane):
    zeroside = WhichSideOfPlane(plane, [0, 0, 0])
    northside = []
    southside = []
    maps = []
    mod = len(pyramids[0])  # mod is the number of vertices of pyramid(with peak)
    for i in range(len(pyramids)):
        pointside = WhichSideOfPlane(plane, pyramids[i][0])
        if pointside == zeroside:
            southside.append(pyramids[i])
        else:
            northside.append(pyramids[i])
    for i in range(len(northside)):
        for j in range(mod):
            if WhichSideOfPlane(plane, northside[i][j]) == zeroside:
                if WhichSideOfPlane(plane, northside[i][(j+1) % (mod-1)]) == zeroside:
                    prior = [northside[i][0], northside[i][j], northside[i][(j+1) % (mod-1)]]
                    next = [northside[i][0], northside[i][(j+1) % (mod-1)], northside[i][(j+2) % (mod-1)]]
                    mappoint = IntersetionOfThreePlanes(plane, prior, next)
                    maps.append(mappoint)


# def SearchEdge(facein, numin, facescan):
#     # SearchEdge znajduje ściany (płaszczyzny) mające wspólną krawędź z daną
#     # facein i numin określają rodzaj i numer danej ściany
#     # facescan określa rodzaj poszukiwanej ściany, tzn. jest wybraną listą ścian
#     # wartościa jest lista zawierająca na przemian 2 wierzchołki i numer ściany:
#     # [[[a1 ,b1, c1], [a2, b2, c2]], n1, [[a3, b3, c3], [a4, b4, c4], n2,...]
#     ret = []
#     for i in range(int(len(facescan))):
#         lc = 0
#         vert = []
#         for j in range(len(facein[0])):
#             for k in range(len(facescan[0])):
#                 lcc = 0
#                 crd = []
#                 for l in range(3):
#                     if facein[numin][j][l] == facescan[i][k][l]:
#                         lcc = lcc + 1
#                         crd.append(facescan[i][k][l])
#                 if lcc == 3:
#                     lc = lc + 1
#                     vert.append(crd)
#         if lc == 2:
#             ret.append(vert)
#             ret.append(i)
#     return ret


def SearchEdge(facein, numin, facescan):
    ret = []
    for i in range(int(len(facescan))):
        lc = 0
        vert = []
        for j in range(len(facein[0])):
            for k in range(len(facescan[0])):
                if EqualityOfPoints(facein[numin][j], facescan[i][k]):
                # if dist(facein[numin][j], facescan[i][k]) < epsilon:
                    lc = lc + 1
                    vert.append(facescan[i][k])
                    break
        if lc == 2:
            if math.fabs(dist(vert[0], vert[1]) - 1) < epsilon:
                ret.append(vert)
                ret.append(i)
    return ret


def SegmIntersect(blue, red):
    # print('blue = ', blue, 'red = ', red)
    lengblue = dist(blue[0], blue[1])
    lengred = dist(red[0], red[1])
    dsttab = []
    dsttab.append(lengblue)
    dsttab.append(lengred)
    dsttab.append(dist(blue[0], red[0]))
    dsttab.append(dist(blue[1], red[1]))
    dsttab.append(dist(blue[0], red[1]))
    dsttab.append(dist(blue[1], red[0]))
    indmax = np.argmax(dsttab)
    if indmax == 0:
        return [red[0], red[1]]
    elif indmax == 1:
        return [blue[0], blue[1]]
    else:
        if lengblue + lengred - dsttab[indmax] > epsilon:
            if indmax == 2:
                return [blue[1], red[1]]
            elif indmax == 3:
                return [blue[0], red[0]]
            elif indmax == 4:
                return [blue[1], red[0]]
            else:  # indmax == 5
                return [blue[0], red[1]]
        else:
            return []


def SearchCut(facein, numin, fsctype, edges, faces):
    # searchcut znajduje ściany (płaszczyzny) przecinające daną
    # facein i numin określają rodzaj i numer danej ściany
    # fsctype określa rodzaj poszukiwanej ściany
    # edges zawiera ściany krawędziowe (ale tylko trzech typów) : [etriangle, epentagon, epentagram]
    # każdy element tej listy jest postaci [[[a1 ,b1, c1], [a2, b2, c2]], n1, [[a3, b3, c3], [a4, b4, c4], n2,...]
    # czyli zawiera krawędź, i numer ściany (danego typu)
    ret = []
    # planein, planeedge i planescan są płaszczyznami, których przecięcie szukamy
    planein = np.zeros((3, 3))
    planeedge = np.zeros((3, 3))
    planescan = np.zeros((3, 3))
    for i in range(3):  # określamy planein przez trzy punkty facein
        for j in range(3):
            planein[i][j] = facein[numin][i][j]
    le = int(len(edges[fsctype])/2)
    for nted in range(len(edges)):  # nted wybiera rodzaj krawędzi
        for ned in range(int(len(edges[nted])/2)):  # ned wybiera kolejna krawędź danego rodzaju
            k = 2 * ned + 1  # wybieramy kolejny numer krawędzi
            k = edges[nted][k]  # k staje się numerem krawędzi przed chwilą odnalexionej
            for i in range(3):  # określamy planeedge przez trzy punkty kolejnej ściany krawędziowej
                for j in range(3):
                    planeedge[i][j] = faces[nted][k][i][j]
            for nsc in range(len(faces[fsctype])):  # nsc wybiera kolejną ścianę potencjalnie poszukiwaną
                et = 0                                  #
                for ne in range(le):                    #
                    if nsc == edges[fsctype][2*ne+1]:   #
                        et = 1                          # odrzucenie ścian krawędziowych
                        break                           #
                if et == 1:                             #
                    continue                            #
                try:
                    for i in range(3):  # określamy planescan przez trzy punkty ściany potencjalnie przecinającej
                        for j in range(3):
                            planescan[i][j] = faces[fsctype][nsc][i][j]
                    p = IntersectionPointOfPlanes([planein, planeedge, planescan])  # obliczenie punktu wspólnego
                    q1 = dist(p, edges[nted][2*ned][0])
                    q2 = dist(p, edges[nted][2*ned][1])
                    if math.fabs(q1 + q2 - 1) <= epsilon:         # zbadanie czy punkt leży na krawędzi
                        ret.append(p)
                        ret.append(nted)
                        ret.append(k)
                        ret.append(nsc)
                        if q1 <= q2:
                            ret.append(q1)
                        else:
                            ret.append(q2)
                except np.linalg.LinAlgError:  # jeśli nie ma (jedynego) punktu przeciecia pomiń tę ścianę
                    continue
    # ret zawiera współrzedne punktu przecięcia, rodzaj krawędzi oraz numery krawędzi i ściany przecinającej
    #     a także nie wiekszą z odległości punktu od końca krawędzi
    # [[a1, b1, c1], et1, en1, ct1, d1, [a2, b2, c2], et2, en2, ct2, d2,...]
    return ret


def EdgeBuild(intype, numface, faces):
    # tworzenie listy ścian krawędziowych dla danej ściany
    edg = []
    for i in range(len(faces)):
        edg.append(SearchEdge(faces[intype], numface, faces[i]))
    return edg


def FaceMap(fintype, finnum, faces, PolyhedronName):
    # FaceMap tworzy listę odcinków, które są częściami wspólnymi danej ściany ze ścianami ją przecinającymi
    # fintype określa rodzaj danej ściany, a finnum jej numer
    # wartością funkcji jest lista, której każdy element skłąda się z czterech danych:
    # [typ ściany przecinającej, numer tej ściany, współrzędne odcinka wspólnego, długość tego odcinka]
    # ret = [[fct, fcn, seg, lngseg], [ , , , ], [ , , , ], ...]
    # print('facemap')
    # print(len(faces[0]), len(faces[1]), len(faces[2]))
    ret = []
    ev1 = EdgeBuild(fintype, finnum, faces)  # krawędzie danej ściany
    # print(ev1)
    z = []  # lista ścian (płaszczyzn) przecinających daną
    for i in range(len(faces)):  # są tylko trzy rodzaje ścian
        z.append(SearchCut(faces[fintype], finnum, i, ev1, faces))
    for n in range(len(z)):
        skip = []
        for m in range(int(len(z[n])/5)):
            nmf = z[n][5*m+3]
            try:
                indx = skip.index(nmf)
                # print('indx', indx)
                continue
            except ValueError:
                pass
            ev2 = EdgeBuild(n, nmf, faces)
            y = []
            for l in range(len(faces)):
                y.append(SearchCut(faces[n], nmf, l, ev2, faces))
            # teraz trzeba sprawdzić czy y zawiera finnum dla fintype
            # a jesli tak to utworzyć listy blue i red
            blue = []
            for j1 in range(int(len(y[fintype])/5)):
                if y[fintype][5*j1+3] == finnum:
                    blue.append(y[fintype][5*j1])
                    blue.append(y[fintype][5*j1+4])
            if blue == []:
                skip.append(nmf)
                continue
            # czyli blue jest niepusty, tworzymy red, który pusty nie jest
            red = []
            for j2 in range(m, int(len(z[n])/5)):
                if z[n][5*j2+3] == nmf:
                    red.append(z[n][5*j2])
                    red.append(z[n][5*j2+4])
            skip.append(nmf)
            red = DuplicatesReduct(red)
            blue = DuplicatesReduct(blue)
            if len(red) <= 2:
                continue
            if len(blue) <= 2:
                continue
            # teraz zamieniamy odległości na typy
            genf = PolygonType(fintype, PolyhedronName)
            genc = PolygonType(n, PolyhedronName)
            for j3 in range(int(len(red)/2)):
                red[2*j3+1] = PointType(genf, red[2*j3+1])
            for j4 in range(int(len(blue)/2)):
                blue[2*j4+1] = PointType(genc, blue[2*j4+1])
            # teraz trzeba uporządkowac punkty, ale tylko dla gwiazd
            if genf >= 7:
                red = PointSort(red)
                # print('redtostella', red)
                red = StellaTrace(red)
                # print('stellared', red)
            else:
                gbg = red.pop(1)
                gbg = red.pop(2)
                red = [red]
            if genc >= 7:
                # print(genc)
                # print(blue)
                blue = PointSort(blue)
                blue = StellaTrace(blue)
                # print('stellablue', blue)
            else:
                gbg = blue.pop(1)
                gbg = blue.pop(2)
                blue = [blue]
            for i1 in range(len(blue)):
                for j5 in range(len(red)):
                    rt = SegmIntersect(blue[i1], red[j5])
                    if len(rt) == 0:
                        continue
                    else:
                        ret.append(rt)
    return ret


def DuplicatesReduct(list):
    newlist = []
    while len(list) > 0:
        i = 2
        while i < len(list):
            if dist(list[0], list[i]) < epsilon:
                del list[i]
                del list[i]
            else:
                i = i+2
        newlist.append(list[0])
        newlist.append(list[1])
        del list[0]
        del list[0]
    return newlist


def PointType(fgen, dist):  # ustala typ punktu w zależności od jego odległości od wierzchołka krawędzi i typu ściany
    if dist < epsilon:
        return 'v'  # vertice
    else:  # dist > 0
        if fgen < 7:  # convex face
            return 'e'  # external
        else:  # concave face
            if fgen == 7:  # pentagram
                if math.fabs(dist - (2 - fi)) < epsilon:
                    return 'b'  # border
                elif dist - (2 - fi) < 0:
                    return 'e'
                else:
                    return 'i'  # internal
            elif fgen == 8:  # octagram
                if math.fabs(dist - (1 - (1/math.sqrt(2)))) < epsilon:
                    return 'b'  # border
                elif dist - (1 - (1/math.sqrt(2))) < 0:
                    return 'e'
                else:
                    return 'i'  # internal
            else:  # decagram
                if math.fabs(dist - (2 * fi - 3)) < epsilon:
                    return 'b'  # border
                elif dist - (2 * fi - 3) < 0:
                    return 'e'
                else:
                    return 'i'  # internal


def PolygonType(t, name):  # funkcja tymczasowa dla gosid, sided, gaquatid, did, seside
    # print('t = ', t)
    if name == 'gosid':
        if t == 0:
            return 1
        else:  # t == 1:
            return 7
    elif name == 'sided':
        if t == 0:
            return 1
        elif t == 1:
            return 3
        else:  # t = 2:
            return 7
    elif name == 'gaquatid':
        if t == 0:
            return 2
        elif t == 1:
            return 4
        else:  # t = 3
            return 9
    elif name == 'did':
        if t == 0:
            return 3
        else:  # t = 1
            return 7
    elif name == 'seside':
        if t == 0:
            return 1
        else:  # t = 1
            return 7
    elif name == 'siddid':
        if t == 0:
            return 1
        elif t == 1:
            return 3
        else:
            return 7
    elif name == 'qrid':
        if t == 0:
            return 1
        elif t == 1:
            return 2
        else:
            return 7
    elif name == 'girsid':
        if t == 0:
            return 1
        else:  # t = 1
            return 7

# funkcja wyznacza zbiór odcinków powstałych przez przecięcie ściany gwiaździstej przez płaszczyznę określone punktami
# przecięcia z krawędziami gwiazdy podanymi w tablicy point, code zawiera kody typów punktów z point
# lista punktów zawartych w point jest uporządkowana w przestrzeni
def StellaTrace(point):
    # print('ST')
    # print('point = ', point)
    l = len(point)-2  # pozycja ostatniego punktu
    # print('l = ', l)
    zero = []  # odcinek zerowy
    s = []
    if point[1] == 'b':  # pierwszy punkt jest graniczny - ślad jest jednym odcinkiem
        return [[point[0], point[l]]]
    if point[1] == 'v' and point[3] == 'v':  # ślad jest dwupunktowy
        return zero
    if point[1] == 'e':  # pierwszy punkt jest zewnętrzny
        k = 1
        # print('k = ', k)
    elif point[1] == 'v':  # pierwszy jest wierzchołkiem ale drugi jest zewnętrzny
        k = 3
#  k jest po prostu indeksem punktu zewnętrznego ('e')
    i = k + 2
    # print('i = ', i)
    while i <= l+1:
        point[k] = 'n'  # 'deaktywacja' początku ('e')
        try:  # znajdź koniec odcinka
            i = point.index('e')
            s.append([point[k-1], point[i-1]])
            # print('s = ', s)
            point[i] = 'n'  # 'deaktywacja' końca ('e')
        except ValueError:
            s.append([point[k-1], point[l]])  # nie ma już końca typu 'e' - to ostatni odcinek
            return s
        try:  # znajdź początek odcinka
            k = point.index('e')
        except ValueError:
            return s  # nie ma więcej odcinków


def PointSort(points):
    dim = int(len(points)/2)
    disttab = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i+1, dim):
            disttab[i][j] = dist(points[2*i], points[2*j])
    indmax = np.argmax(disttab)
    ind = np.unravel_index(indmax, disttab.shape)
    ret = []
    for i in range(dim):
        ret.append([0, 0, 0])
        ret.append(0)
    ret[0] = points[2*ind[0]]
    ret[1] = points[2*ind[0]+1]
    ret[2*(dim-1)] = points[2*ind[1]]
    ret[2*(dim-1)+1] = points[2*ind[1]+1]
    rest = []
    for i in range(ind[0]):
        rest.append(disttab[i][ind[0]])
    for i in range(ind[0]+1, dim):
        if i == ind[1]:
            continue
        rest.append(disttab[ind[0]][i])
    for i in range(dim-2, 0, -1):
        indmax = np.argmax(rest)
        if indmax < ind[0]:
            ret[2*i] = points[2*indmax]
            ret[2*i+1] = points[2*indmax+1]
        elif indmax < ind[1]-1:
            ret[2*i] = points[2*(indmax+1)]
            ret[2*i+1] = points[2*(indmax+1)+1]
        else:
            ret[2*i] = points[2*(indmax+2)]
            ret[2*i+1] = points[2*(indmax+2)+1]
        rest[indmax] = 0  # dezaktywacja
    return ret


def IsElementOfPlane(plane, point):  # sprawdza czy punkt należy do płaszcyzny
    plane = np.transpose(plane)
    # print('plane = ', plane)
    solution = np.linalg.solve(plane, point)
    sum = 0
    for i in range(3):
        sum = sum + solution[i]
    if math.fabs(sum - 1) < epsilon:
        return 1
    else:
        return 0


#
# splane = [[1, -2, 1], [3, 1, -1], [-2, 1, 2]]
# spoint = [2, -8, 2]
# if IsElementOfPlane(splane, spoint):
#     print('tak')
# else:
#     print('nie')
def CoplanarCheck(face, faces):
    # checks if there is a face coplanar with the given
    plane = [faces[face[0]][face[1]][0], faces[face[0]][face[1]][1], faces[face[0]][face[1]][2]]
    plane = np.reshape(plane, (3, 3))
    indx = 0
    ret = []
    for i in range(len(faces)):
        for j in range(len(faces[i])):
            if face[0] == i:
                if face[1] == j:
                    # print('the same face')
                    continue
            for k in range(len(faces[i][j])):
                if IsElementOfPlane(plane, faces[i][j][k]):
                    indx = indx + 1
                else:
                    indx = 0
                    break
            if indx == len(faces[i][j]):
                for l in range(indx-1):
                    ret.append(faces[i][j][l])
                    ret.append(faces[i][j][l+1])
                ret.append(faces[i][j][indx-1])
                ret.append(faces[i][j][0])
                print('coplanar face gen = ', i, 'and num = ', j)
                return ret
    return ret


def CircumradiusOfSnubPoly(config):  # compute circumradius of snub polyhedron
    #  config is the vertex configuration of the snub polyhedron
    p = config[0]
    q = config[2]
    r = config[4]
    density = config[6]

    def ang(vl, pos):
        # print('vlin = ', vl)
        sign = 1
        if pos[0] == 2*pos[1]:  # digon
            return 0
        elif pos[0] < 2*pos[1]:  # retrograde
        #     fract = pos[0]/(pos[0] - pos[1])
            sign = -1
        # else:
        #     fract = pos[0]/pos[1]
        fract = pos[0] / pos[1]
        # print('fract = ', fract)
        vl = 1 - 4 * (1 - math.cos(vl)) * ((math.cos(pi / fract)) ** 2)
        # print('vl = ', vl)
        vl = math.acos(vl)
        vl = sign*vl
        return vl

    def angdiff(vl, pos):
        sign = 1
        if pos[0] == 2 * pos[1]:  # digon
            return 0
        elif pos[0] < 2 * pos[1]:  # retrograde
            fract = pos[0] / (pos[0] - pos[1])
            sign = -1
        else:
            fract = pos[0] / pos[1]
        nom = 4 * ((math.cos(pi / fract)) ** 2) * math.sin(vl)
        denom = math.sqrt((1 - (1 - 4 * (1 - math.cos(vl)) * (math.cos(pi / fract)) ** 2) ** 2))
        ret = sign*nom/denom
        return ret

    def zero(t):
        x = ang(t, p)
        y = ang(t, q)
        z = ang(t, r)
        v = 2 * pi * density - 3 * t - x - y - z
        # v = 4 * pi - 3 * t - x - y - z
        return v

    def zerodiff(t):
        diff = -3 - angdiff(t, p) - angdiff(t, q) - angdiff(t, r)
        return diff

    def newton(first):
        nxt = first - (zero(first) / zerodiff(first))
        # print('next = ', nxt)
        rest = zero(nxt)
        convergence = 0
        while math.fabs(rest) > epsilon:
            nxt = nxt - (zero(nxt) / zerodiff(nxt))
            rest = zero(nxt)
            convergence = convergence + 1
            if convergence > 100:
                print('too many steps')
                break
        print('number of steps = ', convergence)
        return nxt

    solution = newton(1.329)
    ro = 1 / (2 * math.sin(solution / 2))
    print('base of pyramid radius = ', ro)
    R = 1 / (2 * math.sqrt(1 - ro ** 2))
    print('circumradius = ', R)
    return R

# symbol = [[3, 1], [3, 1], [2, 1], [3, 1], [5, 2], [3, 1]]
# symbol = [[3, 1], [3, 1], [2, 1], [3, 1], [5, 3], [3, 1]]
# symbol = [[5, 3], [3, 1], [2, 1], [3, 1], [5, 1], [3, 1]]
# symbol = [[2, 1], [3, 1], [5, 1], [3, 1], [5, 3], [3, 1]]
# symbol = [[3, 1], [3, 1], [5, 1], [3, 1], [5, 3], [3, 1]]
# symbol = [[5, 1], [3, 1], [2, 1], [3, 1], [5, 2], [3, 1]]
# symbol = [[3, 1], [3, 1], [3, 1], [3, 1], [5, 2], [3, 1]]
# radius = CircumradiusOfSnubPoly(symbol)




