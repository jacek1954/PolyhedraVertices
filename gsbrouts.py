import numpy as np
import math
import datagosid

epsilon = 0.0000000000001
fi = (math.sqrt(5) + 1)/2


def wyznacznik(macierz, wymiar):
    # oblicza wyznacznik macierzy
    if wymiar == 1:
        return macierz[0]
    bez1wiersza = np.delete(macierz, 0, 0)
    ret = 0
    for i in range(wymiar):
        minormac = np.delete(bez1wiersza, i, 1)
        ret = ret + ((-1)**i)*macierz[0][i]*wyznacznik(minormac, wymiar - 1)
    return ret


def punkt(pst):
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
            t[i][j] = -wyznacznik(s, 3)
    for i in range(3):
        t[i][2] = -t[i][2]
    # t jest obliczone, tworzymy układ równań i rozwiązujemy
    prawa = np.delete(t, (1, 2, 3), 1)  # prawa jest wektorem po prawej stronie układu równań
    tablica = np.delete(t, 0, 1)  # tablica jest macierzą lewej strony ukłądu równań
    w = wyznacznik(tablica, 3)
    if math.fabs(w) < epsilon:
        raise np.linalg.LinAlgError
    solution = np.zeros((3, 1))
    solution = np.linalg.solve(tablica, prawa)  # linalg rozwiązuje układ równań
    r = np.reshape(solution, 3)
    return r


def dist(p, q):
    # oblicza odległoość punktów p i q
    w = 0
    for i in range(3):
        w = w + (p[i] - q[i])**2
    return math.sqrt(w)





class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, point):
        return self.x == point.x and self.y == point.y and self.z == self.z

    def getCordList(self):
        return [self.x, self.y, self.z]

    def getDistance(self, point):
        # oblicza odległoość punktów p i q
        w = 0
        cords = self.getCordList()
        cordsPoint = point.getCordList()

        for i in range(len(cords)):
            w = w + (cords[i] - cordsPoint[i]) ** 2
        return math.sqrt(w)

class Triangle:
    def __init__(self, a: Point, b: Point, c: Point):
        self.vert = (a, b, c)


a = Point(1, 2, 3)
b = Point(1, 2, 3)
c = [Point(3, 4, 5), Point(3, 4, 5)]

# print('aaa', Point)
# print('aaa', a)
# print('aaa', c[0] == c[1])


def SearchEdge(facein, numin, facescan):
    # SearchEdge znajduje ściany (płaszczyzny) mające wspólną krawędź z daną
    # facein i numin określają rodzaj i numer danej ściany
    # facescan określa rodzaj poszukiwanej ściany, tzn. jest wybraną listą ścian
    # wartościa jest lista zawierająca na przemian 2 wierzchołki i numer ściany:
    # [[[a1 ,b1, c1], [a2, b2, c2]], n1, [[a3, b3, c3], [a4, b4, c4], n2,...]
    # [
    #   [Point, Point],
    #   [Point, Point]
    # ]
    ret = []
    # print('SE', numin, len(facein), len(facescan))
    for i in range(int(len(facescan))):
        lc = 0
        vert = []
        for j in range(len(facein[0])):
            for k in range(len(facescan[0])):
                lcc = 0
                crd = []
                for l in range(3):
                    if facein[numin][j][l] == facescan[i][k][l]:
                        lcc = lcc + 1
                        crd.append(facescan[i][k][l])
                if lcc == 3:
                    lc = lc + 1
                    vert.append(crd)
        if lc == 2:
            ret.append(vert)
            ret.append(i)
    return ret


def SegmIntersect(blue, red):
    print('SIblue = ', blue)
    print('SIred = ', red)
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


def SearchCut(facein, numin, fsctype, edges):
    # searchcut znajduje ściany (płaszczyzny) przecinające daną
    # facein i numin określają rodzaj i numer danej ściany
    # fsctype określa rodzaj poszukiwanej ściany
    # edges zawiera ściany krawędziowe (ale tylko trzech typów) : [etriangle, epentagon, epentagram]
    # każdy element tej listy jest postaci [[[a1 ,b1, c1], [a2, b2, c2]], n1, [[a3, b3, c3], [a4, b4, c4], n2,...]
    # czyli zawiera krawędź, i numer ściany (danego typu)
    ret = []
    # planein, planeedge i planescan są płaszczyznami, których przecięcie szukamy
    # print(numin, fsctype)
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
                    planeedge[i][j] = datagosid.face[nted][k][i][j]
            for nsc in range(len(datagosid.face[fsctype])):  # nsc wybiera kolejną ścianę potencjalnie poszukiwaną
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
                            planescan[i][j] = datagosid.face[fsctype][nsc][i][j]
                    p = punkt([planein, planeedge, planescan])  # obliczenie punktu wspólnego płaszczyzn
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


def EdgeBuild(intype, numface):
    print('EB', intype, numface)
    # tworzenie listy ścian krawędziowych dla danej ściany
    # tu szukane są tylko trójkąty, pięciokąty i pentagramy:
    edgtriangle = SearchEdge(datagosid.face[intype], numface, datagosid.face[0])
    edgpentagram = SearchEdge(datagosid.face[intype], numface, datagosid.face[1])
    edg = [edgtriangle, edgpentagram]
    # print('edg = ', edg)
    return edg


def FaceMap(fintype, finnum):  # tu trzeba uwzglednić zmiany w SearchCut
    # FaceMap tworzy listę odcinków, które są częściami wspólnymi danej ściany ze ścianami ją przecinającymi
    # fintype określa rodzaj danej ściany, a finnum jej numer
    # wartością funkcji jest lista, której każdy element skłąda się z czterech danych:
    # [typ ściany przecinającej, numer tej ściany, współrzędne odcinka wspólnego, długość tego odcinka]
    # ret = [[fct, fcn, seg, lngseg], [ , , , ], [ , , , ], ...]
    ret = []
    print(fintype, finnum)
    ev1 = EdgeBuild(fintype, finnum)  # krawędzie danej ściany
    # print(ev1[0])
    # print(len(ev1[0]))
    # print(ev1[1])
    # print(ev1[2])
    z = []  # lista ścian (płaszczyzn) przecinających daną
    for i in range(len(datagosid.face)):  # są tylko trzy rodzaje ścian
        z.append(SearchCut(datagosid.face[fintype], finnum, i, ev1))
    # print('FM')
    # print(len(z[0])/5)
    # print(z[0])
    # print(len(z[1])/5)
    # print(z[1])
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
            print('n = ', n, 'm = ', m, 'nmf = ', nmf)
            ev2 = EdgeBuild(n, nmf)
            # print('ev2 = ', ev2)
            y = []
            for l in range(len(datagosid.face)):
                y.append(SearchCut(datagosid.face[n], nmf, l, ev2))
            # teraz trzeba sprawdzić czy y zawiera finnum dla fintype
            # a jesli tak to utworzyć listy blue i red
            # print('y = ', y)
            blue = []
            for j1 in range(int(len(y[fintype])/5)):
                if y[fintype][5*j1+3] == finnum:
                    blue.append(y[fintype][5*j1])
                    blue.append(y[fintype][5*j1+4])
            if blue == []:
                # print('blue jest pusty')
                skip.append(nmf)
                continue
            # czyli blue jest niepusty, tworzymy red, który pusty nie jest
            # print('blue = ', blue)
            red = []
            for j2 in range(m, int(len(z[n])/5)):
                if z[n][5*j2+3] == nmf:
                    red.append(z[n][5*j2])
                    red.append(z[n][5*j2+4])
            skip.append(nmf)
            red = DuplicatesReduct(red)
            blue = DuplicatesReduct(blue)
            print('Redblue = ', blue)
            print('Redred = ', red)
            if len(red) <= 2:
                continue
            if len(blue) <= 2:
                continue
            # teraz zamieniamy odległości na typy
            genf = PolygonType(fintype)
            genc = PolygonType(n)
            # print('genf = ', genf,'genc = ', genc)
            for j3 in range(int(len(red)/2)):
                red[2*j3+1] = PointType(genf, red[2*j3+1])
            for j4 in range(int(len(blue)/2)):
                blue[2*j4+1] = PointType(genc, blue[2*j4+1])
            # teraz trzeba uporządkowac punkty, ale tylko dla gwiazd
            if genf >= 7:
                red = PointSort(red)
                print('redsort = ', red)
                red = StellaTrace(red)
                print('stellared = ', red)
            else:
                gbg = red.pop(1)
                gbg = red.pop(2)
                red = [red]
                print('Nosortred = ', red)
            if genc >= 7:
                # print(genc)
                # print(blue)
                blue = PointSort(blue)
                print('bluesort', blue)
                blue = StellaTrace(blue)
                print('stellablue', blue)
            else:
                gbg = blue.pop(1)
                gbg = blue.pop(2)
                blue = [blue]
                print('Nosortblue = ', blue)
            for i1 in range(len(blue)):
                for j5 in range(len(red)):
                    rt = SegmIntersect(blue[i1], red[j5])
                    print('SIreturn = ', rt)
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
    # print(fgen, dist)
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


def PolygonType(t):  # funkcja tymczasowa dla gosid
    if t == 0:
        return 1
    elif t == 1:
        return 7


# funkcja wyznacza zbiór odcinków powstałych przez przecięcie ściany gwiaździstej przez płaszczyznę określone punktami
# przecięcia z krawędziami gwiazdy podanymi w tablicy point, code zawiera kody typów punktów z point
# lista punktów zawartych w point jest uporządkowana w przestrzeni
def StellaTrace(point):
    print('Stellain = ', point)
    l = len(point)-1  # pozycja ostatniego punktu
    zero = []  # odcinek zerowy
    s = []
    if point[1] == 'b':  # pierwszy punkt jest graniczny - ślad jest jednym odcinkiem
        return [point[1], point[l]]
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
    while i <= l:
        point[k] = 'n'  # 'deaktywacja' początku ('e')
        try:  # znajdź koniec odcinka
            i = point.index('e')
            s.append([point[k-1], point[i-1]])
            point[i] = 'n'  # 'deaktywacja' końca ('e')
        except ValueError:
            s.append([point[k-1], point[l]])  # nie ma już końca typu 'e' - to ostatni odcinek
            return s
        try:  # znajdź początek odcinka
            k = point.index('e')
        except ValueError:
            return s  # nie ma więcej odcinków


def PointSort(points):
    # print('points = ', points)
    dim = int(len(points)/2)
    disttab = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i+1, dim):
            disttab[i][j] = dist(points[2*i], points[2*j])
    indmax = np.argmax(disttab)
    ind = np.unravel_index(indmax, disttab.shape)
    # print(disttab)
    # print(indmax)
    # print(ind)
    ret = []
    for i in range(dim):
        ret.append([0, 0, 0])
        ret.append(0)
    ret[0] = points[2*ind[0]]
    ret[1] = points[2*ind[0]+1]
    ret[2*(dim-1)] = points[2*ind[1]]
    ret[2*(dim-1)+1] = points[2*ind[1]+1]
    # print('ret = ', ret)
    rest = []
    for i in range(ind[0]):
        rest.append(disttab[i][ind[0]])
    for i in range(ind[0]+1, dim):
        if i == ind[1]:
            continue
        rest.append(disttab[ind[0]][i])
    # print('rest = ', rest)
    for i in range(dim-2, 0, -1):
        # print('i = ', i)
        indmax = np.argmax(rest)
        # print('indmax = ', indmax)
        if indmax < ind[0]:
            ret[2*i] = points[2*indmax]
            ret[2*i+1] = points[2*indmax+1]
            # print('ret = ', ret)
        elif indmax < ind[1]-1:
            ret[2*i] = points[2*(indmax+1)]
            ret[2*i+1] = points[2*(indmax+1)+1]
            # print('ret = ', ret)
        else:
            ret[2*i] = points[2*(indmax+2)]
            ret[2*i+1] = points[2*(indmax+2)+1]
            # print('ret = ', ret)
        rest[indmax] = 0  # dezaktywacja
    return ret


