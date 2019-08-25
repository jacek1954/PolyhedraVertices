import numpy as np
import math
import danesided

epsilon = 0.0000000000001


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


def searchedge(facein, numin, facescan):
    # searchedge znajduje ściany (płaszczyzny) mające wspólną krawędź z daną
    # facein i numin określają rodzaj i numer danej ściany
    # facescan określa rodzaj poszukiwanej ściany, tzn. jest wybraną listą ścian
    # wartościa jest lista zawierająca na przemian 2 wierzchołki i numer ściany:
    # [[[a1 ,b1, c1], [a2, b2, c2]], n1, [[a3, b3, c3], [a4, b4, c4], n2,...]
    ret = []
    for i in range(len(facescan)):
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


def segmintersect(evin, evct, fintype, finnum, fctype, fcnum):
    # oblicza współrzędne odcinka będącego częścią wspólną danych ścian
    # fintype, finnum, fctype, fcnum określają typ i numer danych ścian
    # evin i evct  zawierają współrzedne punktu przecięcia oraz numery krawędzi i ściany przecinającej
    # [[ina1, inb1, inc1], ine1, inct1, [ina2, inb2, inc2], ine2, inct2, ...]
    # [[cta1, ctb1, ctc1], cte1, ctct1, [cta2, ctb2, ctc2], cte2, ctct2, ...]
    blue = []
    red = []
    for i in range(int(len(evin[fctype])/4)):
        if evin[fctype][4*i+3] == fcnum:
            blue.append(evin[fctype][4*i])
    for i in range(int(len(evct[fintype])/4)):
        if evct[fintype][4*i+3] == finnum:
            red.append(evct[fintype][4*i])
    if len(red) == 0:  # odcinek red nie istnieje
        return [[0, 0, 0], [0, 0, 0]]
    lengblue = dist(blue[1], blue[0])
    lengred = dist(red[0], red[1])
    dsttab = []
    dsttab.append(lengblue)
    dsttab.append(lengred)
    dsttab.append(dist(blue[0], red[0]))
    dsttab.append(dist(blue[1], red[1]))
    dsttab.append(dist(blue[0], red[1]))
    dsttab.append(dist(blue[1], red[0]))
    indmax = np.argmax(dsttab)
    res = []
    if indmax == 0:
        res.append(red[0])
        res.append(red[1])
    elif indmax == 1:
        res.append(blue[0])
        res.append(blue[1])
    elif indmax == 2:
        if lengblue+lengred-dist(blue[0], red[0]) >= epsilon:  # badamy, czy iloczyn red i blue jest niepusty
            res.append(blue[1])
            res.append(red[1])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 3:
        if lengblue+lengred-dist(blue[1], red[1]) >= epsilon:
            res.append(blue[0])
            res.append(red[0])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 4:
        if lengblue+lengred-dist(blue[0], red[1]) >= epsilon:
            res.append(blue[1])
            res.append(red[0])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 5:
        if lengblue+lengred-dist(blue[1], red[0]) >= epsilon:
            res.append(blue[0])
            res.append(red[1])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    return res


def searchcut(facein, numin, fsctype, edges):
    # searchcut znajduje ściany (płaszczyzny) przecinające daną
    # facein i numin określają rodzaj i numer danej ściany
    # fsctype określa rodzaj poszukiwanej ściany, tzn. fsctype jest numerem wybranej listy ścian
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
                    planeedge[i][j] = danesided.face[nted][k][i][j]
            for nsc in range(len(danesided.face[fsctype])):  # nsc wybiera kolejną ścianę potencjalnie poszukiwaną
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
                            planescan[i][j] = danesided.face[fsctype][nsc][i][j]
                    p = punkt([planein, planeedge, planescan])  # obliczenie punktu wspólnego płaszczyzn
                    q1 = dist(p, edges[nted][2*ned][0])
                    q2 = dist(p, edges[nted][2*ned][1])
                    if math.fabs(q1 + q2 - 1) <= epsilon:         # zbadanie czy punkt leży na krawędzi
                        ret.append(p)
                        ret.append(nted)
                        ret.append(k)
                        ret.append(nsc)
                except np.linalg.LinAlgError:  # jeśli nie ma (jedynego) punktu przeciecia pomiń tę ścianę
                    continue
    # ret zawiera współrzedne punktu przecięcia, rodzaj krawędzi oraz numery krawędzi i ściany przecinającej
    # [[a1, b1, c1], et1, en1, ct1, [a2, b2, c2], et2, en2, ct2, ...]
    return ret


def edgebuild(intype, numface):
    # tworzenie listy ścian krawędziowych dla danej ściany
    # tu szukane są tylko trójkąty, pięciokąty i pentagramy:
    edgtriangle = searchedge(danesided.face[intype], numface, danesided.face[0])
    edgpentagon = searchedge(danesided.face[intype], numface, danesided.face[1])
    edgpentagram = searchedge(danesided.face[intype], numface, danesided.face[2])
    edg = [edgtriangle, edgpentagon, edgpentagram]
    return edg


def facemap(fintype, finnum):  # tu trzeba uwzglednić zmiany w searchcut
    # facemap tworzy listę odcinków, które są częściami wspólnymi danej ściany ze ścianami ją przecinającymi
    # fintype określa rodzaj danej ściany, a finnum jej numer
    # wartością funkcji jest lista, której każdy element skłąda się z czterech danych:
    # [typ ściany przecinającej, numer tej ściany, współrzędne odcinka wspólnego, długość tego odcinka]
    # ret = [[fct, fcn, seg, lngseg], [ , , , ], [ , , , ], ...]
    ret = []
    ev1 = edgebuild(fintype, finnum)  # krawędzie danej ściany
    z = []  # lista ścian (płaszcyzn) przecinających daną
    for i in range(3):  # są tylko trzy rodzaje ścian
        z.append(searchcut(danesided.face[fintype], finnum, i, ev1))
    for n in range(3):
        fctype = n
        fcnumtab = []  # tablica numerów ścian przecinających
        for m in range(int(len(z[n])/4)):
            nmf = z[n][4*m+3]
            count = 1
            for k in range(4*(m+1)+3, int(len(z[n])), 4):
                if z[n][k] == nmf:
                    count = count + 1
            if count >= 2:  # czyli płaszczyzna przecina ścianę w co najmniej dwóch krawędziach!
                fcnumtab.append(nmf)
        # tworzenie tablicy odcinków dla elementów fcnumtab
        for m in range(len(fcnumtab)):
            fcnum = fcnumtab[m]
            ev2 = edgebuild(fctype, fcnum)
            y = []
            for i in range(3):
                y.append(searchcut(danesided.face[fctype], fcnum, i, ev2))
            seg = segmintersect(z, y, fintype, finnum, fctype, fcnum)
            seg = np.reshape(seg, (2, 3))
            lngseg = dist(seg[1], seg[0])
            if lngseg < epsilon:
                continue
            ret.append(fctype)
            ret.append(fcnum)
            ret.append(seg)
            ret.append(lngseg)
    return ret
