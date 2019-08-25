import numpy as np
import sbrouts
import math
import danesided


def searchcut(facein, numin, fsctype, edges):
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
                    p = sbrouts.punkt([planein, planeedge, planescan])  # obliczenie punktu wspólnego płaszczyzn
                    q1 = sbrouts.dist(p, edges[nted][2*ned][0])
                    q2 = sbrouts.dist(p, edges[nted][2*ned][1])
                    if math.fabs(q1 + q2 - 1) <= sbrouts.epsilon:         # zbadanie czy punkt leży na krawędzi
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


def segmintersect(evin, evct, fintype, finnum, fctype, fcnum):
    # oblicza współrzędne odcinka będącego częścią wspólną danych ścian
    # fintype, finnum, fctype, fcnum określają typ i numer danych ścian
    # evin i evct  współrzedne punktu przecięcia, rodzaj krawędzi oraz numery krawędzi i ściany przecinającej
    #     a także nie wiekszą z odległości punktu od końca krawędzi
    # [[a1, b1, c1], et1, en1, ct1, d1, [a2, b2, c2], et2, en2, ct2, d2,...]
    blue = []
    red = []
    for i in range(int(len(evin[fctype])/5)):  # blue zawiera punkty przecięcia ściany "in" przez płaszczyznę "ct"
        if evin[fctype][5*i+3] == fcnum:       # zakładamy, że te pumkty istnieją
            blue.append(evin[fctype][5*i])
            blue.append(evin[fctype][5*i+4])
    for i in range(int(len(evct[fintype])/5)):  # red zawiera punkty przecięcia ściany "ct" przez płaszczyznę "in"
        if evct[fintype][5*i+3] == finnum:      # ale takie punkty moga nie istnieć!!!
            red.append(evct[fintype][5*i])
            red.append(evct[fintype][5*i+4])
    if len(red) == 0:  # red nie istnieje
        return [[0, 0, 0], [0, 0, 0]]
    lengblue = sbrouts.dist(blue[1], blue[0])
    lengred = sbrouts.dist(red[0], red[1])
    dsttab = []
    dsttab.append(lengblue)
    dsttab.append(lengred)
    dsttab.append(sbrouts.dist(blue[0], red[0]))
    dsttab.append(sbrouts.dist(blue[1], red[1]))
    dsttab.append(sbrouts.dist(blue[0], red[1]))
    dsttab.append(sbrouts.dist(blue[1], red[0]))
    indmax = np.argmax(dsttab)
    res = []
    if indmax == 0:
        res.append(red[0])
        res.append(red[1])
    elif indmax == 1:
        res.append(blue[0])
        res.append(blue[1])
    elif indmax == 2:
        if lengblue+lengred-sbrouts.dist(blue[0], red[0]) >= sbrouts.epsilon:  # badamy, czy iloczyn red i blue jest niepusty
            res.append(blue[1])
            res.append(red[1])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 3:
        if lengblue+lengred-sbrouts.dist(blue[1], red[1]) >= sbrouts.epsilon:
            res.append(blue[0])
            res.append(red[0])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 4:
        if lengblue+lengred-sbrouts.dist(blue[0], red[1]) >= sbrouts.epsilon:
            res.append(blue[1])
            res.append(red[0])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    elif indmax == 5:
        if lengblue+lengred-sbrouts.dist(blue[1], red[0]) >= sbrouts.epsilon:
            res.append(blue[0])
            res.append(red[1])
        else:
            return [[0, 0, 0], [0, 0, 0]]
    return res


# funkcja wyznacza zbiór odcinków powstałych przez przecięcie ściany gwiaździstej przez płaszczyznę określone punktami
# przecięcia z krawędziami gwiazdy podanymi w tablicy point, code zawiera kody typów punktów z point
# lista punktów zawartych w point jest uporządkowana w przestrzeni
def stellatrace(point, code):
    l = len(point)-1
    zero = [[0, 0, 0], [0, 0, 0]]  # odcinek punktowy
    if len(point) <= 1:
        return zero
    s = []
    if code[0] == 'b':
        return [point[0], point[l]]
    if code[0] == 'v':
        if code[1] != 'e':
            if code[1] == 'v':
                return zero
            else:
                return [point[0], point[l]]
    if code[0] == 'e':
        k = 0
    elif code[0] == 'v':
        k = 1
    else:
        print('niedopuszczalny znak typu punktu')
        return zero
    i = k + 1
    while i <= l:
        code[k] = 'n'
        try:  # znajdź koniec odcinka
            i = code.index('e')
            s.append([point[k], point[i]])
            code[i] = 'n'
        except ValueError:
            s.append([point[k], point[l]])
            return s
        try:  # znajdź początek odcinka
            k = code.index('e')
        except ValueError:
            return s


p = [0, 1, 2, 3, 4, 5]
c = ['s', 's', 's', 'i', 'i', 's']
seg = stellatrace(p, c)
print(seg)
print(c)
p = [0, 1, 2]
c = ['e', 'b', 'e']
seg = stellatrace(p, c)
print(seg)
print(c)
p = [0, 1, 2, 3, 4]
c = ['e', 'i', 'i', 'e', 'v']
seg = stellatrace(p, c)
print(seg)
print(c)
p = [0, 1, 2, 3]
c = ['e', 'b', 'e', 'v']
seg = stellatrace(p, c)
print(seg)
print(c)
p = [0, 1]
c = ['v', 'v']
seg = stellatrace(p, c)
print(seg)
print(c)

