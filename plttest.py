# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

plt.axis([0, 1, 0, 1], 'scaled')
crd = open('danetest', 'r')
sig = '0123456789.-'

# czyta kolejną liczbę
def readstr(f, patt):
    st = ''
    ps = -1
    while ps == -1:
        sg = f.read(1)
        ps = patt.find(sg)
    while ps >= 0:
        st = st + sg
        sg = f.read(1)
        ps = patt.find(sg)
    return st

# rysuje odcinek o danych końcach
def drawline(p, q):
    plt.plot([p[0], q[0]], [p[1], q[1]])
    return 0


# wczytanie rodzaju ściany
s = readstr(crd, sig)
d = int(s)

# rysowanie pierwszego boku ściany
x1 = float(readstr(crd, sig))
y1 = float(readstr(crd, sig))
x2 = float(readstr(crd, sig))
y2 = float(readstr(crd, sig))
x0 = x1
y0 = y1
n = drawline([x1, y1], [x2, y2])
# rysowanie kolejnych boków
for i in range(d-2):
    x1 = x2
    y1 = y2
    x2 = float(readstr(crd, sig))
    y2 = float(readstr(crd, sig))
    n = drawline([x1, y1], [x2, y2])
# rysowanie ostatniego boku
x1 = x2
y1 = y2
x2 = x0
y2 = y0
n = drawline([x1, y1], [x2, y2])

# rysowanie mapy cięć
s =  readstr(crd, sig)
while s != '':
    x1 = float(s)
    y1 = float(readstr(crd, sig))
    x2 = float(readstr(crd, sig))
    y2 = float(readstr(crd, sig))
    n = drawline([x1, y1], [x2, y2])
    s =  readstr(crd, sig)
crd.close()
plt.show()
