import numpy as np
import sbrouts
import danesided
import movetoplane


dimgf = [3, 5, 5]  # liczba wierzchołków ścian danego typu

gface = 0  # wybór rodzaju ściany
nface = 0  # i numeru

# tworzenie zbioru punktów do rysowania siatki - wspórzędne wierzchołków ściany i odcinków z facemap
rct = sbrouts.facemap(gface, nface)
pctcoord = []
for i in range(len(danesided.face[gface][nface])):  # współrzędne wierzchołków ściany
    vrt = danesided.face[gface][nface][i]
    vrt = np.reshape(vrt, 3)
    pctcoord.append(vrt)
for i in range(int(len(rct)/4)):  # współrzędne odcinków mapy
    pctcoord.append(rct[4*i+2][0])
    pctcoord.append(rct[4*i+2][1])

# przeniesienie pctcoord na płaszzczyznę xy - jako convcoord
convcoord = movetoplane.movetoxyplane(pctcoord)

# zaokrąglenie danych do 5 miejsc po przecinku
for i in range(len(convcoord[0])):
    for j in range(3):
        convcoord[0][i][j] = round(convcoord[0][i][j], 5)

# przejście do 2 wymiaru - wyrzucenie niepotrzebnej trzeciej współrzędnej (zerowej)
for i in range(len(convcoord[0])):
    convcoord[0][i] = np.delete(convcoord[0][i], 2, 0)

# tworzenie pliku ze współrzędnymi
templatedata = open('danetest', 'w')
recf = str(dimgf[gface])
templatedata.write(recf)  # zapisujemu na początku pliku liczbę wierzchołków ściany
for i in range(len(convcoord[0])):
    srecf = str(convcoord[0][i])
    templatedata.write(srecf)
templatedata.close()
