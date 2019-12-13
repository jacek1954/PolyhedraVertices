import numpy as np
import gsbrouts
import datagosid
import movetoplane


dimgf = [3, 5]  # liczba wierzchołków ścian danego typu

gface = 0  # wybór rodzaju ściany
nface = 44  # i numeru

# tworzenie zbioru punktów do rysowania siatki - wspórzędne wierzchołków ściany i odcinków z FaceMap
rct = gsbrouts.FaceMap(gface, nface)
print('rct = ', rct)
pctcoord = []
for i in range(len(datagosid.face[gface][nface])):  # współrzędne wierzchołków ściany
    vrt = datagosid.face[gface][nface][i]
    vrt = np.reshape(vrt, 3)
    pctcoord.append(vrt)
for i in range((len(rct))):  # współrzędne odcinków mapy
    pctcoord.append(rct[i][0])
    pctcoord.append(rct[i][1])

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
recf = str(dimgf[gface])
templatedata.write(recf)  # zapisujemu na początku pliku liczbę wierzchołków ściany
for i in range(len(convcoord[0])):
    srecf = str(convcoord[0][i])
    templatedata.write(srecf)
templatedata.write('\n')
templatedata.close()
