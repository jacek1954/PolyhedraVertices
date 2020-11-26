import numpy as np
import sbrouts
import movetoplane
import coordinates


def BowersName():
    # name = 'siddid'
    name = 'sided'
    # name = 'qrid'
    # name = 'gaquatid'
    # name = 'girsid'
    return name


def FaceId():
    gen = 0
    num = 0
    return [gen, num]


if BowersName() == 'sided':
    faces = coordinates.sided()
elif BowersName() == 'gosid':
    faces = coordinates.gosid()
elif BowersName() == 'gaquatid':
    faces = coordinates.gaquatid()
elif BowersName() == 'did':
    faces = coordinates.did()
elif BowersName() == 'seside':
    faces = coordinates.seside()
elif BowersName() == 'siddid':
    faces = coordinates.siddid()
elif BowersName() == 'qrid':
    faces = coordinates.qrid()
elif BowersName() == 'girsid':
    faces = coordinates.girsid()
else:
    faces = []


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

