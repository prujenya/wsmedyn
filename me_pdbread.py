#!/usr/bin/env python
# -*- coding: utf-8 -*-

#__author__ = 'jenya'


#type = 'CA'
    #raw_input('Input Type for search:')
#for line in open('1AKI.pdb'):#
    #list = line.split()
    #id = list[0]
    #if id == 'ATOM':
        #if type in list[2] :
			#print "Nomer  	|TYPE   | Ostatok | Chain |    Atom_N   |   X	     Y	     Z		|"
			#print " "
			#print list[1], "	|", list[2], "	|", list[3], "	  |", list[4], "	  |","	", list[5], "	|", list[6],"", list[7],"", list[8],"", "	|"
			#print " "

################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'artyom'

import math

atoms_XYZ = []

atoms_X = []
atoms_Y = []
atoms_Z = []

dist_mtx = []

type = 'CA'
pdbname = '1g2r' #'1vku'#'1g2r' #'2efv'#''1mn8' #'helix'#'4x5m'

if len(pdbname) == 4:
    for line in open('filez/'+ pdbname + '.pdb'):#1AKI#1g2r  # 1vku
        list = line.split()
        id = list[0]
        if id == 'ATOM':
            if type in list[2]:
                atoms_X.append(float(list[6]))
                atoms_Y.append(float(list[7]))
                atoms_Z.append(float(list[8]))
            #atoms_XYZ.append((list[6], list[7], list[8]))

#atoms_XYZ[1] = map(int, atoms_XYZ)
    atoms_XYZ.append(atoms_X)
    atoms_XYZ.append(atoms_Y)
    atoms_XYZ.append(atoms_Z)

    length = len(atoms_X)

    dist_mtx = [[math.sqrt((atoms_X[i]-atoms_X[j])**2+(atoms_Y[i]-atoms_Y[j])**2+(atoms_Z[i]-atoms_Z[j])**2) for j in range(length)] for i in range(length)]
else:
    length = 0

#print ("матрица расстояний : \n")
#for rr in dist_mtx:
 #   print rr

#print atoms_XYZ