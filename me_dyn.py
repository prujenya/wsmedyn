#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'jenya'


import random
#import copy
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
#import graphlab

import me_contacts
import me_energy
import me_io_data as io

from me_pdbread import length

from me_pdbread import pdbname
import me_extension as ext


#from me_contacts import wsme_contacts
#from me_energy import me_energy

#graphlab.canvas.set_target()
start = time.time()
# universal gas constant
R =1.986 # cal/(K*mol)
epsilon0 = -3.0*550 # cal/mol
dS = -2.6 #(3.9+1.3)/2 #-3.1  # cal/(K*mol)
# M is the number of time steps
M = 550# 550 # 750 #750 # 1600 # 200 #3000

# N is the number of repeated residues
Nr = 85

Delta = 3


# Temperature
T0 = 350 # 300 # K

#TempList = [20, 40, 60, 80] # random contact matrix
TempList = [285, 320, 350, 380] ### pdb
#TempList = [380, 400, 420, 440]
#TempList = [170, 200, 230, 260] ### hairpin



# "spin" flip probability
P_flip = .5
e_list = np.zeros(shape=(len(TempList), M)) #empty
m_list = np.zeros(shape=(len(TempList), M))
contacts_list = np.zeros(shape=(len(TempList), M))




listime = np.array(range(1, M))

print ("теперь ка пожалуй сгенерим матрицу : \n")

if len(pdbname) == 4:

    Nr = length
    matrix = me_contacts.wsme_contacts_pdb()##############

elif pdbname == 'helix':

    matrix = me_contacts.wsme_contacts_helix(Nr,Delta)

elif pdbname == 'hairpin':

    matrix = me_contacts.wsme_contacts_hairpin(Nr)

elif pdbname == 'random':

    matrix = me_contacts.wsme_contacts(Nr)


entropies0 = [0 for i in range(Nr)]
entropies = [dS for i in range(Nr)]

contacts_max = -me_energy.me_energy(matrix, [1 for i in range(Nr)], entropies0, T0, -1)
    # [[randint(0, 1) for j in range(10)] for i in range(10)]

print("maximal number of native contacts = ", contacts_max, "\n")

for r in matrix:
    print(r)


for T0 in TempList:

    idx = TempList.index(T0)

    np.append(e_list, [])
    np.append(m_list, [])
    np.append(contacts_list, [])


# вот создал список
    list0 = np.array([random.randint(0, 1) for i in range(Nr)], dtype='int64')
    list1 = np.empty(shape=(Nr,), dtype='int64')
    print("вот наш начальный список без Б    - ", list0, "\n")

# initialize lists: energy, "helicity",contacts
    try:
        """
        np.insert(e_list[idx], 0, me_energy.me_energy(matrix, list0, entropies, T0, epsilon0))
        np.insert(m_list[idx], 0, sum(list0)*Nr**(-1))
        np.insert(contacts_list[idx], 0, -contacts_max**(-1)*me_energy.me_energy(matrix, list0, entropies0, T0, -1))
        """
        e_list[idx, 0] = me_energy.me_energy(matrix, list0, entropies, T0, epsilon0)
        m_list[idx, 0] = sum(list0)*Nr**(-1)
        contacts_list[idx, 0] = -contacts_max**(-1)*me_energy.me_energy(matrix, list0, entropies0, T0, -1)
    except Exception as inst:
        print(type(inst))
        print('idx =', idx)



    for t in listime:

        while True:

            i = random.randint(0, Nr-1)
            #print "это time - ", t, "\n"

            np.copyto(list1, list0) #copy.copy(list0)

            #print "это i - ", i, "\n"
            #print "это i-ый эллемент - ", list0[i], "\n"

            if list1[i] == 1:
                list1[i] = 0
            elif list1[i] == 0:
                list1[i] = 1

            e0 = me_energy.me_energy(matrix, list0, entropies, T0, epsilon0)
            e1 = me_energy.me_energy(matrix, list1, entropies, T0, epsilon0)
            delta_e = e1 - e0

            try:
                a = math.exp(delta_e/(R*T0))
            except OverflowError:
                a = float('inf')

            P_flip = 1/(1+a)



            # Здесь надо добавить проверку для delta_e и ее обработку

            if P_flip > random.random(): # >
                np.copyto(list0, list1) #copy.copy(list1)
                """
                np.copyto(e_list[idx], np.insert(e_list[idx], t, me_energy.me_energy(matrix, list0, entropies, T0, epsilon0)))
                np.insert(m_list[idx], t, sum(list0)*Nr**(-1))
                np.insert(contacts_list[idx], t, -contacts_max**(-1)*me_energy.me_energy(matrix, list0, entropies0, T0, -1))
                """
                e_list[idx, t] = me_energy.me_energy(matrix, list0, entropies, T0, epsilon0)
                m_list[idx, t] = sum(list0) * Nr ** (-1)
                contacts_list[idx, t] = -contacts_max ** (-1) * me_energy.me_energy(matrix, list0, entropies0, T0, -1)

                break
            else:
                """
                np.insert(e_list[idx], t, me_energy.me_energy(matrix, list0, entropies, T0, epsilon0))
                np.insert(m_list[idx], t, sum(list0) * Nr ** (-1))
                np.insert(contacts_list[idx], t, -contacts_max ** (-1) * me_energy.me_energy(matrix, list0, entropies0, T0,-1))
                """
                e_list[idx, t] = me_energy.me_energy(matrix, list0, entropies, T0, epsilon0)
                m_list[idx, t] = sum(list0) * Nr ** (-1)
                contacts_list[idx, t] = -contacts_max ** (-1) * me_energy.me_energy(matrix, list0, entropies0, T0, -1)

                continue

    print("exec. time: ", str(time.time() - start))
    print("temperature: ", str(T0))
    print("вот конечная конформация - ",  list0, "\n")
    print("final energy = ", me_energy.me_energy(matrix, list0, entropies, T0, epsilon0), "\n")
    print("final moment = ", sum(list0)*Nr**(-1), "\n")
    print("final number of contacts = ", -me_energy.me_energy(matrix, list0, entropies0, T0,-1), "\n")

    ###print("listime length = ", listime.size, "\n") #########???
    ###print("e_list length = ", len(e_list[idx]), "\n") #########???
    ###print("contacts matrix: ", matrix)

    plt.figure(1)

    plt.Text(text='tempr=')

    ex = plt.subplot(223)
    ex.set_title("energy")
    plt.plot(e_list[idx], linewidth=2.5, label=T0,)

    mx = plt.subplot(221)
    mx.set_title("residues")
    plt.plot(m_list[idx], linewidth=2.5, label=T0)


    cx = plt.subplot(222)
    cx.set_title("contacts")
    plt.plot(contacts_list[idx], linewidth=2.5, label=T0)


plt.legend(loc="lower right")
plt.show()
