#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Artyom"

import random
#import math
import matplotlib.pyplot as plt
import me_contacts
import me_energy
import sys
sys.stdout = open('me_dyn_log.txt', 'w')
from me_pdbread import *
from datetime import *
print("Run little bottom program, Run!")
tstart = datetime.now()
print tstart

#from me_contacts import wsme_contacts
#from me_energy import me_energy

#graphlab.canvas.set_target()

# M is the number of time steps
M = 2000

# N is the number of repeated residues
Nr = length


# Temperature
T0 = .25

# "spin" flip probability
P_flip = .5
e_list = []
m_list = []

# вот создал список
list0 = [random.randint(0, 1) for i in range(Nr)]
print "вот наш начальный список без Б    - ", list0, "\n"

listime = range(1, M)

print ("теперь ка пожалуй сгенерим матрицу : \n")

matrix = me_contacts.wsme_contacts(Nr)

#matrix = me_contacts.wsme_contacts_pdb()

    # [[randint(0, 1) for j in range(10)] for i in range(10)]

for r in matrix:
    print r

#print ("\nА сейчас мы случайным образом выберим число из нашей матрицы: \n")
#i = random.randint(0, Nr-1)
#j = random.randint(0, Nr-1)

#print "matrix[", i, "][", j, "]=", matrix[i][j], "\n"

print "initial energy = ", me_energy.me_energy(matrix, list0), "\n"
print "initial moment = ", sum(list0)*Nr**(-1), "\n"
#print ("Run little bottom program, Run!"),tstart
# вот выбирается случайный элемент из всего списка и меняется
e_list.append(me_energy.me_energy(matrix, list0))
m_list.append(sum(list0)*Nr**(-1))

for t in listime:

    while True:

        i = random.randint(0, Nr-1)
        #print "это time - ", t, "\n"

        list1 = list0

        #print "это i - ", i, "\n"
        #print "это i-ый эллемент - ", list0[i], "\n"

        if list1[i] == 1:
            list1[i] = 0
        elif list1[i] == 0:
            list1[i] = 1

        e0 = me_energy.me_energy(matrix, list0)
        e1 = me_energy.me_energy(matrix, list1)
        delta_e = e1 - e0


        P_flip = 1/(1+math.exp(delta_e/T0))

    # Здесь надо добавить проверку для delta_e и ее обработку

        if P_flip > random.random():
            list0 = list1
            e_list.append(me_energy.me_energy(matrix, list0))
            m_list.append(sum(list0)*Nr**(-1))
            break
        else:
            continue


print "energy list = ", e_list
print "moment list = ", m_list
print ("вот конечная конформация - "),  list0, ("\n")
print "final energy = ", me_energy.me_energy(matrix, list0), "\n"
print "final moment = ", sum(list0)*Nr**(-1), "\n"

tend = datetime.now()
delta = tend - tstart
#timeng = int(delta.total_seconds() * 1000)
print("Do you remember the bottom program, that was standing and defending itself? Well, that has changed - "), tend
print "Current run taked = ", delta, "sec."
plt.figure(1)

plt.subplot(211)
plt.plot(e_list)

plt.subplot(212)
plt.plot(m_list,color='r')


plt.show()