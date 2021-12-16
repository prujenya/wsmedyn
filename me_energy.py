#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'jenya'

from functools import reduce

import me_extension as ext

persist = 1

def me_energy(contacts, config, entropies, T, epsilon):
    mysum = 0

    L = len(contacts)

    for i in range(L-1):  # L-1
        for j in range(i + 1, L):  # L
            temp = config[i-1:j]  # config[i-1:j]  ??? Предполагается, что это часть списка с i-го по j-й элемент. так ли это???
            mysum = mysum + epsilon * contacts[i][j] * reduce(lambda res, x: res * x, temp, 1)

    for i in range(L):
        mysum = mysum - T * entropies[i] * config[i]

    alpha = 0
    if -1 != epsilon:
        alpha = 1

    return mysum + alpha*ext.me_extension_energy(10, persist, T, config)
