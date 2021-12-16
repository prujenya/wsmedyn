#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'jenya'

from functools import reduce
import numpy as np


def me_extension_energy(force, persist, T, config):
    my_length = contour_length(config)
    a = force * persist / T
    inv_ext = inv_extension([4, -(9 + a), 2 * (3 + a), -a])

    return my_length*inv_ext


def contour_length(config):
    mysum = 0

    L = len(config)

    for i in range(L-1):
        for j in range(i + 1, L):
            temp = config[i-1:j]  # ??? Предполагается, что это часть списка с i-го по j-й элемент. так ли это???
            mysum = mysum + config[i - 1] * reduce(lambda res, x: res * x, temp, 1) * config[j]

    return mysum

def inv_extension(coeff):
    myroots = np.roots(coeff)
    myroot = np.max(myroots[myroots.imag == 0])

    print('my root=', myroot)

    return myroot

