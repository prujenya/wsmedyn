#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'jenya'

from me_pdbread import atoms_XYZ
import math

cut_length = float(4)

# pythran export wsme_contacts(int)
def wsme_contacts(nn):

    from random import randint
    import numpy

    z = 3
    matrix = numpy.zeros(shape=(nn,nn))

    for j in range(nn):
        for i in range(nn):
            if randint(1, .5*nn) == 1:
                matrix[i][j] = 1
            else:
                matrix[i][j] = 0

    #matrix = [[randint(0, 1) for j in range(nn)] for i in range(nn)]

    return matrix


def wsme_contacts_pdb():

    from me_pdbread import length
    from me_pdbread import dist_mtx
    import numpy

    matrix = numpy.zeros(shape=(length,length))

    #matrix = []

    for i in range(length):
        for j in range(length):

            if dist_mtx[i][j] < cut_length:
                matrix[i][j]=1
            else:
                matrix[i][j]=0

    return matrix



def wsme_contacts_hairpin(L):


    import numpy

    matrix = numpy.zeros(shape=(L,L))

    #matrix = []

    for i in range(L):
        for j in range(L):

            if i+j == L:
                matrix[i][j]=1
            else:
                matrix[i][j]=0

    return matrix


def wsme_contacts_helix(L, Delta):

    import numpy

    matrix = numpy.zeros(shape=(L,L))

    #matrix = []

    for i in range(L):
        for j in range(L):

            if math.fabs(j-i) == Delta:
                matrix[i][j]=1
            else:
                matrix[i][j]=0

    return matrix