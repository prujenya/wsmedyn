__author__ = 'jenya'

import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)


def me_createDF(energylst, residuelst, contactlst):#

    data = {'energy': energylst,'residues': residuelst,'contacts': contactlst} #

    return pd.DataFrame(data)


