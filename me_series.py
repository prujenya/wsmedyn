__author__ = 'jenya'

import statsmodels as sm
import patsy
import pandas as pd
from me_pdbread import pdbname


df = pd.read_csv('filez/output/' + '4x5m_T=320' + '.csv')
#contacts,energy,residues

ct = df['contacts']
en = df['energy']
rs = df['residues']

#print ct
#print en
#print rs


#print df.describe()
print ct.describe()
print en.describe()
print rs.describe()


df1 = pd.DataFrame([ct, ct])

print(df1)

