import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

data = pd.read_csv('MinnesotaJordan.tsv', sep='\t')

# To metric
data['Q'] /= 3.28**3
data['Stage'] /= 3.28


#h_bank = 7. # meters
b = 70. # meters
S = 1E-4

"""
def _manning(h, k_Qbank, P_Qbank, stage_depth_Q_offset):
    Q_ch = b * h**(5/3.) * S**(1/2.)
    if h > h_bank:
        Q_ob = k_Qbank * (h - h_bank)**(P_Qbank)
    else:
        Q_ob = 0
    return Q_ch + Q_ob + stage_depth_Q_offset
"""

def _manning(h, n, k_Qbank, P_Qbank, stage_depth_Q_offset, h_bank):
    return ( b * h**(5/3.) * S**(1/2.) / n ) \
           + (h > h_bank) * k_Qbank * (h - h_bank)**(P_Qbank * (h > h_bank)) \
           + stage_depth_Q_offset

#_h = np.arange(0.,10.)
#_Q = _manning(_h, 1E2, 8/3., 0.)

#plt.plot(_Q, _h, 'k.')

popt, pcov = curve_fit( _manning, data['Stage'], data['Q'] )#,
#                        p0 = [0.025, 1E4, 8/3., 0., 1.],
#                        bounds=[[1E-5, 0., 5/3., -500., 1E-3],
#                                [1., 1E8, 10., 500., 1E2]] )
popt, pcov = curve_fit( _manning, data['Stage'], data['Q'],
                        p0 = [0.025, 1E4, 8/3., 0., 1.],
                        bounds=[[1E-5, 0., 8/3., -500., 1E-3],
                                [1., 1E8, 10., 500., 1E2]] )

_h = np.arange(0.,10.1, 0.1)
plt.plot(data['Stage'], data['Q'], 'k.')
plt.plot(_h, _manning(_h, *popt))
b *= 2
plt.plot(_h, _manning(_h, *popt))


def _pl(h, k, P):
    return k * h**P
    
popt, pcov = curve_fit( _pl, data['Stage'], data['Q'])

def _lin(h, k):
    return k * h

popt, pcov = curve_fit( _lin, data['Stage'], data['Q'])

