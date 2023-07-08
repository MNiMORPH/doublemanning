import argparse
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from functools import partial

def _manning(h, n, k_Qbank, P_Qbank, stage_depth_Q_offset, h_bank, channelwidth: float, slope: float, use_Rh=True):
    """
    Returns discharge given flow depth, 
    * h: Input. Stage.
    * n: Manning's n
    * k_Qbank, P_Qbank: Power-law fititng parameters for floodplain hypsometry,
                        thereby controlling flow depth beyond the channel, that 
                        also fold in the value of Manning's n and slope
                        (k_Qbank) and the 5/3 exponent to h (P_Qbank)
    * stage_depth_Q_offset: Q(h=0), meant to solve for the offset between
                            bed elevation and flow stage
    * h_bank: Streambank elevation. This might be known a priori, but it can
              also be solved here as a function of the inflection in the
              rating-curve data
    * use R_h: Compute hydraulic radius for the channel and use this to compute
               the parameters
    """
    if use_Rh:
        Rh = h*channelwidth/(2*h+channelwidth)
        Q_ch = channelwidth * Rh**(5/3.) * slope**(1/2.) / n
    else:
        Q_ch = channelwidth * h**(5/3.) * slope**(1/2.) / n
    _ob = (h > h_bank)
    Q_fp = _ob * k_Qbank * (h-h_bank)**(P_Qbank * _ob)
    return Q_ch + Q_fp + stage_depth_Q_offset

def makemanning(channelwidth, slope):
    return partial(_manning, channelwidth=channelwidth, slope=slope)


parser = argparse.ArgumentParser(description='stores the name of your data file, the delimiter which separates your data, your channel width, and slope.')

parser.add_argument('filename', type=str, help='specify the name of the file containing your data')
parser.add_argument('delimiter', type=str, help='specify the type of delimiter your data is separated by')
parser.add_argument('-c', '--channelwidth', type=float, default=70, help='specify the width of your channel')
parser.add_argument('-s', '--slope', type=float, default=1E-4, help='specify your slope')
parser.add_argument('-H', '--use_depth', action='store_true', default=True, help='Use flow depth instead of hydraulic radius.')
args = parser.parse_args()
if args.delimiter=='tab':
    args.delimiter='\t'
elif args.delimiter=='comma':
    args.delimiter=','
elif args.delimiter=='semicolon':
    args.delimiter=';'
    
data = pd.read_csv(args.filename, sep=args.delimiter)

# To metric -- because of USA units here
print(data.columns)
data['Q'] /= 3.28**3
data['Stage'] /= 3.28

# popt = optimization parameters, pcor = covariance matrix
popt, pcov = curve_fit( makemanning(args.channelwidth, args.slope), data['Stage'], data['Q'] )

flow_params = { "Manning's n": [popt[0]],
                "Overbank flow coefficient": [popt[1]],
                "Overbank flow power-law exponent": [popt[2]]
              }

outparams = pd.DataFrame.from_dict(flow_params)

outparams.to_csv('flow_params_MinnesotaJordan.csv', index=False)

_h = np.arange(0.,10.1, 0.1)
plt.plot(data['Stage'], data['Q'], 'k.')
plt.plot(_h, makemanning(args.channelwidth, args.slope)(_h, *popt))
plt.plot(_h, makemanning(2*args.channelwidth, args.slope)(_h, *popt))
plt.show()
