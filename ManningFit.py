#! /usr/bin/python3

import argparse
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from functools import partial
import warnings
import sys

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

def makemanning(channelwidth, slope, use_Rh):
    return partial(_manning, channelwidth=channelwidth, slope=slope, use_Rh=use_Rh)


parser = argparse.ArgumentParser(description='stores the name of your data file, the delimiter which separates your data, your channel width, and slope.')

parser.add_argument('-f', '--configfile', type=str, help='configuration YAML file name')
parser.add_argument('-d', '--datafile', type=str, help='file with two columns: Q, stage')
parser.add_argument('--delimiter', type=str, default='\t', help='"tab", "comma", or "semicolon"')
parser.add_argument('-c', '--channel_width', type=float, default=None, help='river-channel width')
parser.add_argument('-H', '--channel_depth', type=float, default=None, help='river-channel depth (not flow depth)')
parser.add_argument('-s', '--slope', type=float, default=None, help='channel slope')
parser.add_argument('--use_depth', action='store_true', default=False, help='Use flow depth instead of hydraulic radius.')
parser.add_argument('--us_units', action='store_true', default=False, help='Convert imported data from cfs and feet')
parser.add_argument('--plot', default=False, action='store_true', help='Plot h-Q relationship')

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)
    
if args.configfile is not None:
    warnings.warn( "Configfile not yet configured. The irony." )
    parser.print_help()
    sys.exit(0)
else:
    try:
        data = pd.read_csv(args.datafile, sep=args.delimiter)
    except:
        parser.print_help()
        sys.exit(0)

if args.delimiter=='tab':
    args.delimiter='\t'
elif args.delimiter=='comma':
    args.delimiter=','
elif args.delimiter=='semicolon':
    args.delimiter=';'

# To metric, if needed
if args.us_units:
    print(data.columns)
    print(args.use_depth)
    data['Q'] /= 3.28**3
    data['Stage'] /= 3.28

# popt = optimization parameters, pcor = covariance matrix
popt, pcov = curve_fit( makemanning(args.channel_width, args.slope, not args.use_depth), data['Stage'], data['Q'] )

flow_params = { "Manning's n": [popt[0]],
                "Overbank flow coefficient": [popt[1]],
                "Overbank flow power-law exponent": [popt[2]]
              }

outparams = pd.DataFrame.from_dict(flow_params)

outparams.to_csv('flow_params_MinnesotaJordan.csv', index=False)

if args.plot:
    _h = np.arange(0.,10.1, 0.1) # Fixed for now
    plt.plot(data['Stage'].to_list(), data['Q'].to_list(), 'k.')
    plt.plot(_h, makemanning(args.channel_width, args.slope, not args.use_depth)(_h, *popt))
    #plt.plot(_h, makemanning(2*args.channelwidth, args.slope)(_h, *popt))
    plt.show()
