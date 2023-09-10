#! /usr/bin/python3

# Written by A. Wickert

import argparse
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from functools import partial
import warnings
import sys

#############
# FUNCTIONS #
#############

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

def calib_manning(channeldepth, channelwidth, slope, use_Rh):
    return partial( _manning, h_bank=channeldepth, channelwidth=channelwidth,
                    slope=slope, use_Rh=use_Rh )

def calib_manning_depth(channelwidth, slope, use_Rh):
    return partial( _manning, channelwidth=channelwidth, slope=slope,
                    use_Rh=use_Rh )

def calib_manning_depth_width(slope, use_Rh):
    return partial( _manning, slope=slope, use_Rh=use_Rh )


##########
# PARSER #
##########

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
parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Plot h-Q relationship')

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
    if args.verbose:
        print( "" )
        print("Data columns :", data.columns)
    data['Q'] /= 3.28**3
    data['Stage'] /= 3.28

if args.use_depth:
    if args.verbose:
        print( "" )
        print( "Using flow depth instead of hydraulic radius" )


# popt = optimization parameters, pcor = covariance matrix
if args.channel_width is not None and args.channel_depth is not None:
    ncalib = 0 # Number of calibrated geometries: width, depth
    popt, pcov = curve_fit( calib_manning(             args.channel_depth,
                                                       args.channel_width,
                                                       args.slope,
                                                       not args.use_depth ),
                            data['Stage'], data['Q'] )
elif args.channel_width is not None:
    ncalib = 1
    popt, pcov = curve_fit( calib_manning_depth(       args.channel_width,
                                                       args.slope,
                                                       not args.use_depth ),
                            data['Stage'], data['Q'] )
elif args.channel_depth is not None:
    ncalib = 1
    sys.exit("Not set up to calibrate an unknown channel width with a known "+
             "channel depth.")
else:
    ncalib = 2
    popt, pcov = curve_fit( calib_manning_depth_width( args.slope,
                                                       not args.use_depth ),
                            data['Stage'], data['Q'] )

flow_param_names = [ "Manning's n",
                     "Overbank flow coefficient",
                     "Overbank flow power-law exponent",
                     "Stage depth Q offset",
                     "Bank height",
                     "Channel width" ]

_param_sd = np.diag(pcov)**2
flow_params = {}
flow_param_SDs = {}
for i in range(4+ncalib):
    flow_params[flow_param_names[i]] = [popt[i]]
    flow_param_SDs[flow_param_names[i]] = _param_sd[i]

_param_sd = np.diag(pcov)**2

if args.verbose:
    print( "" )
    print( "PARAMETER VALUES" )
    for key in flow_params:
        print( key, ":", flow_params[key] )
    print( "" )
    print( "PARAMETER STANDARD DEVIATIONS" )
    for key in flow_params:
        print( key, ":", flow_param_SDs[key] )

#outparams = pd.DataFrame.from_dict(flow_params)

#outparams.to_csv('flow_params_MinnesotaJordan.csv', index=False)

# Add a trailing blank line if we've been verbose
if args.verbose:
    print( "" )

if args.plot:
    _h = np.arange(0.,10.1, 0.1) # Fixed for now
    plt.plot(data['Stage'].to_list(), data['Q'].to_list(), 'k.')
    if args.channel_width is not None and args.channel_depth is not None:
        plt.plot(_h, calib_manning(args.channel_depth, args.channel_width, args.slope, not args.use_depth)(_h, *popt))
    elif args.channel_width is not None:
        plt.plot(_h, calib_manning_depth(args.channel_width, args.slope, not args.use_depth)(_h, *popt))
    elif args.channel_depth is not None:
        sys.exit("Not set up to calibrate an unknown channel width with a known "+
                 "channel depth.")
    else:
        plt.plot(_h, calib_manning_depth_width(args.slope, not args.use_depth)(_h, *popt))
        pass
    #plt.plot(_h, makemanning(2*args.channelwidth, args.slope)(_h, *popt))
    plt.show()

