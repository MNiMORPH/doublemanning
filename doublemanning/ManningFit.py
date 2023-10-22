#! /usr/bin/python3

# Written by A. Wickert

import argparse
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit, fsolve
from matplotlib import pyplot as plt
from functools import partial
import warnings
import sys
from sklearn.metrics import mean_squared_error
import yaml


#############
# FUNCTIONS #
#############

def _manning(stage, n, k_Qbank, P_Qbank, stage_depth_h_offset, h_bank,
              channelwidth: float, slope: float, use_Rh=True):
    """
    Returns discharge given flow depth, 
    * stage: Input. Flow stage.
    * n: Manning's n
    * k_Qbank, P_Qbank: Power-law fititng parameters for floodplain hypsometry,
                        thereby controlling flow depth beyond the channel, that 
                        also fold in the value of Manning's n and slope
                        (k_Qbank) and the 5/3 exponent to h (P_Qbank)
    * stage_depth_h_offset: stage(Q=0); this is the stage at which h=0
    * h_bank: Streambank elevation. This might be known a priori, but it can
              also be solved here as a function of the inflection in the
              rating-curve data
    * use R_h: Compute hydraulic radius for the channel and use this to compute
               the parameters
    """
    h = stage - stage_depth_h_offset
    if use_Rh:
        Rh = h*channelwidth/(2* np.minimum(h, h_bank) + channelwidth)
        Q_ch = np.sign(Rh) * \
                  channelwidth * np.abs(Rh)**(5/3.) * slope**(1/2.) / n
    else:
        Q_ch = np.sign(h) * \
                  channelwidth * np.abs(h)**(5/3.) * slope**(1/2.) / n
    _ob = (h > h_bank)
    Q_fp = _ob * k_Qbank * (h-h_bank)**(P_Qbank * _ob)
    return Q_ch + Q_fp

def calib_manning(channeldepth, channelwidth, slope, use_Rh):
    return partial( _manning, h_bank=channeldepth, channelwidth=channelwidth,
                    slope=slope, use_Rh=use_Rh )

def calib_manning_depth(channelwidth, slope, use_Rh):
    return partial( _manning, channelwidth=channelwidth, slope=slope,
                    use_Rh=use_Rh )

def calib_manning_depth_width(slope, use_Rh):
    return partial( _manning, slope=slope, use_Rh=use_Rh )

################
# MAIN PROGRAM #
################

def main():

    ##########
    # PARSER #
    ##########

    parser = argparse.ArgumentParser( description=
              'Pass channel and flow characteristics to obtain a '+
              '"Double Manning" -- Manning\'s Equation (channel) '+
              ' + generic power-law (floodplain) stage--discharge '+
              'relationship.'
              )

    parser.add_argument('-y', '--configfile', type=str,
                            help='YAML file from which all inputs are read.')
    parser.add_argument('-f', '--datafile', type=str,
                            help='file with two columns: Discharge, Stage')
    parser.add_argument('--delimiter', type=str, default='\t',
                            help='"tab", "comma", or "semicolon"')
    parser.add_argument('-b', '--channel_width', type=float, default=None,
                            help='river-channel width')
    parser.add_argument('-H', '--channel_depth', type=float, default=None,
                            help='river-channel depth (not flow depth)')
    parser.add_argument('-s', '--slope', type=float, default=None,
                            help='channel slope')
    parser.add_argument('-o', '--outfile', default=None,
                            help='Stores fit parameters.')
    parser.add_argument('--use_depth', action='store_true', default=False,
                            help='Use flow depth instead of hydraulic radius.')
    parser.add_argument('--us_units', action='store_true', default=False,
                            help='Convert imported data from cfs and feet')
    parser.add_argument('--plot', default=False, action='store_true',
                            help='Plot stage-discharge relationship')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                            help='Plot stage-discharge relationship')

    # Parse args if anything is passed.
    # If nothing is passed, then print help and exit.
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    if args.configfile is not None:
        try:
            with open(args.configfile, "r") as yamlfile:
                yconf = yaml.load(yamlfile, Loader=yaml.FullLoader)
        except:
            print("\nCould not read from", args.configfile, "\n")
            sys.exit(2)

        # Data
        datafile = yconf['data']['filename']
        delimiter = yconf['data']['delimiter']
        us_units = yconf['data']['us-units']
        
        # Channel
        # None-type --> Include as free variable rather than specifying
        try:
            channel_width = float(yconf['channel']['width'])
        except:
            channel_width = None
        try:
            channel_depth = float(yconf['channel']['depth'])
        except:
            channel_depth = None
        slope = float(yconf['channel']['slope'])
        use_depth = yconf['channel']['use_depth']
        
        # Output
        outfile = yconf['output']['outfile']
        plotflag = yconf['output']['plot']
        verboseflag = yconf['output']['verbose']
    else:
        # Data
        datafile = args.datafile
        delimiter = args.delimiter
        us_units = args.us_units
        
        # Channel
        channel_width = args.channel_width
        channel_depth = args.channel_depth
        slope = args.slope
        use_depth = args.use_depth
        
        # Output
        outfile = args.outfile
        plotflag = args.plot
        verboseflag = args.verbose
        
        # Indicate that there is no YAML file
        yconf = None


    ###############
    # IMPORT DATA #
    ###############

    # Change delimiter string into its appropriate character
    if delimiter=='tab':
        delimiter='\t'
    elif delimiter=='comma':
        delimiter=','
    elif delimiter=='semicolon':
        delimiter=';'

    # Import data
    try:
        data = pd.read_csv(datafile, sep=delimiter)
    except:
        print("\nCould not read from", datafile, "\n")
        sys.exit(2)    

    # To metric, if needed
    if us_units:
        if verboseflag:
            print( "" )
            print("Data columns :", data.columns)
        data['Discharge'] /= 3.28**3
        data['Stage'] /= 3.28


    ##############################
    # DEPTH OR HYDRAULIC RADIUS? #
    ##############################

    if use_depth:
        if verboseflag:
            print( "" )
            print( "Using flow depth instead of hydraulic radius" )


    ###################################################
    # CURVE FIT: STAGE-DISCHARGE DATA WITH 2X MANNING #
    ###################################################

    # Default bounds
    mannings_n_bounds = (0, 0.2)
    floodplain_coeff_bounds = (0, np.inf)
    floodplain_exponent_bounds = (0, np.inf)
    stage_offset_bounds = (-np.inf, np.inf)
    channel_depth_bounds = (0, np.inf)
    channel_width_bounds = (0, np.inf)
    
    # User-set bounds
    # Possible only via YAML file -- too much to really do with simple
    # CLI parsing
    if yconf is not None:
        try:
            mannings_n_bounds = yconf['bounds']['mannings_n_bounds']
        except:
            pass
        try:
            floodplain_coeff_bounds = \
                yconf['bounds']['floodplain_coeff_bounds']
        except:
            pass
        try:
            floodplain_exponent_bounds = \
                yconf['bounds']['floodplain_exponent_bounds']
        except:
            pass
        try:
            stage_offset_bounds = yconf['bounds']['stage_offset_bounds']
        except:
            pass
        try:
            channel_depth_bounds = yconf['bounds']['channel_depth_bounds']
        except:
            pass
        try:
            channel_width_bounds = yconf['bounds']['channel_width_bounds']
        except:
            pass
    
    # Combine these together
    bounds = [ mannings_n_bounds,
               floodplain_coeff_bounds,
               floodplain_exponent_bounds,
               stage_offset_bounds,
               channel_depth_bounds,
               channel_width_bounds ]
               
    # popt = optimization parameters, pcor = covariance matrix
    if channel_width is not None and channel_depth is not None:
        ncalib = 0 # Number of calibrated geometries: width, depth
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4+ncalib], _bounds[1][:4+ncalib])
        # Create the curve fit
        popt, pcov = curve_fit( calib_manning(             channel_depth,
                                                           channel_width,
                                                           slope,
                                                           not use_depth ),
                                data['Stage'], data['Discharge'],
                                bounds=_bounds )
    elif channel_width is not None:
        ncalib = 1
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4+ncalib], _bounds[1][:4+ncalib])
        # Create the curve fit
        popt, pcov = curve_fit( calib_manning_depth(       channel_width,
                                                           slope,
                                                           not use_depth ),
                                data['Stage'], data['Discharge'],
                                bounds=_bounds )
    elif channel_depth is not None:
        ncalib = 1
        sys.exit("Not set up to calibrate an unknown channel width with "+
                 "a known channel depth.")
    else:
        ncalib = 2
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4+ncalib], _bounds[1][:4+ncalib])
        # Create the curve fit
        popt, pcov = curve_fit( calib_manning_depth_width( slope,
                                                           not use_depth ),
                                data['Stage'], data['Discharge'],
                                bounds=_bounds )


    ################
    # COMPUTE RMSE #
    ################

    if channel_width is not None and channel_depth is not None:
        Q_predicted = calib_manning( channel_depth, channel_width,
                                         slope, not use_depth) \
                                         ( data['Stage'], *popt )
    elif channel_width is not None:
        Q_predicted = calib_manning_depth( channel_width, slope,
                                                not use_depth ) \
                                                ( data['Stage'], *popt)
    elif channel_depth is not None:
        sys.exit("Not set up to calibrate an unknown channel width with a known "+
                 "channel depth.")
    else:
        Q_predicted = calib_manning_depth_width( slope,
                                                     not use_depth ) \
                                                     ( data['Stage'], *popt)
    print( Q_predicted - data['Discharge'] )

    # Maybe add this as a plotting option, eventually
    #plt.hist( Q_predicted - data['Discharge'] )
    #plt.show()

    rmse = mean_squared_error( data['Discharge'], Q_predicted, squared=False)

    if verboseflag:
        print( "Fit RMSE [m^3/s]", ":", rmse )


    ##########################
    # PARAMETER DICTIONARIES #
    ##########################

    flow_param_names = [ "Manning's n",
                         "Overbank flow coefficient",
                         "Overbank flow power-law exponent",
                         "Stage at Q = 0 [m]",
                         "Bank height [m]",
                         "Channel width [m]" ]

    _param_sd = np.diag(pcov)**2
    flow_params = {}
    flow_param_SDs = {}
    for i in range(4+ncalib):
        flow_params[flow_param_names[i]] = [popt[i]]
        flow_param_SDs["SD: "+flow_param_names[i]] = _param_sd[i]
    rmse_dict = { "Fit RMSE [m^3/s]": rmse }

    _param_sd = np.diag(pcov)**2

    if verboseflag:
        print( "" )
        print( "PARAMETER VALUES" )
        for key in flow_params:
            print( key, ":", flow_params[key] )
        print( "" )
        print( "PARAMETER STANDARD DEVIATIONS" )
        for key in flow_param_SDs:
            print( key, ":", flow_param_SDs[key] )

    if outfile is not None:
        outarray = {**flow_params, **flow_param_SDs, **rmse_dict}
        outparams = pd.DataFrame.from_dict(outarray)
        outparams.to_csv(outfile, index=False)

    # Add a trailing blank line if we've been verbose
    if verboseflag:
        print( "" )


    ############
    # PLOTTING #
    ############

    if plotflag:
        print( flow_params['Stage at Q = 0 [m]'] )
        print( flow_params['Bank height [m]'] )
        _h = np.arange(0., 2*(flow_params['Bank height [m]'][0] + flow_params['Stage at Q = 0 [m]'][0]), 0.01) # Fixed for now
        plt.plot(data['Stage'].to_list(), data['Discharge'].to_list(), 'k.')
        _xlim = plt.xlim()
        _ylim = plt.ylim()
        if channel_width is not None and channel_depth is not None:
            plt.plot(_h, calib_manning(channel_depth, channel_width, slope, not use_depth)(_h, *popt))
        elif channel_width is not None:
            plt.plot(_h, calib_manning_depth(channel_width, slope, not use_depth)(_h, *popt))
        elif channel_depth is not None:
            sys.exit("Not set up to calibrate an unknown channel width with a known "+
                     "channel depth.")
        else:
            plt.plot(_h, calib_manning_depth_width(slope, not use_depth)(_h, *popt))
        #plt.plot(_h, makemanning(2*args.channelwidth, slope)(_h, *popt))
        plt.xlabel('Stage [m]')
        plt.ylabel('Discharge [m$^3$ s$^{-1}$]')
        plt.xlim((_xlim[0], _xlim[-1]))
        plt.ylim((0, _ylim[-1]))
        plt.show()


################
# ACCESS POINT #
################

if __name__ == "__main__":
    main()

