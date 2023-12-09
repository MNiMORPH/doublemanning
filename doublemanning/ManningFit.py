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
        Rh = h * channelwidth / (2 * np.minimum(h, h_bank) + channelwidth)
        Q_ch = np.sign(Rh) * \
            channelwidth * np.abs(Rh)**(5 / 3.) * slope**(1 / 2.) / n
    else:
        Q_ch = np.sign(h) * \
            channelwidth * np.abs(h)**(5 / 3.) * slope**(1 / 2.) / n
    _ob = (h > h_bank)
    Q_fp = _ob * k_Qbank * (h - h_bank)**(P_Qbank * _ob)
    return Q_ch + Q_fp


def calib_manning(channeldepth, channelwidth, slope, use_Rh):
    return partial(_manning, h_bank=channeldepth, channelwidth=channelwidth,
                   slope=slope, use_Rh=use_Rh)


def calib_manning_depth(channelwidth, slope, use_Rh):
    return partial(_manning, channelwidth=channelwidth, slope=slope,
                   use_Rh=use_Rh)


def calib_manning_depth_width(slope, use_Rh):
    return partial(_manning, slope=slope, use_Rh=use_Rh)

################
# MAIN PROGRAM #
################


def main():

    ##########
    # PARSER #
    ##########

    parser = argparse.ArgumentParser(
        description='Pass channel and flow characteristics to obtain a ' +
        '"Double Manning" -- Manning\'s Equation (channel) ' +
        ' + generic power-law (floodplain) stage--discharge ' +
        'relationship.')

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

    # Parse configuration file.
    # Assign variables except those for optimization bounds and plotting.
    if args.configfile is not None:
        try:
            with open(args.configfile, "r") as yamlfile:
                yconf = yaml.load(yamlfile, Loader=yaml.FullLoader)
        except BaseException:
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
        except BaseException:
            channel_width = None
        try:
            channel_depth = float(yconf['channel']['depth'])
        except BaseException:
            channel_depth = None
        slope = float(yconf['channel']['slope'])
        use_depth = yconf['channel']['use_depth']

        # Output
        outfile = yconf['output']['outfile']
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

        # Plotting
        display_plot_flag = args.plot
        plot_save_path = None  # No option here

        # Output
        outfile = args.outfile
        verboseflag = args.verbose

        # Indicate that there is no YAML file
        yconf = None

    ###############
    # IMPORT DATA #
    ###############

    # Change delimiter string into its appropriate character
    if delimiter == 'tab':
        delimiter = '\t'
    elif delimiter == 'comma':
        delimiter = ','
    elif delimiter == 'semicolon':
        delimiter = ';'

    # Import data
    try:
        data = pd.read_csv(datafile, sep=delimiter)
    except BaseException:
        print("\nCould not read from", datafile, "\n")
        sys.exit(2)

    # To metric, if needed
    if us_units:
        if verboseflag:
            print("")
            print("Data columns :", data.columns)
        data['Discharge'] /= 3.28**3
        data['Stage'] /= 3.28

    ##############################
    # DEPTH OR HYDRAULIC RADIUS? #
    ##############################

    if use_depth:
        if verboseflag:
            print("")
            print("Using flow depth instead of hydraulic radius")

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
        except BaseException:
            pass
        try:
            floodplain_coeff_bounds = \
                yconf['bounds']['floodplain_coeff_bounds']
        except BaseException:
            pass
        try:
            floodplain_exponent_bounds = \
                yconf['bounds']['floodplain_exponent_bounds']
        except BaseException:
            pass
        try:
            stage_offset_bounds = yconf['bounds']['stage_offset_bounds']
        except BaseException:
            pass
        try:
            channel_depth_bounds = yconf['bounds']['channel_depth_bounds']
        except BaseException:
            pass
        try:
            channel_width_bounds = yconf['bounds']['channel_width_bounds']
        except BaseException:
            pass

    # Combine these together
    bounds = [mannings_n_bounds,
              floodplain_coeff_bounds,
              floodplain_exponent_bounds,
              stage_offset_bounds,
              channel_depth_bounds,
              channel_width_bounds]

    # popt = optimization parameters, pcor = covariance matrix
    if channel_width is not None and channel_depth is not None:
        ncalib = 0  # Number of calibrated geometries: width, depth
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4 + ncalib], _bounds[1][:4 + ncalib])
        # Create the curve fit
        popt, pcov = curve_fit(calib_manning(channel_depth,
                                             channel_width,
                                             slope,
                                             not use_depth),
                               data['Stage'], data['Discharge'],
                               bounds=_bounds)
    elif channel_width is not None:
        ncalib = 1
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4 + ncalib], _bounds[1][:4 + ncalib])
        # Create the curve fit
        popt, pcov = curve_fit(calib_manning_depth(channel_width,
                                                   slope,
                                                   not use_depth),
                               data['Stage'], data['Discharge'],
                               bounds=_bounds)
    elif channel_depth is not None:
        ncalib = 1
        sys.exit("Not set up to calibrate an unknown channel width with " +
                 "a known channel depth.")
    else:
        ncalib = 2
        # Bounds for curve fit
        _bounds = np.array(bounds).transpose()
        _bounds = (_bounds[0][:4 + ncalib], _bounds[1][:4 + ncalib])
        # Create the curve fit
        popt, pcov = curve_fit(calib_manning_depth_width(slope,
                                                         not use_depth),
                               data['Stage'], data['Discharge'],
                               bounds=_bounds)

    ################
    # COMPUTE RMSE #
    ################

    if channel_width is not None and channel_depth is not None:
        Q_predicted = calib_manning(channel_depth, channel_width,
                                    slope, not use_depth)(data['Stage'], *popt)
    elif channel_width is not None:
        Q_predicted = calib_manning_depth(channel_width, slope,
                                          not use_depth)(data['Stage'], *popt)
    elif channel_depth is not None:
        sys.exit("Not set up to calibrate an unknown channel width with a " +
                 "known channel depth.")
    else:
        Q_predicted = calib_manning_depth_width(
            slope, not use_depth)(
            data['Stage'], *popt)

    #print( Q_predicted - data['Discharge'] )
    # Maybe add this as a plotting option, eventually
    #plt.hist( Q_predicted - data['Discharge'] )
    # plt.show()

    rmse = mean_squared_error(data['Discharge'], Q_predicted, squared=False)

    if verboseflag:
        print("Fit RMSE [m^3/s]", ":", rmse)

    ##########################
    # PARAMETER DICTIONARIES #
    ##########################

    flow_param_names = ["Manning's n",
                        "Floodplain discharge coefficient",
                        "Floodplain discharge exponent",
                        "Stage at Q = 0 [m]",
                        "Bank height [m]",
                        "Channel width [m]",
                        "Channel slope"]

    _param_sd = np.diag(pcov)**2
    flow_params = {}
    flow_param_SDs = {}
    for i in range(len(flow_param_names)):
        # First see if it was solved for
        try:
            flow_params[flow_param_names[i]] = [popt[i]][0]
        # If not, see if it was a default-passed value
        # (This should always be the case for slope)
        except BaseException:
            if flow_param_names[i] == "Bank height [m]":
                try:
                    flow_params[flow_param_names[i]] = channel_depth[0]
                except BaseException:
                    flow_params[flow_param_names[i]] = channel_depth
            elif flow_param_names[i] == "Channel width [m]":
                try:
                    flow_params[flow_param_names[i]] = channel_width[0]
                except BaseException:
                    flow_params[flow_param_names[i]] = channel_width
            elif flow_param_names[i] == "Channel slope":
                flow_params[flow_param_names[i]] = slope
            else:
                raise Exception("No candidate for output.")

        # First see if it was solved for
        try:
            flow_param_SDs["SD: " + flow_param_names[i]] = _param_sd[i]
        # If not, see if it was a default-passed value
        # (This should always be the case for slope)
        except BaseException:
            if flow_param_names[i] == "Bank height [m]":
                flow_param_SDs["SD: " + flow_param_names[i]] = 0
            elif flow_param_names[i] == "Channel width [m]":
                flow_param_SDs["SD: " + flow_param_names[i]] = 0
            elif flow_param_names[i] == "Channel slope":
                flow_param_SDs["SD: " + flow_param_names[i]] = 0
            else:
                raise Exception("No candidate for output.")

    rmse_dict = {"Fit RMSE [m^3/s]": rmse}

    _param_sd = np.diag(pcov)**2

    if verboseflag:
        print("")
        print("PARAMETER VALUES")
        for key in flow_params:
            print(key, ":", flow_params[key])
        print("")
        print("PARAMETER STANDARD DEVIATIONS")
        for key in flow_param_SDs:
            print(key, ":", flow_param_SDs[key])

    use_depth_dict = {"Use flow depth instead of Rh": use_depth}

    if outfile is not None:
        outdict = {**flow_params,
                   **flow_param_SDs,
                   **rmse_dict,
                   **use_depth_dict}
        outparams = pd.DataFrame(outdict, index=[0])
        outparams.to_csv(outfile, index=False)

    # Add a trailing blank line if we've been verbose
    if verboseflag:
        print("")

    ############
    # PLOTTING #
    ############

    # Parse plotting options
    # Possible only via YAML file -- too much to really do with simple
    # CLI parsing
    if yconf is not None:
        display_plot_flag = yconf['plotting']['show']
        display_negative_rating_curve = yconf['plotting'] \
                                             ['display_negative_rating_curve']
        try:
            plot_save_path = yconf['plotting']['savepath']
        except BaseException:
            plot_save_path = None
        try:
            plot_xlim_stage_min = float(yconf['plotting']['stage_min'])
        except BaseException:
            plot_xlim_stage_min = None
        try:
            plot_xlim_stage_max = float(yconf['plotting']['stage_max'])
        except BaseException:
            plot_xlim_stage_max = None
        try:
            plot_ylim_discharge_min = float(yconf['plotting']['discharge_min'])
        except BaseException:
            plot_ylim_discharge_min = None
        try:
            plot_ylim_discharge_max = float(yconf['plotting']['discharge_max'])
        except BaseException:
            plot_ylim_discharge_max = None
        # Flags for plot annotations
        plot_stage_offset_hash_bottom = \
            yconf['plotting']['stage_offset_hash_bottom']
        plot_stage_offset_hash_top = \
            yconf['plotting']['stage_offset_hash_top']
        plot_stage_offset_dotted_line = \
            yconf['plotting']['stage_offset_dotted_line']
        plot_bank_height_hash_bottom = \
            yconf['plotting']['bank_height_hash_bottom']
        plot_bank_height_hash_top = \
            yconf['plotting']['bank_height_hash_top']
        plot_bank_height_dotted_line = \
            yconf['plotting']['bank_height_dotted_line']

    # Otherwise, set the yaml-generated variables to None
    # to ensure that the code can still run
    else:
        plot_save_path = None
        plot_xlim_stage_min = None
        plot_xlim_stage_max = None
        plot_ylim_discharge_min = None
        plot_ylim_discharge_max = None
        plot_stage_offset_hash_bottom = True
        plot_stage_offset_hash_top = False
        plot_stage_offset_dotted_line = False
        plot_bank_height_hash_bottom = False
        plot_bank_height_hash_top = False
        plot_bank_height_dotted_line = False

    # Actually plotting now.

    if display_plot_flag or plot_save_path is not None:

        # First, set up the data on the plot and use this to find the
        # suggested axis limits for the full data range
        plt.figure(figsize=(4.5, 3.5))
        plt.plot(data['Stage'].to_list(), data['Discharge'].to_list(), 'k.')

        _xlim = list(plt.xlim())
        _ylim = list(plt.ylim())

        if plot_xlim_stage_min:
            _xlim[0] = plot_xlim_stage_min
        if plot_xlim_stage_max:
            _xlim[-1] = plot_xlim_stage_max
        if plot_ylim_discharge_min:
            _ylim[0] = plot_ylim_discharge_min
        if plot_ylim_discharge_max:
            _ylim[-1] = plot_ylim_discharge_max

        # Next, process stage and discharge to create a reasonable
        # array for plotting the data

        # Create an array of stages for plotting
        _zs = np.linspace(_xlim[0], _xlim[-1], 200)

        # Obtain _Q
        if channel_width is not None and channel_depth is not None:
            _Q = calib_manning(channel_depth,
                               channel_width,
                               slope,
                               not use_depth)(_zs, *popt)
        elif channel_width is not None:
            _Q = calib_manning_depth(channel_width,
                                     slope,
                                     not use_depth)(_zs, *popt)
        elif channel_depth is not None:
            sys.exit("Not set up to calibrate an unknown channel width using "+
                     "a known channel depth.")
        else:
            _Q = calib_manning_depth_width(slope, not use_depth)(_zs, *popt)

        # Do we show negative (nonphysical) rating-curve outcomes?
        # Could also do this using stage_offset, but using Q feels intuitive
        if display_negative_rating_curve is False:
            _zs = _zs[_Q > 0]
            _Q = _Q[_Q > 0]

        # Plot the rating curve
        plt.plot(_zs, _Q, linewidth=4, color='.3', alpha=.7)

        # And set up the plot's display options
        plt.xlabel('Stage [m]', fontsize=14)
        plt.ylabel('Discharge [m$^3$ s$^{-1}$]', fontsize=14)
        plt.xlim((_xlim[0], _xlim[-1]))
        plt.ylim((_ylim[0], _ylim[-1]))

        # Plot hashes or lines where Q=0 (depth = 0) and at bankfull
        if plot_stage_offset_hash_bottom:
            _xlocs = [flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[0]
            _ymax = (_ylim[-1] - _ylim[0]) / 50. + _ylim[0]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=3,
                       zorder=-10)
        if plot_stage_offset_hash_top:
            _xlocs = [flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[-1] - (_ylim[-1] - _ylim[0]) / 50.
            _ymax = _ylim[-1]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=3,
                       zorder=-10)
        if plot_stage_offset_dotted_line:
            _xlocs = [flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[0]
            _ymax = _ylim[-1]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=1,
                       linestyles='dotted', zorder=-10)
        if plot_bank_height_hash_bottom:
            _xlocs = [flow_params['Bank height [m]'] +
                      flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[0]
            _ymax = (_ylim[-1] - _ylim[0]) / 50. + _ylim[0]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=3,
                       zorder=-10)
        if plot_bank_height_hash_top:
            _xlocs = [flow_params['Bank height [m]'] +
                      flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[-1] - (_ylim[-1] - _ylim[0]) / 50.
            _ymax = _ylim[-1]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=3,
                       zorder=-10)
        if plot_bank_height_dotted_line:
            _xlocs = [flow_params['Bank height [m]'] +
                      flow_params['Stage at Q = 0 [m]']]
            _ymin = _ylim[0]
            _ymax = _ylim[-1]
            plt.vlines(_xlocs, _ymin, _ymax, color='.8', linewidth=1,
                       linestyles='dotted', zorder=-10)

        plt.tight_layout()

        if plot_save_path is not None:
            plt.savefig(plot_save_path)
        if display_plot_flag:
            plt.show()


################
# ACCESS POINT #
################
if __name__ == "__main__":
    main()
