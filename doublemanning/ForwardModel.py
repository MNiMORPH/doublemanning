#! /usr/bin/python3

# Written by A. Wickert

from scipy.optimize import fsolve
import pandas as pd
import argparse
import sys

class ForwardModel( object ):
    """
    Use Manning's equation to obtain either
      * flow depth from river discharge
      or
      * discharge from flow depth
    using a conversion from ManningFit.py outputs
    """
    def __init__(self, use_Rh=True):
        # Default to using hyraulic radius and not just depth
        self.use_Rh = use_Rh
    
    #################################
    # SET DOUBLE-MANNING PARAMETERS #
    #################################

    def set_n(self, _var):
        """
        Set Manning's n
        """
        self.n = _var

    def set_k(self, _var):
        """
        Set floodplain coefficient
        """
        self.k = _var

    def set_P(self, _var):
        """
        Set floodplain power
        """
        self.P = _var
        
    def set_stage_offset(self, _var):
        """
        Set stage offset: depth = stage - stage_offset
        Therefore, this is the stage where Q = 0.
        """
        self.stage_offset = _var

    def set_h_bank(self, _var):
        """
        Set bank height = channel depth
        """
        self.h_bank = _var

    def set_b(self, _var):
        """
        Set channel width
        """
        self.b = _var

    def set_S(self, _var):
        """
        Set channel slope
        """
        self.S = _var
        
    def set_use_Rh(self, _var):
        """
        Set the flag to use hydraulic radius instead of flow depth.
        This should typically be true.
        """
        self.use_Rh = _var
        
    # Set all of the parameters at once from directly passed information
    def set_all_params(self, n, k, P, stage_offset, h_bank, b, S, use_Rh):
        """
        If you wish to set all of the parameters without an input file
        """
        self.set_n(n)
        self.set_k(k)
        self.set_P(P)
        self.set_stage_offset(stage_offset)
        self.set_h_bank(h_bank)
        self.set_b(b)
        self.set_S(S)
        self.set_use_Rh(use_Rh)
    
    # Set all of the parameters from the output file from a double-Manning fit
    def set_parameters_from_DoubleManning_fit(self, _csv_path):
        """
        Reads all parameters from a CSV file provided as output from the
        double-Manning calculation.
        These include:
          * Manning's n
          * Floodplain k
          * Floodplain P
          * Stage offset (stage at Q=0)
          * Bank heights
          * Channel width
          * Channel slope
          * Flag to compute using depth instead of hydraulic radius
        """
        try:
            params = pd.read_csv( _csv_path )
        except:
            print("\nCould not read params from file:\n", _csv_path, "\n")
            sys.exit(2)
        self.set_n( params["Manning's n"] )
        self.set_k( params["Floodplain flow coefficient"] )
        self.set_P( params["Floodplain flow power-law exponent"] )
        self.set_stage_offset( params["Stage at Q = 0 [m]"] )
        self.set_h_bank( params["Bank height [m]"] )
        self.set_b( params["Channel width [m]"] )
        self.set_S( params["Channel slope"] )
        self.set_use_Rh( not params["Use flow depth instead of Rh"].bool() )

    def set_Q(self, _var):
        """
        Set river discharge, if this is to be an input (and depth or stage
        the output).
        This can be a scalar or an array.
        """
        self.Q = _var

    def set_stage(self, _var):
        """
        Set river stage, if this is to be an input (and discharge the output).
        This can be a scalar or an array.
        """
        self.Q = _var

    def set_depth(self, _var):
        """
        Set flow depth, if this is to be an input (and discharge the output).
        This differs from setting stage in assumning that stage_offset = 0.
        This can be a scalar or an array.
        """
        self.Q = _var

    #######################
    # BMI: SHARED METHODS #
    #######################

    def initialize(self, paramfile):
        """
        Inspired by the CSDMS BMI, but taking a CSV instead of a YAML
        based on the output from ManningFit
        This could eventually be from a YAML, but I will save that possibility
        for a future.
        
        :param paramfile: path to file with double-Manning fit parameters
        """
        self.set_parameters_from_DoubleManning_fit( paramfile )

    #####################################
    # Compute flow depth from discharge #
    #####################################

    def _stage_from_discharge_rootfinder( self, stage ):
        """
        Returns a function whose root gives the river stage at the discharge
        set by the class variable self.Q
        """
        # flow depth
        h = stage - self.stage_offset
        # Does the flow go overbank?
        ob = h > self.h_bank
        if self.use_Rh:
            _r = h*self.b / (2*h + self.b)
        else:
            _r = h
        return self.b/self.n * _r**(5/3.) * self.S**0.5 \
                  + ob*self.k*(h-self.h_bank)**(ob*self.P) - self.Q

    def depth_from_discharge(self, Q=None):
        if Q is not None:
            self.Q = Q
        if Q == 0:
            return 0
        else:
            return fsolve( self._stage_from_discharge_rootfinder, 1. )[0]

    ################################
    # Compute stage from discharge #
    ################################

    def stage_from_discharge(self, Q=None):
        return self.depth_from_discharge(Q) + self.stage_offset

    ################################
    # Compute discharge from stage #
    ################################
    
    def discharge_from_stage( self, _stage=None ):
        # flow depth
        if self.stage is None:
            if _stage is None:
                raise Exception("Stage must be set.")
            else: self.stage = _stage # unnecessary but does no harm
        h = self.stage - self.stage_offset
        # Does the flow go overbank?
        ob = h > self.h_bank
        if self.use_Rh:
            _r = h*self.b / (2*h + self.b)
        else:
            _r = h
        self.Q =  self.b/self.n * _r**(5/3.) * self.S**0.5 \
                    + ob*self.k*(h-self.h_bank)**(ob*self.P)
        return self.Q


class FlowDepthDoubleManning( ForwardModel ):
    """
    Compute depth from discharge.
    """
    def __init__(self, Q):
        self.Q = Q

    def update(self):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.h = self.depth_from_discharge(Q)
        return self.h

    def run(self):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update(Q)

    def finalize(self):
        print( self.h.iloc[0] )

class StageDoubleManning( ForwardModel ):
    """
    Compute stage from discharge.
    """
    def __init__(self, Q):
        self.Q = Q

    def update(self):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.stage = self.stage_from_discharge(Q)
        return self.stage

    def run(self):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update(Q)

    def finalize(self):
        print( self.stage.iloc[0] )

class DischargeDoubleManning( ForwardModel ):
    """
    Compute discharge from stage
    """
    def __init__(self, stage):
        self.stage = stage

    def update(self):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.Q = self.discharge_from_stage()
        return self.Q

    def run(self):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update()
    
    def finalize(self):
        print( self.Q.iloc[0] )

################
# MAIN PROGRAM #
################

def main():

    ##########
    # PARSER #
    ##########

    parser = argparse.ArgumentParser( description=
              'Return stage or discharge based on a double-Manning fit. '+
              'All values are SI (mks).'
              )

    parser.add_argument('-p', '--paramfile', type=str,
                            help='CSV file for double-Manning parameters.')
    parser.add_argument('-s', '--stage', type=float, default=None,
                            help='Calculate discharge from this stage.')
    parser.add_argument('-H', '--depth', type=float, default=None,
                            help='Calculate discharge from this flow depth.')
    parser.add_argument('-Q', '--discharge', type=float, default=None,
                            help='Calculate stage from this discharge.')

    # Parse args if anything is passed.
    # If nothing is passed, then print help and exit.
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    if sum( ( args.stage is not None,
              args.depth is not None,
              args.discharge is not None) ) > 1:
        print("\nSelect only one of s, h, Q.\n")
        sys.exit(2)

    if args.stage:
        m2 = DischargeDoubleManning( args.stage )
    elif args.depth:
        m2 = DischargeDoubleManning( args.depth )
    elif args.discharge:
        m2 = StageDoubleManning( args.discharge )

    m2.initialize( args.paramfile )
    m2.run()
    m2.finalize()
    
################
# ACCESS POINT #
################

if __name__ == "__main__":
    main()

