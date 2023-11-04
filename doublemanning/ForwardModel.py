from scipy.optimize import fsolve

class ForwardModel( object )
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

    def set_S(self, _input):
        """
        Set channel slope
        """
        self.S = _var
        
    # Set all of the parameters at once from directly passed information
    def set_all_params(self, n, k, P, stage_offset, h_bank, b, S):
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
        """
        params = pd.read_csv(_csv_path)
        self.set_n( params["Manning's n"] )
        self.set_k( params["Floodplain flow coefficient"] )
        self.set_P( params["Floodplain flow power-law exponent"] )
        self.set_stage_offset( params["Stage at Q = 0 [m]"] )
        self.set_h_bank( params["Bank height [m]"] )
        self.set_b( params["Channel width [m]"] )
        self.set_S( params["Channel slope"] )

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

    def initialize(self, paramfile)
        """
        Inspired by the CSDMS BMI, but taking a CSV instead of a YAML
        based on the output from ManningFit
        This could eventually be from a YAML, but I will save that possibility
        for a future.
        
        :param paramfile: path to file with double-Manning fit parameters
        """
        self.set_parameters_from_DoubleManning_fit( paramfile )

    def finalize(self):
        pass

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
    
    def discharge_from_stage( self, stage=None ):
        # flow depth
        self.stage = stage # unnecessary but does no harm
        h = stage - self.stage_offset
        # Does the flow go overbank?
        ob = h > self.h_bank
        if self.use_Rh:
            _r = h*self.b / (2*h + self.b)
        else:
            _r = h
        return self.b/self.n * _r**(5/3.) * self.S**0.5 \
                  + ob*self.k*(h-self.h_bank)**(ob*self.P)


class FlowDepthDoubleManning( ForwardModel ):
    """
    Compute depth from discharge.
    """

    def update(self, Q=None):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.h = self.depth_from_discharge(Q)
        return self.h

    def run(self, Q=None):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update(Q)

class StageDoubleManning( ForwardModel ):
    """
    Compute stage from discharge.
    """

    def update(self, Q=None):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.stage = self.stage_from_discharge(Q)
        return self.stage

    def run(self, Q=None):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update(Q)

class DischargeDoubleManning( ForwardModel ):
    """
    Compute discharge from stage
    """

    def update(self, stage=None):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        """
        self.Q = self.compute_depth(Q)
        self.stage = self.h + self.stage_offset
        return self.stage

    def run(self, stage=None):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        (i.e., BMI)
        Same as "update" step
        """
        return self.update(stage)
    
    
