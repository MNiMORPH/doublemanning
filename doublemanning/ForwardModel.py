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

class FlowDepthDoubleManning( object ):
    def flow_depth_from_Manning_discharge( self, stage ):
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

    def compute_depth(self, Q=None):
        if Q is not None:
            self.Q = Q
        if Q == 0:
            return 0
        else:
            return fsolve( self.flow_depth_from_Manning_discharge, 1. )[0]

    def initialize(self, n, k, P, stage_offset, h_bank, b, S):
        self.set_n(n)
        self.set_k(k)
        self.set_P(P)
        self.set_stage_offset(stage_offset)
        self.set_h_bank(h_bank)
        self.set_b(b)
        self.set_S(S)

    def update(self, Q=None):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        """
        self.h = self.compute_depth(Q)
        return self.h

    def run(self, Q=None):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        Same as "update" step
        """
        return self.update(Q)

    def finalize(self):
        pass

