from scipy.optimize import fsolve

class FlowDepthDoubleManning( object ):
    """
    Use Manning's equation to obtain flow depth from river discharge,
    using a conversion from ManningFit.py outputs
    """
    def __init__(self):
        pass
        
    def set_n(self, _var):
        self.n = _var

    def set_k(self, _var):
        self.k = _var

    def set_P(self, _var):
        self.P = _var

    def set_h_bank(self, _var):
        self.h_bank = _var

    def set_b(self, _var):
        self.b = _var
        
    def set_S(self, _var):
        self.S = _var

    def set_Q(self, _var):
        self.Q = _var

    def flow_depth_from_Manning_discharge( self, h ):
        ob = h > self.h_bank
        return self.b/self.n * h**(5/3.) * self.S**0.5 \
                  + ob*self.k*(h-self.h_bank)**(ob*self.P) - self.Q

    def compute_depth(self, Q=None):
        if Q is not None:
            self.Q = Q
        if Q == 0:
            return 0
        else:
            return fsolve( self.flow_depth_from_Manning_discharge, 1. )[0]
    
    def initialize(self, n, k, P, h_bank, b, S):
        self.set_n(n)
        self.set_k(k)
        self.set_P(P)
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

