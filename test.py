from depth_from_discharge import FlowDepthDoubleManning
import pandas as pd

flow_params = pd.read_csv('flow_params_MinnesotaJordan.csv')

h = FlowDepthDoubleManning()

h.set_n(flow_params["Manning's n"])
h.set_k(flow_params["Overbank flow coefficient"])
h.set_P(flow_params["Overbank flow power-law exponent"])

h.set_h_bank(50)
h.set_b(70)
h.set_S(1E-4)

h.set_Q(100)

#h.flow_depth_from_Manning_discharge(3)

print( h.run() )
print( h.run(200) )
print( h.run(1000) )
