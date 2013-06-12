from mcmc_analysis import Ensemble
import sys

filename = sys.argv[1]
chop = int(sys.argv[2])
ensemble = Ensemble(filename, chop)
ensemble.Analyze()
print ensemble.Mode_Chi_Squared()
if len(sys.argv) == 4:
    if sys.argv[3] == 'plot':
        ensemble.Make_Plots()
    if sys.argv[3] == 'Flux_Calc':
        ensemble.Flux_Calc()
    if sys.argv[3] == 'R_eff':
        ensemble.Get_Rad_eff()