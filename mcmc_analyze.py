from mcmc_analysis import Ensemble
import sys

filename = sys.argv[1]
chop = int(sys.argv[2])
ensemble = Ensemble(filename, chop)
ensemble.Analyze()
#ensemble.Mode_Chi_Squared()
if sys.argv[3] != 'no_plots':
    ensemble.Make_Plots()
if sys.argv[3] == 'Flux_Calc':
    ensemble.Flux_Calc()