from mcmc_analysis import Ensemble
import sys

filename = sys.argv[1]
chop = int(sys.argv[2])
ensemble = Ensemble(filename, chop)
ensemble.Analyze()