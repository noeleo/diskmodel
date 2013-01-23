'''
I wrote this little program for the sole purpose of calculating the uncertainty of our flux measurement.
'''

import math
import numpy
import matplotlib.pyplot as plt
import sys
from disk import Disk

#Read Data
if sys.argv[2] == 'ensemble':
    f = open(sys.argv[1],'r')
    print 'Reading in data from', sys.argv[1], 'as an ensemble run.'
    walkers = []
    trials = []
    innerRadsteps = []
    outerRadsteps = []
    grainSizesteps = []
    diskMasssteps = []
    powerLawsteps = []
    grainEfficiencysteps = []
    beltMasssteps =[]
    chisteps = []
    acceptance = []
    line = f.readline()
    while line != '':
        info = line.split()
        trials.append(int(info[0]))
        walkers.append(int(info[1]))
        innerRadsteps.append(float(info[2]))
        outerRadsteps.append(float(info[3]))
        grainSizesteps.append(float(info[4])) #Note that this is still log.
        diskMasssteps.append(float(info[5]))  #This one too.
        powerLawsteps.append(float(info[6]))
        grainEfficiencysteps.append(float(info[7]))
        beltMasssteps.append(float(info[8]))
        chisteps.append(float(info[9]))
        acceptance.append(int(info[10])) #There's no more chi2; it wasn't very interesting.
        line = f.readline()
    f.close()

#Read Data
if sys.argv[2] == 'mh':
    f = open(sys.argv[1],'r')
    print 'Reading in data from', sys.argv[1], 'as a Metropolis-Hastings run.'
    innerRadsteps = []
    outerRadsteps = []
    grainSizesteps = []
    diskMasssteps = []
    powerLawsteps = []
    grainEfficiencysteps = []
    beltMasssteps =[]
    chisteps = []
    acceptance = []
    line = f.readline()
    while line != '':
        info = line.split()
        innerRadsteps.append(float(info[1]))
        outerRadsteps.append(float(info[2]))
        grainSizesteps.append(float(info[3])) #Note that this is still log.
        diskMasssteps.append(float(info[4]))  #This one too.
        powerLawsteps.append(float(info[5]))
        grainEfficiencysteps.append(float(info[6]))
        beltMasssteps.append(float(info[7]))
        chisteps.append(float(info[8]))
        acceptance.append(int(info[10])) #Position 9 contains chi2 (as opposed to what I called chi1).
        line = f.readline()
    f.close()

#Calculate total flux for each model
print 'Calculating list of fluxes...'
Disk = Disk(71.65, 75.23, 0.0452, 8.463e-5, 0.5, 0.387, 1.255e-6) #Those values are just placeholders.
totalflux = []
for i in xrange(len(innerRadsteps)):
    Disk.changeParameters(innerRadsteps[i], outerRadsteps[i], 10**grainSizesteps[i], 10**diskMasssteps[i], powerLawsteps[i], grainEfficiencysteps[i], 10**beltMasssteps[i])
    flux = Disk.calculateFlux(1.3e-3)*1.3e-3/2.998e8 #Because calculateFlux returns Jy*Hz
    #print 'step =', i, 'flux =', flux, 'Jy'
    totalflux.append(flux)
    
chop = int(math.ceil(len(acceptance)*0.15)) #Ignore the first 15% of the chain.    
print 'Total Flux: ', 'Mean =', numpy.average(totalflux[chop:]), 'Median =', numpy.median(totalflux[chop:]), 'STD =', numpy.std(totalflux[chop:])
chibest = min(chisteps)
print 'The best value for chi-squared was', chibest, 'and the median value for the accepted chi-squared was', numpy.median(chisteps)
bestmodels = []
for i in xrange(len(chisteps)): 
    if chisteps[i]==chibest:
        bestmodels.append(i+1) #+1 because chisteps should have 1 fewer entry, namely the first one
print 'This value is seen at step(s):', bestmodels
top = bestmodels[0] 
print 'For that first model,', 'Total Flux =', totalflux[top]

print 'Getting sigma from posterior distribution...'
fluxchop=totalflux[chop:]
deltas=[abs(x-totalflux[top]) for x in fluxchop]
deltas.sort()
print 'Sigma from posterior = ', deltas[int(0.6827*len(deltas))]