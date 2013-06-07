import math
import numpy
import matplotlib.pyplot as plt
import sys
from disk_as import Disk
from matplotlib.patches import Rectangle

'''
This is the presentation part of the MCMC code using data saved into a text file.

Ex:  python mcmcreader.py MCMC_Chains/mcmc_ensemble1113 'ensemble'
'''

if len(sys.argv) < 3:
    print 'Error:  I need more arguments.  Give me something like this:'
    print "python mcmcreader.py MCMC_Chains/mcmc_ensemble1113 'ensemble'"
    exit()

if sys.argv[2] != 'ensemble' and sys.argv[2] != 'mh':
    print 'Error:  Give me a mode as sys.argv[2].  Either "ensemble" or "mh".'
    exit()

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
    print 'Data read!  Beginning to make boxes...'

#Mean and Standard Deviation
chop = int(math.ceil(len(acceptance)*0.10)) #Ignore the first 10% of the chain.
length = len(acceptance)
print 'Number of Steps =', length
print 'Inner Radius: ', 'Mean =', numpy.average(innerRadsteps[chop:]), 'Median =', numpy.median(innerRadsteps[chop:]), 'STD =', numpy.std(innerRadsteps[chop:])
print 'log(Grain Size): ', 'Mean =', numpy.average(grainSizesteps[chop:]), 'Median =', numpy.median(grainSizesteps[chop:]), 'STD =', numpy.std(grainSizesteps[chop:])
print 'log(Disk Mass): ', 'Mean =', numpy.average(diskMasssteps[chop:]), 'Median =', numpy.median(diskMasssteps[chop:]),'STD =', numpy.std(diskMasssteps[chop:])
print 'Grain Efficiency: ', 'Mean =', numpy.average(grainEfficiencysteps[chop:]), 'Median =', numpy.median(grainEfficiencysteps[chop:]),'STD =', numpy.std(grainEfficiencysteps[chop:])
print 'log(Belt Mass): ', 'Mean =', numpy.average(beltMasssteps[chop:]), 'Median =', numpy.median(beltMasssteps[chop:]),'STD =', numpy.std(beltMasssteps[chop:])
print 'Acceptance Ratio =', float(numpy.sum(acceptance))/len(acceptance)
print '-------------------------------------------------'

#Locate Best Fit
chibest = min(chisteps)
print 'The best value for chi-squared was', chibest, 'and the median value for the accepted chi-squared was', numpy.median(chisteps)
bestmodels = []
if sys.argv[2] == 'mh':
    for i in xrange(len(chisteps)): 
        if chisteps[i]==chibest:
            bestmodels.append(i+1) #+1 because chisteps should have 1 fewer entry, namely the first one
    print 'This value is seen at step(s):', bestmodels
if sys.argv[2] == 'ensemble':
    for i in xrange(len(chisteps)): 
        if chisteps[i]==chibest:
            bestmodels.append(i)
    n_walkers = walkers[max(walkers)]+1 #+1 just because I want to count from 1, whereas the mcmc code counts from 0.  Doesn't matter.
    print 'This value is seen at trial(s)', [(x - (x % n_walkers))/n_walkers for x in bestmodels], 'of walker(s)', [x % n_walkers for x in bestmodels]
top = bestmodels[0] 

#Calculate Uncertainties.  

IR_deltas=[abs(x-innerRadsteps[top]) for x in innerRadsteps[chop:]]
IR_deltas.sort()
innerRadsigma = IR_deltas[int(0.6827*len(IR_deltas))]

OR_deltas=[abs(x-outerRadsteps[top]) for x in outerRadsteps[chop:]]
OR_deltas.sort()
outerRadsigma = OR_deltas[int(0.6827*len(OR_deltas))]

GS_deltas=[abs(x-grainSizesteps[top]) for x in grainSizesteps[chop:]]
GS_deltas.sort()
grainSizesigma = GS_deltas[int(0.6827*len(GS_deltas))]

DM_deltas=[abs(x-diskMasssteps[top]) for x in diskMasssteps[chop:]]
DM_deltas.sort()
diskMasssigma = DM_deltas[int(0.6827*len(DM_deltas))]

PL_deltas=[abs(x-powerLawsteps[top]) for x in powerLawsteps[chop:]]
PL_deltas.sort()
powerLawsigma = PL_deltas[int(0.6827*len(PL_deltas))]

GE_deltas=[abs(x-grainEfficiencysteps[top]) for x in grainEfficiencysteps[chop:]]
GE_deltas.sort()
grainEfficiencysigma = GE_deltas[int(0.6827*len(GE_deltas))]

BM_deltas=[abs(x-beltMasssteps[top]) for x in beltMasssteps[chop:]]
BM_deltas.sort()
beltMasssigma = BM_deltas[int(0.6827*len(BM_deltas))]

print 'For that first model...' 
print 'Inner Radius =', innerRadsteps[top], '+/-', innerRadsigma
print 'Outer Radius =', outerRadsteps[top], '+/-', outerRadsigma
print 'log(Grain Size) =', grainSizesteps[top], '+/-', grainSizesigma
print 'log(Disk Mass) =', diskMasssteps[top], '+/-', diskMasssigma
print 'Power Law =', powerLawsteps[top],'+/-', powerLawsigma
print 'Grain Efficiency =', grainEfficiencysteps[top], '+/-', grainEfficiencysigma
print 'log(Belt Mass) =', beltMasssteps[top], '+/-', beltMasssigma

#Calculate Autocorrelation Time for Each parameter
if sys.argv[3] == 'acor':
    def C_f(array,T):
        avg = numpy.average(array)
        M = len(array)
        return sum([(array[T+m]-avg)*(array[m]-avg) for m in range(M-T)])/(M-T)
    
    def tau(array):
        return 1 + 2*sum([C_f(array,T) for T in [i+1 for i in range(len(array)-1)]])/C_f(array,0)
    
    print 'Calculating autocorrelation times.'
    taus = []
    for param in [innerRadsteps,grainSizesteps,diskMasssteps,grainEfficiencysteps,beltMasssteps]:
        taus.append(tau(param))
    print 'Autocorrelation times for R_in, a, M_D, beta, and M_B:', taus
    exit()

#Implement IDL's smooth function
def smooth(array, number):
    remainder = len(array) % number
    smooth_array = [0 for dummy in array]
    for i in range(0,len(array)-remainder,number):
        smooth_avg = numpy.mean(array[i:i+number])
        for j in range(i,i+number):
            smooth_array[j] = smooth_avg
    for i in range(0,remainder):
        smooth_avg = numpy.mean(array[len(array)-remainder:])
        for j in range(len(array)-remainder,len(array)):
            smooth_array[j] = smooth_avg
    return smooth_array
    return array

#Plot the chain

plt.figure(1, figsize=(9,8))

smoothstep = length/100

plt.subplot(321)
plt.plot(range(len(chisteps)), smooth(chisteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'$\chi ^2$', fontsize=16)
#plt.title(r'$\chi ^2$')
#plt.ylim(456090,456130)

plt.subplot(322)
plt.plot(range(len(innerRadsteps)), smooth(innerRadsteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'$R_{in}$ (AU)', fontsize=16)
#plt.title(r'$R_{in}$')

plt.subplot(323)
plt.plot(range(len(grainSizesteps)), smooth(grainSizesteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'log(a) ($\mu m$)', fontsize=16)
#plt.title(r'log(a)')

plt.subplot(324)
plt.plot(range(len(diskMasssteps)), smooth(diskMasssteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'log($M_D$) ($M_{\oplus}$)', fontsize=16)
#plt.title(r'log($M_D$)')

plt.subplot(325)
plt.plot(range(len(grainEfficiencysteps)), smooth(grainEfficiencysteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'$\beta$', fontsize=16)
#plt.title(r'$\beta$')

plt.subplot(326)
plt.plot(range(len(beltMasssteps)), smooth(beltMasssteps, smoothstep))
plt.xlabel('Steps', fontsize=16)
plt.ylabel(r'log($M_B$) ($M_{\oplus}$)', fontsize=16)
#plt.title(r'log($M_B$)')

plt.subplots_adjust(wspace=0.6, hspace=0.4)

#plt.savefig('chain_'+str(sys.argv[1])+'.png')

#Histograms
plt.figure(2, figsize=(9,7))

histobars = 20

#INNER RADIUS
minIR = 55.0
maxIR = 80.0
innerRadsteps2 = [i for i in innerRadsteps[chop:] if i > minIR and i < maxIR]
weighter2ir = [1./length for i in innerRadsteps2]

ax1 = plt.subplot(321)
n1, bins1, patches1 = ax1.hist(innerRadsteps2, histobars, weights=weighter2ir, normed=0, facecolor='green', alpha=0.75)
ax1.set_xlabel(r'$R_{in}$ [AU]', fontsize=16)
ax1.set_ylabel('Fraction', fontsize=16)
ax1.grid(True)
plt.axvline(x=innerRadsteps[top], ymin=0, ymax=100, color='k', linewidth=3)
plt.axvline(x=62.1, ymin=0, ymax=100, color='c', linewidth=3)
plt.axvline(x=60.4, ymin=0, ymax=100, color='c', linewidth=3)

#Tick Customization
ticks1 = numpy.arange(55,81,5)
plt.xticks(ticks1)
plt.yticks(numpy.arange(0,0.21,0.05))

print 'Mode of Inner Radius =', bins1[numpy.where(numpy.array(n1)==max(numpy.array(n1)))], ', where bins are of size', bins1[1]-bins1[0]

#GRAIN SIZE
minGS = 0.1
maxGS = 0.65
grainSizesteps2 = [i for i in grainSizesteps[chop:] if i > minGS and i < maxGS]
weighter2gs = [1./length for i in grainSizesteps2]

ax2 = plt.subplot(322)
n2, bins2, patches2 = ax2.hist(grainSizesteps2, histobars, weights=weighter2gs, normed=0, facecolor='green', alpha=0.75)
ax2.set_xlabel(r'log(a [$\mu m$])', fontsize=16)
ax2.set_ylabel('Fraction', fontsize=16)
ax2.grid(True)
plt.axvline(x=grainSizesteps[top], ymin=0, ymax=100, color='k', linewidth=3)

#Tick Customization
ticks2 = numpy.arange(0.1,0.71,0.1)
plt.xticks(ticks2)
plt.yticks(numpy.arange(0,0.17,0.04))

print 'Mode of Grain Size =', bins2[numpy.where(numpy.array(n2)==max(numpy.array(n2)))], ', where bins are of size', bins2[1]-bins2[0]

#DISK MASS
minDM = -3.4
maxDM = -2.5
diskMasssteps2 = [i for i in diskMasssteps[chop:] if i > minDM and i < maxDM]
weighter2dm = [1./length for i in diskMasssteps2]

ax3 = plt.subplot(323)
n3, bins3, patches3 = ax3.hist(diskMasssteps2, histobars, weights=weighter2dm, normed=0, facecolor='green', alpha=0.75)
ax3.set_xlabel(r'log($M_D$ [$M_{\oplus}$])', fontsize=16)
ax3.set_ylabel('Fraction', fontsize=16)
ax3.grid(True)
plt.axvline(x=diskMasssteps[top], ymin=0, ymax=100, color='k', linewidth=3)

#Tick Customization
ticks3 = numpy.arange(-3.4,-2.4,0.2)
plt.xticks(ticks3)
plt.yticks(numpy.arange(0,0.17,0.04))

print 'Mode of Disk Mass =', bins3[numpy.where(numpy.array(n3)==max(numpy.array(n3)))], ', where bins are of size', bins3[1]-bins3[0]

#GRAIN EFFICIENCY
minbeta = 0.2
maxbeta = 0.65
grainEfficiencysteps2 = [i for i in grainEfficiencysteps[chop:] if i > minbeta and i < maxbeta]
weighter2ge = [1./length for i in grainEfficiencysteps2]

ax4 = plt.subplot(324)
n4, bins4, patches4 = ax4.hist(grainEfficiencysteps2, histobars, weights=weighter2ge, normed=0, facecolor='green', alpha=0.75)
ax4.set_xlabel(r'$\beta$', fontsize=16)
ax4.set_ylabel('Fraction', fontsize=16)
ax4.grid(True)
plt.axvline(x=grainEfficiencysteps[top], ymin=0, ymax=100, color='k', linewidth=3)

#Tick Customization
ticks4 = numpy.arange(0.2,0.71,0.1)
plt.xticks(ticks4)
plt.yticks(numpy.arange(0,0.21,0.05))

print 'Mode of Grain Efficiency =', bins4[numpy.where(numpy.array(n4)==max(numpy.array(n4)))], ', where bins are of size', bins4[1]-bins4[0]

#BELT MASS
minBM = -6.15
maxBM = -5.7
beltMasssteps2 = [i for i in beltMasssteps[chop:] if i > minBM and i < maxBM]
weighter2bm = [1./length for i in beltMasssteps2]

ax5 = plt.subplot(325)
n5, bins5, patches5 = ax5.hist(beltMasssteps2, histobars, weights=weighter2bm, normed=0, facecolor='green', alpha=0.75)
ax5.set_xlabel(r'log($M_B$ [$M_{\oplus}$])', fontsize=16)
ax5.set_ylabel('Fraction', fontsize=16)
ax5.grid(True)
plt.axvline(x=beltMasssteps[top], ymin=0, ymax=100, color='k', linewidth=3)

#Tick Customization
ticks5 = numpy.arange(-6.2,-5.65,0.1)
plt.xticks(ticks5)
plt.yticks(numpy.arange(0,0.17,0.04))

print 'Mode of Belt Mass =', bins5[numpy.where(numpy.array(n5)==max(numpy.array(n5)))], ', where bins are of size', bins5[1]-bins5[0]

plt.subplots_adjust(wspace=0.4, hspace=0.7)

#SED
plt.figure(3)
disk = Disk(innerRadsteps[top], outerRadsteps[top], 10**grainSizesteps[top], 10**diskMasssteps[top], powerLawsteps[top], grainEfficiencysteps[top], 10**beltMasssteps[top])
disk.plotSED()
print 'SED chi-squared =', disk.computeChiSquared()
print 'T_avg =', disk.calculateGrainTemperature(0.5*(innerRadsteps[top] + outerRadsteps[top])*1.496e11)

plt.show()