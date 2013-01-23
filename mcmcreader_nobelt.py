import math
import numpy
import matplotlib.pyplot as plt
import sys
from disk import Disk

'''
This is the presentation part of the MCMC code using data saved into a text file.
'''

#Read Data
f = open(sys.argv[1],'r')
innerRadsteps = []
outerRadsteps = []
grainSizesteps = []
diskMasssteps = []
powerLawsteps = []
grainEfficiencysteps = []
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
    chisteps.append(float(info[7]))
    acceptance.append(float(info[9])) #Position 8 contains chi2.
    line = f.readline()
f.close()

#Mean and Standard Deviation
chop = int(math.ceil(len(acceptance)*0.5)) #Ignore the first 50% of the chain.
print 'Number of Steps =', len(acceptance)
print 'Inner Radius: ', 'Mean =', numpy.average(innerRadsteps[chop:]), 'Median =', numpy.median(innerRadsteps[chop:]), 'STD =', numpy.std(innerRadsteps[chop:])
print 'log(Grain Size): ', 'Mean =', numpy.average(grainSizesteps[chop:]), 'Median =', numpy.median(grainSizesteps[chop:]), 'STD =', numpy.std(grainSizesteps[chop:])
print 'log(Disk Mass): ', 'Mean =', numpy.average(diskMasssteps[chop:]), 'Median =', numpy.median(diskMasssteps[chop:]),'STD =', numpy.std(diskMasssteps[chop:])
print 'Grain Efficiency: ', 'Mean =', numpy.average(grainEfficiencysteps[chop:]), 'Median =', numpy.median(grainEfficiencysteps[chop:]),'STD =', numpy.std(grainEfficiencysteps[chop:])
print 'Acceptance Rate =', float(numpy.sum(acceptance))/len(acceptance)

#Locate Best Fit
chibest = min(chisteps)
print 'The best value for chi-squared was', chibest, 'and the median value for the accepted chi-squared was', numpy.median(chisteps)
bestmodels = []
for i in xrange(len(chisteps)): 
    if chisteps[i]==chibest:
        bestmodels.append(i+1)
print 'This value is seen at step(s):', bestmodels
firstmodel = bestmodels[0]
if acceptance[firstmodel] == 0:
    top = firstmodel
elif acceptance[firstmodel] == 1:
    top = firstmodel - 1
print 'For that first model,', 'Inner Radius =', innerRadsteps[top], 'Outer Radius =', outerRadsteps[top], 'Grain Size =', 10**grainSizesteps[top], 'Disk Mass =', 10**diskMasssteps[top], 'Power Law =', powerLawsteps[top], 'Grain Efficiency =', grainEfficiencysteps[top]

#Plot the chain
plt.figure(1)

plt.subplot(321)
plt.plot(range(len(chisteps)), chisteps)
plt.xlabel('Steps')
plt.ylabel(r'$\chi ^2$')
plt.title(r'$\chi ^2$')

plt.subplot(323)
plt.plot(range(len(innerRadsteps)), innerRadsteps)
plt.xlabel('Steps')
plt.ylabel(r'Inner Radius (AU)')
plt.title(r'Inner Radius')

plt.subplot(324)
plt.plot(range(len(grainSizesteps)), grainSizesteps)
plt.xlabel('Steps')
plt.ylabel(r'log(Grain Size) (microns)')
plt.title(r'log(Grain Size)')

plt.subplot(325)
plt.plot(range(len(diskMasssteps)), diskMasssteps)
plt.xlabel('Steps')
plt.ylabel(r'log(Disk Mass) (Earth Masses)')
plt.title(r'log(Disk Mass)')

plt.subplot(326)
plt.plot(range(len(grainEfficiencysteps)), grainEfficiencysteps)
plt.xlabel('Steps')
plt.ylabel(r'Grain Efficiency')
plt.title(r'Grain Efficiency')

plt.subplots_adjust(wspace=0.5, hspace=0.5)
#plt.savefig('chain.ps')

#Histograms
plt.figure(2)

ax1 = plt.subplot(221)
n1, bins1, patches1 = ax1.hist(innerRadsteps[chop:], 50, normed=1, facecolor='green', alpha=0.75)
ax1.set_xlabel('Inner Radius (AU)')
ax1.set_ylabel('Probability')
ax1.set_title('Inner Radius')
ax1.grid(True)

ax2 = plt.subplot(222)
n2, bins2, patches2 = ax2.hist(grainSizesteps[chop:], 50, normed=1, facecolor='green', alpha=0.75)
ax2.set_xlabel('log(Grain Size) (microns)')
ax2.set_ylabel('Probability')
ax2.set_title('log(Grain Size)')
ax2.grid(True)

ax3 = plt.subplot(223)
n3, bins3, patches3 = ax3.hist(diskMasssteps[chop:], 50, normed=1, facecolor='green', alpha=0.75)
ax3.set_xlabel('log(Disk Mass) (Earth Masses)')
ax3.set_ylabel('Probability')
ax3.set_title('log(Disk Mass)')
ax3.grid(True)

ax4 = plt.subplot(224)
n4, bins4, patches4 = ax4.hist(grainEfficiencysteps[chop:], 50, normed=1, facecolor='green', alpha=0.75)
ax4.set_xlabel('Grain Efficiency')
ax4.set_ylabel('Probability')
ax4.set_title('Grain Efficiency')
ax4.grid(True)

plt.subplots_adjust(wspace=0.5, hspace=0.5)

#plt.savefig('histo.ps')

plt.figure(3)
disk = Disk(innerRadsteps[top], outerRadsteps[top], 10**grainSizesteps[top], 10**diskMasssteps[top], powerLawsteps[top], grainEfficiencysteps[top])
disk.plotSED()
print 'SED chi-squared =', disk.computeChiSquared()

plt.show()