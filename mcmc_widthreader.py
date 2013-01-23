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
widthsteps = []
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
    widthsteps.append(float(info[2]))
    grainSizesteps.append(float(info[3])) #Note that this is still log.
    diskMasssteps.append(float(info[4]))  #This one too.
    powerLawsteps.append(float(info[5]))
    grainEfficiencysteps.append(float(info[6]))
    beltMasssteps.append(float(info[7]))
    chisteps.append(float(info[8]))
    acceptance.append(float(info[10])) #Position 9 contains chi2.
    line = f.readline()
f.close()

#Mean and Standard Deviation
chop = int(math.ceil(len(acceptance)*0.15)) #Ignore the first 15% of the chain.
length = len(acceptance)
print 'Number of Steps =', length
print 'Inner Radius: ', 'Mean =', numpy.average(innerRadsteps[chop:]), 'Median =', numpy.median(innerRadsteps[chop:]), 'STD =', numpy.std(innerRadsteps[chop:])
print 'Width: ', 'Mean =', numpy.average(widthsteps[chop:]), 'Median =', numpy.median(widthsteps[chop:]), 'STD =', numpy.std(widthsteps[chop:])
print 'log(Grain Size): ', 'Mean =', numpy.average(grainSizesteps[chop:]), 'Median =', numpy.median(grainSizesteps[chop:]), 'STD =', numpy.std(grainSizesteps[chop:])
print 'log(Disk Mass): ', 'Mean =', numpy.average(diskMasssteps[chop:]), 'Median =', numpy.median(diskMasssteps[chop:]),'STD =', numpy.std(diskMasssteps[chop:])
print 'Grain Efficiency: ', 'Mean =', numpy.average(grainEfficiencysteps[chop:]), 'Median =', numpy.median(grainEfficiencysteps[chop:]),'STD =', numpy.std(grainEfficiencysteps[chop:])
print 'log(Belt Mass): ', 'Mean =', numpy.average(beltMasssteps[chop:]), 'Median =', numpy.median(beltMasssteps[chop:]),'STD =', numpy.std(beltMasssteps[chop:])
print 'Acceptance Rate =', float(numpy.sum(acceptance))/len(acceptance)
print '-------------------------------------------------'

#Locate Best Fit
chibest = min(chisteps)
print 'The best value for chi-squared was', chibest, 'and the median value for the accepted chi-squared was', numpy.median(chisteps)
bestmodels = []
for i in xrange(len(chisteps)): 
    if chisteps[i]==chibest:
        bestmodels.append(i+1)  #+1 because chisteps should have 1 fewer entry, namely the first one
print 'This value is seen at step(s):', bestmodels
top = bestmodels[0] 
print 'For that first model,', 'Inner Radius =', innerRadsteps[top], 'Width =', widthsteps[top],\
 'Grain Size =', 10**grainSizesteps[top], 'Disk Mass =', 10**diskMasssteps[top], 'Power Law =', powerLawsteps[top],\
  'Grain Efficiency =', grainEfficiencysteps[top], 'Belt Mass =', 10**beltMasssteps[top]

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
plt.ylim(456090,456130)

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

weighter = [1./length for dummy in chisteps[chop:]]
histobars = 20

ax1 = plt.subplot(321)
n1, bins1, patches1 = ax1.hist(innerRadsteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax1.set_xlabel(r'$R_{in}$ [AU]', fontsize=16)
ax1.set_ylabel('Fraction', fontsize=16)
#ax1.set_title(r'$R_{in}$')
ax1.grid(True)
plt.axvline(x=innerRadsteps[top], ymin=0, ymax=100, color='k', linewidth=3)

ax2 = plt.subplot(322)
n2, bins2, patches2 = ax2.hist(grainSizesteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax2.set_xlabel(r'log(a [$\mu m$])', fontsize=16)
ax2.set_ylabel('Fraction', fontsize=16)
#ax2.set_title(r'log(a)')
ax2.grid(True)
plt.axvline(x=grainSizesteps[top], ymin=0, ymax=100, color='k', linewidth=3)
'''
ticks2 = numpy.arange(-2.2,-0.8,0.3)
plt.xticks(ticks2)
plt.yticks(numpy.arange(0,0.1,0.02))
'''

ax3 = plt.subplot(323)
n3, bins3, patches3 = ax3.hist(diskMasssteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax3.set_xlabel(r'log($M_D$ [$M_{\oplus}$])', fontsize=16)
ax3.set_ylabel('Fraction', fontsize=16)
#ax3.set_title(r'log($M_D$)')
ax3.grid(True)
plt.axvline(x=diskMasssteps[top], ymin=0, ymax=100, color='k', linewidth=3)
'''
plt.yticks(numpy.arange(0,0.1,0.02))
'''

ax4 = plt.subplot(324)
n4, bins4, patches4 = ax4.hist(grainEfficiencysteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax4.set_xlabel(r'$\beta$', fontsize=16)
ax4.set_ylabel('Fraction', fontsize=16)
#ax4.set_title(r'$\beta$')
ax4.grid(True)
plt.axvline(x=grainEfficiencysteps[top], ymin=0, ymax=100, color='k', linewidth=3)
'''
ticks4 = numpy.arange(0.31,0.45,0.03)
plt.xticks(ticks4)
plt.yticks(numpy.arange(0,0.1,0.02))
'''

ax5 = plt.subplot(325)
#Customization, since this last plot was a problem just for mcmc0704
'''
minBM = -6.1
beltMasssteps2 = [i for i in beltMasssteps[chop:] if i > minBM]
weighter2 = [1./length for i in beltMasssteps2]
'''
n5, bins5, patches5 = ax5.hist(beltMasssteps[chop:], histobars, normed=0, weights=weighter, facecolor='green', alpha=0.75)
ax5.set_xlabel(r'log($M_B$ [$M_{\oplus}$])', fontsize=16)
ax5.set_ylabel('Fraction', fontsize=16)
#ax5.set_title(r'log($M_B$)')
ax5.grid(True)
plt.axvline(x=beltMasssteps[top], ymin=0, ymax=100, color='k', linewidth=3)
'''
ticks5 = numpy.arange(-6.1,-5.7,0.1)
plt.xticks(ticks5)
plt.xticks(numpy.arange(-6.1,-5.7,0.1))
plt.yticks(numpy.arange(0,0.1,0.02))
'''
plt.subplots_adjust(wspace=0.4, hspace=0.7)
plt.suptitle('Probability Distributions from MCMC', fontsize=18)


'''
#Plot the chain

plt.figure(1, figsize=(9,8))

smoothstep = length/100

plt.subplot(331)
plt.plot(range(len(chisteps)), smooth(chisteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'$\chi ^2$')
plt.title(r'$\chi ^2$')
plt.ylim(456090,456130)

plt.subplot(332)
plt.plot(range(len(innerRadsteps)), smooth(innerRadsteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'$R_{in}$ (AU)')
plt.title(r'$R_{in}$')

plt.subplot(333)
plt.plot(range(len(grainSizesteps)), smooth(grainSizesteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'log(a) ($\mu m$)')
plt.title(r'log(a)')

plt.subplot(334)
plt.plot(range(len(diskMasssteps)), smooth(diskMasssteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'log($M_D$) ($M_{\oplus}$)')
plt.title(r'log($M_D$)')

plt.subplot(335)
plt.plot(range(len(grainEfficiencysteps)), smooth(grainEfficiencysteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'$\beta$')
plt.title(r'$\beta$')

plt.subplot(336)
plt.plot(range(len(beltMasssteps)), smooth(beltMasssteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'log($M_B$) ($M_{\oplus}$)')
plt.title(r'log($M_B$)')

plt.subplot(337)
plt.plot(range(len(widthsteps)), smooth(widthsteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'$W_D$ (AU)')
plt.title(r'$W_D$')

plt.subplot(338)
plt.plot(range(len(powerLawsteps)), smooth(powerLawsteps, smoothstep))
plt.xlabel('Steps')
plt.ylabel(r'p')
plt.title(r'p')

plt.subplots_adjust(wspace=0.6, hspace=0.4)

#plt.savefig('chain_'+str(sys.argv[1])+'.png')

#Histograms
plt.figure(2, figsize=(9,7))

weighter = [1./length for dummy in chisteps[chop:]]
histobars = 20

ax1 = plt.subplot(331)
n1, bins1, patches1 = ax1.hist(innerRadsteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax1.set_xlabel(r'$R_{in}$ [AU]')
ax1.set_ylabel('Fraction')
ax1.set_title(r'$R_{in}$')
ax1.grid(True)
plt.axvline(x=innerRadsteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax2 = plt.subplot(332)
n2, bins2, patches2 = ax2.hist(grainSizesteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax2.set_xlabel(r'log(a [$\mu m$])')
ax2.set_ylabel('Fraction')
ax2.set_title(r'log(a)')
ax2.grid(True)
plt.axvline(x=grainSizesteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax3 = plt.subplot(333)
n3, bins3, patches3 = ax3.hist(diskMasssteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax3.set_xlabel(r'log($M_D$ [$M_{\oplus}$])')
ax3.set_ylabel('Fraction')
ax3.set_title(r'log($M_D$)')
ax3.grid(True)
plt.axvline(x=diskMasssteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax4 = plt.subplot(334)
n4, bins4, patches4 = ax4.hist(grainEfficiencysteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax4.set_xlabel(r'$\beta$')
ax4.set_ylabel('Fraction')
ax4.set_title(r'$\beta$')
ax4.grid(True)
plt.axvline(x=grainEfficiencysteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax5 = plt.subplot(335)
n5, bins5, patches5 = ax5.hist(beltMasssteps[chop:], histobars, normed=0, weights=weighter, facecolor='green', alpha=0.75)
ax5.set_xlabel(r'log($M_B$ [$M_{\oplus}$])')
ax5.set_ylabel('Fraction')
ax5.set_title(r'log($M_B$)')
ax5.grid(True)
plt.axvline(x=beltMasssteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax6 = plt.subplot(336)
n6, bins6, patches6 = ax6.hist(widthsteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax6.set_xlabel(r'$W_D$ (AU)')
ax6.set_ylabel('Fraction')
ax6.set_title(r'$W_D$')
ax6.grid(True)
plt.axvline(x=widthsteps[top], ymin=0, ymax=100, color='r', linewidth=2)

ax7 = plt.subplot(338)
n7, bins7, patches7 = ax7.hist(powerLawsteps[chop:], histobars, weights=weighter, normed=0, facecolor='green', alpha=0.75)
ax7.set_xlabel(r'p')
ax7.set_ylabel('Fraction')
ax7.set_title(r'p')
ax7.grid(True)
plt.axvline(x=powerLawsteps[top], ymin=0, ymax=100, color='r', linewidth=2)

plt.subplots_adjust(wspace=0.4, hspace=0.7)
plt.suptitle('Probability Distributions from MCMC', fontsize=16)

#plt.savefig('histo_'+str(sys.argv[1])+'.png')
'''

plt.figure(3)
disk = Disk(innerRadsteps[top], widthsteps[top]+innerRadsteps[top], 10**grainSizesteps[top], 10**diskMasssteps[top], powerLawsteps[top], grainEfficiencysteps[top], 10**beltMasssteps[top])
disk.plotSED()
print 'SED chi-squared =', disk.computeChiSquared()

plt.show()