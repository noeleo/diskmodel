import math
import random
import numpy
import os
import matplotlib.pyplot as plt
import sys
from disk import Disk
from visgen import VisibilityGenerator

'''
Currently fixing width to 65 AU
'''

#Important Disk Parameters for the FITS file
im_width = 512
theta_i = 84.3
theta_pa = 70.3

#Set Length of Chain
stop = 10**5

#Make File to Save
if (len(sys.argv) == 1):
    print 'You are not saving, fyi.  If you want to save, the file name is sys.argv[1].'
if (len(sys.argv) == 2):
    save = open(sys.argv[1],'w')
    
#Step 1:  Initialize the Chain.  These values are from a chi-by-eye
innerRad1 = 45
width1 = 65  #Set equal to innerRad x 1.05
grainSize1 = -1  #Will be raised to the 10th power
diskMass1 = -3.7  #Will be raised to the 10th power
powerLaw1 = 0.5  #Currently fixed
grainEfficiency1 = 0.323
beltMass1 = -5.7 #Will be raised to the 10th power

betaIR = 8
#betaWD = 1
betaGS = 0.1
betaDM = 0.25
#betaPL = 
betaGE = 0.08 
betaBM = 1.0 #Just a guess

#Keep Track of Chain
innerRadsteps = [innerRad1]
widthsteps = [width1]
grainSizesteps = [grainSize1] 
diskMasssteps = [diskMass1]  
powerLawsteps = [powerLaw1]
grainEfficiencysteps = [grainEfficiency1]
beltMasssteps = [beltMass1]
chisteps = []
acceptance = []

disk = Disk(innerRad1, innerRad1 + width1, 10**grainSize1, 10**diskMass1, powerLaw1, grainEfficiency1, 10**beltMass1)
vis = VisibilityGenerator(im_width, theta_i, theta_pa, 'mcmc.fits')

#Step 2:  Generate a Trial State by randomly selecting a parameter to vary
count = 0
while count < stop:    
    choice = random.randint(1,5)
    if choice == 1:
        innerRad2 = random.gauss(innerRad1, betaIR)
        while innerRad2 <= 0:
            innerRad2 = random.gauss(innerRad1, betaIR)
    else:
        innerRad2=innerRad1
    if choice == 2:
        grainSize2 = random.gauss(grainSize1, betaGS)
    else:
        grainSize2=grainSize1
    if choice == 3:
        diskMass2 = random.gauss(diskMass1, betaDM)
    else:
        diskMass2 = diskMass1
    powerLaw2 = powerLaw1 #Currently fixed
    if choice == 4:
        grainEfficiency2 = random.gauss(grainEfficiency1, betaGE)
        while grainEfficiency2 <= 0:
            grainEfficiency2 = random.gauss(grainEfficiency1, betaGE)
    else:
        grainEfficiency2 = grainEfficiency1
    if choice == 5:
        beltMass2 = random.gauss(beltMass1, betaBM)
    else:
        beltMass2 = beltMass1
    if choice == 6:
        width2 = random.gauss(width1, betaWD)
        while width2 <= 0:
            width2 = random.gauss(width1, betaWD)
    else:
        width2 = width1
        
    #Step 3:  Compute Chi-Squared for nth and (n+1)th states  
    disk.changeParameters(innerRad1, innerRad1+width1, 10**grainSize1, 10**diskMass1, powerLaw1, grainEfficiency1, 10**beltMass1)
    sedchi1 = disk.computeChiSquared()
    vischi1 = vis.computeChiSquared(disk)
    chi1 = sedchi1 + vischi1
    
    disk.changeParameters(innerRad2, innerRad2+width2, 10**grainSize2, 10**diskMass2, powerLaw2, grainEfficiency2, 10**beltMass2)
    sedchi2 = disk.computeChiSquared()
    vischi2 = vis.computeChiSquared(disk)
    chi2 = sedchi2 + vischi2
    
    #Step 4:  Calculate f'(x)/f(x_n)
    prob = math.e**(-0.5*(chi2 - chi1))
    
    #Step 5:  Draw a random number between 0 and 1
    dice = random.random()
    
    #Step 6:  Determine whether to keep trial state
    alphatest = min(prob,1)
    if dice <= alphatest:
        innerRad1 = innerRad2
        width1 = width2
        grainSize1 = grainSize2
        diskMass1 = diskMass2
        #powerLaw1 = powerLaw2
        grainEfficiency1 = grainEfficiency2
        beltMass1 = beltMass2
        accept = 1
    else:
        accept = 0
    
    innerRadsteps.append(innerRad1)
    widthsteps.append(width1)
    grainSizesteps.append(grainSize1)
    diskMasssteps.append(diskMass1)
    powerLawsteps.append(powerLaw1)
    grainEfficiencysteps.append(grainEfficiency1)
    beltMasssteps.append(beltMass1)
    chisteps.append(chi1) 
    acceptance.append(accept)
    
    #Step 7:  Increase counter by 1 and save
    count = count + 1
    print 'Step =', count
    print 'IR =', innerRad1, 'Width =', width1, 'log(GS) =', grainSize1, 'log(DM) =', diskMass1, 'GE =', grainEfficiency1, 'log(BM) =', beltMass1, 'chi1 =', chi1, 'chi2 =', chi2, "accept =", accept
    if (len(sys.argv) == 2):
        save.write(str(count)+' '+str(innerRad1)+' '+str(width1)+' '+str(grainSize1)+' '+str(diskMass1)+' '+str(powerLaw1)+' '+str(grainEfficiency1)+' '+str(beltMass1)+' '+str(chi1)+' '+str(chi2)+' '+str(accept)+'\n')
        
#Mean and Standard Deviation
chop = int(math.ceil(stop*0.1)) #Ignore the first 10% of the chain.
print 'Number of Steps =', stop
print 'Now ignoring the first 10% of the chain...'
print 'Inner Radius: ', 'Mean =', numpy.average(innerRadsteps[chop:]), 'STD =', numpy.std(innerRadsteps[chop:])
print 'Width:  ', 'Mean =', numpy.average(widthsteps[chop:]), 'STD =', numpy.std(widthsteps[chop:])
print 'log(Grain Size): ', 'Mean =', numpy.average(grainSizesteps[chop:]), 'STD =', numpy.std(grainSizesteps[chop:])
print 'log(Disk Mass): ', 'Mean =', numpy.average(diskMasssteps[chop:]), 'STD =', numpy.std(diskMasssteps[chop:])
print 'Grain Efficiency: ', 'Mean =', numpy.average(grainEfficiencysteps[chop:]), 'STD =', numpy.std(grainEfficiencysteps[chop:])
print 'log(Belt Mass): ', 'Mean =', numpy.average(beltMasssteps[chop:]), 'STD =', numpy.std(beltMasssteps[chop:])
print 'Acceptance Rate =', float(numpy.sum(acceptance))/float(stop)
if (len(sys.argv) == 2):
    save.close()

#For display, use mcmcreader.py
'''
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
plt.ylabel(r'Grain Emissivity (Units?)')
plt.title(r'Grain Emissivity')

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
ax4.set_xlabel('Grain Efficiency (Units?)')
ax4.set_ylabel('Probability')
ax4.set_title('Grain Efficiency')
ax4.grid(True)

plt.subplots_adjust(wspace=0.5, hspace=0.5)

#plt.savefig('histo.ps')
#plt.show()
'''