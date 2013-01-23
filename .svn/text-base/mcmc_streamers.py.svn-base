import math
import random
import numpy
import sys
from disk import Disk
from visgen_streamers import VisibilityGenerator

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
    
#Step 1:  Initialize the Chain.  
flux1 = 0.0066
streamerpart1 = 0.5
diskpart1 = 1.0 - streamerpart1

betaFX = 0.001
betaSP = 0.15

#Keep Track of Chain
FXsteps = [flux1]
SPsteps = [streamerpart1]
DPsteps = [diskpart1]
chisteps = []
acceptance = []

disk = Disk(70.5, 74.0, 0.00358, 7.865e-6, 0.5, 0.323, 2e-8)
#The long decimals here are to normalize the disk and streamer fluxes to 1 before multiplying.
#This means that the disk flux will equal flux1*diskpart1 and the streamer flux will equal flux1*streamerpart1.
diskscale1 = flux1*diskpart1*151.805886744
streamerscale1 = flux1*streamerpart1*113.35516812
vis = VisibilityGenerator(im_width, theta_i, theta_pa, diskscale1, streamerscale1, 'mcmc_streamers.fits')

#Step 2:  Generate a Trial State by randomly selecting a parameter to vary
count = 0
while count < stop:    
    choice = random.randint(1,2)
    if choice == 1:
        flux2 = random.gauss(flux1, betaFX)
        while flux2 < 0:
            flux2 = random.gauss(flux1, betaFX)
    else:
        flux2 = flux1
    if choice == 2:
        streamerpart2 = random.gauss(streamerpart1, betaSP)
        while streamerpart2 < 0 or streamerpart2 > 1:
            streamerpart2 = random.gauss(streamerpart1, betaSP)
        diskpart2 = 1 - streamerpart2
    else:
        streamerpart2 = streamerpart1
        diskpart2 = diskpart1
    
    #Step 3:  Compute Chi-Squared for nth and (n+1)th states
    diskscale1 = flux1*diskpart1*151.805886744
    streamerscale1 = flux1*streamerpart1*113.35516812  
    vis.changeParameters(diskscale1, streamerscale1)
    chi1 = vis.computeChiSquared(disk)
    
    diskscale2 = flux2*diskpart2*151.805886744
    streamerscale2 = flux2*streamerpart2*113.35516812
    vis.changeParameters(diskscale2, streamerscale2)
    chi2 = vis.computeChiSquared(disk)
    
    #Step 4:  Calculate f'(x)/f(x_n)
    prob = math.e**(-0.5*(chi2 - chi1))
    
    #Step 5:  Draw a random number between 0 and 1
    dice = random.random()
    
    #Step 6:  Determine whether to keep trial state
    alphatest = min(prob,1)
    if dice <= alphatest:
        flux1 = flux2
        streamerpart1 = streamerpart2
        diskpart1 = diskpart2
        accept = 1
    else:
        accept = 0
    
    FXsteps.append(flux1)
    SPsteps.append(streamerpart1)
    DPsteps.append(diskpart1)
    chisteps.append(chi1) 
    acceptance.append(accept)
    
    #Step 7:  Increase counter by 1 and save
    count = count + 1
    print 'Step =', count
    print 'Total Flux =', flux1, 'Streamer Contribution =', streamerpart1, 'chi1 =', chi1, 'chi2 =', chi2, 'accept =', accept
    if (len(sys.argv) == 2):
        save.write(str(count)+' '+str(flux1)+' '+str(streamerpart1)+' '+str(diskpart1)+' '+str(chi1)+' '+str(chi2)+' '+str(accept)+'\n')
        
#Mean and Standard Deviation
chop = int(math.ceil(stop*0.1)) #Ignore the first 10% of the chain.
print 'Number of Steps =', stop
print 'Now ignoring the first 10% of the chain...'
print 'Total Flux: ', 'Mean =', numpy.average(FXsteps[chop:]), 'STD =', numpy.std(FXsteps[chop:])
print 'Streamer Contribution: ', 'Mean =', numpy.average(SPsteps[chop:]), 'STD =', numpy.std(SPsteps[chop:])
print 'Acceptance Rate =', float(numpy.sum(acceptance))/float(stop)
if (len(sys.argv) == 2):
    save.close()