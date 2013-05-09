'''
This is a code to implement an affine-invariant MCMC algorithm.  This works well when parameters are degenerate.
Basically, you have an ensemble of "walkers" that each evolve over time.  To determine a trial state for one of 
the walkers, you randomly choose a different walker and move toward it.
'''

import math
import random
import numpy
import sys
from disk_as import Disk
from visgen import VisibilityGenerator

#Important Disk Parameters for the FITS file needed by visgen.
im_width = 512
theta_i = 84.3
theta_pa = 70.3

#Make File to Save
if (len(sys.argv) == 1):
    print 'Error: No file name given.'
    exit()
if (len(sys.argv) == 2):
    save = open(sys.argv[1],'w')
    print 'Saving data to', sys.argv[1]

#Set number of walkers and then number of trials for each.
n_walkers = 100
n_trials = 1000
a = 2.2           #This is just an adjustable parameter that Foreman-Mackey et al., 2012 suggest should be 2.  My accept ratio was too high.
    
#Initialize the Ensemble.  These values are from a chi-by-eye
innerRad0 = 68.
outerRad0 = innerRad0*1.05  #Set equal to innerRad x 1.05
grainSize0 = 0.5  #Will be raised to the 10th power
diskMass0 = -2.7  #Will be raised to the 10th power
powerLaw0 = 1.0  #Currently fixed
grainEfficiency0 = 0.5
beltMass0 = -6.0 #Will be raised to the 10th power

#These widths are roughly equal to the stddev from mcmc_ensemble0423.
betaIR = 2.0
betaOR = 2.0     #We don't actually use this, since the outer radius is fixed to 1.05 times the inner radius.
betaGS = 0.05
betaDM = 0.1
betaPL = 0.01   #Never used this; that number is a placeholder.
betaGE = 0.05 
betaBM = 0.05

#Create a disk and visibility generator
disk = Disk(innerRad0, outerRad0, 10.**grainSize0, 10.**diskMass0, powerLaw0, grainEfficiency0, 10.**beltMass0)
vis = VisibilityGenerator(im_width, theta_i, theta_pa, 'mcmc.fits')

chi0 = 0.0       #A placeholder. 
acceptance0 = 1  #This is a dimension used to track the acceptance probability, which we want to be \in (0.2,0.5)

parameters = [innerRad0, outerRad0, grainSize0, diskMass0, powerLaw0, grainEfficiency0, beltMass0, chi0, acceptance0] #Each walker looks like this.
betas = [betaIR, betaOR, betaGS, betaDM, betaPL, betaGE, betaBM]                                                      #Gaussian Widths
ensemble = [[[x for x in parameters]] for y in range(n_walkers)]                                                      #An ensemble of walkers.

#Distribute the walkers.
print "Distributing walkers and calculating initial chi^2.  This could take many minutes."
for walker in range(n_walkers):
    for param in [0,2,3,5,6]:  #R_in, a, M_D, \beta, M_B
        ensemble[walker][0][param] += betas[param]*random.gauss(0,1)

    #Check that inner radius and grain efficiency are positive.
    while ensemble[walker][0][0] <=0:
        ensemble[walker][0][0] = parameters[0] + betas[0]*random.gauss(0,1)
    while ensemble[walker][0][5] <=0:
        ensemble[walker][0][5] = parameters[5] + betas[5]*random.gauss(0,1)    
    
    #Set the outer radius based on inner radius.
    ensemble[walker][0][1] = ensemble[walker][0][0]*1.05  
    
    #Calculate chi^2 for each walker.
    disk.changeParameters(ensemble[walker][0][0], ensemble[walker][0][1], 10.**(ensemble[walker][0][2]), 10.**(ensemble[walker][0][3]), \
                          ensemble[walker][0][4], ensemble[walker][0][5], 10.**(ensemble[walker][0][6]))
    sedchi = disk.computeChiSquared()
    vischi = vis.computeChiSquared(disk)
    ensemble[walker][0][7] = sedchi + vischi
    print "Walker", walker+1, "of", n_walkers, "established."
print "Initial walkers created!"

'''
The ensemble is a list of lists of lists.  It has "n_walkers" walkers (models), each of which will be appended to contain "n_trials" trials.
'''

#Advance the walkers!
for trial in range(n_trials):
    for walker in range(n_walkers):  
        print "Trial", trial+1, "of", n_trials, ", walker", walker+1, "of", n_walkers, ":"
        #Find a random other walker.
        random_other = int(random.random()*n_walkers)
        while random_other == walker:
            random_other = int(random.random()*n_walkers)
        zz = ((a - 1.)*random.uniform(0.,1.) + 1.)**2/a               
        
        #Propose a step along the line between them.
        parameter_difference = [ensemble[random_other][trial][param]-ensemble[walker][trial][param] for param in range(len(ensemble[walker][trial])-2)] #-2 to exclude chi2, acceptance
        proposition = [ensemble[walker][trial][param] + zz*parameter_difference[param] for param in range(len(ensemble[walker][trial])-2)]
        
        while proposition[0] <= 0 or proposition[5] <= 0:
            zz = ((a - 1.)*random.uniform(0.,1.) + 1.)**2/a
            parameter_difference = [ensemble[random_other][trial][param]-ensemble[walker][trial][param] for param in range(len(ensemble[walker][trial])-2)] #-2 to exclude chi2, acceptance
            proposition = [ensemble[walker][trial][param] + zz*parameter_difference[param] for param in range(len(ensemble[walker][trial])-2)]
        print "Proposition:", proposition
        #Calculate chi^2.
        disk.changeParameters(proposition[0], proposition[1], 10.**(proposition[2]), 10.**(proposition[3]), \
                              proposition[4], proposition[5], 10.**(proposition[6]))
        sedchi = disk.computeChiSquared()
        vischi = vis.computeChiSquared(disk)
        prop_chi = sedchi + vischi
        
        #From that, get probability of acceptance.  Here, we feel like using the log.
        log_prob_new = -0.5*prop_chi
        log_prob_old = -0.5*ensemble[walker][trial][7]
        log_prob_diff = (len(proposition) - 1)*math.log(zz) + log_prob_new - log_prob_old
        
        #Draw a random number.
        dice = random.random()
        while dice == 0:
            dice = random.random()
        log_dice = math.log(dice)
        
        if log_dice <= log_prob_diff:
            #Accept new state.
            proposition.append(prop_chi)
            proposition.append(1)
            ensemble[walker].append(proposition)
        else:
            #Keep old state.
            ensemble[walker].append(ensemble[walker][trial])
            ensemble[walker][trial+1][8] = 0
        print 'R_in =', ensemble[walker][trial+1][0], 'log(a) =', ensemble[walker][trial+1][2], 'log(M_D) =', ensemble[walker][trial+1][3], \
        'Beta =', ensemble[walker][trial+1][5], 'log(M_B) =', ensemble[walker][trial+1][6], 'Chi^2 =', ensemble[walker][trial+1][7], \
        'Accept =', ensemble[walker][trial+1][8]
        if (len(sys.argv) == 2):
            save.write(str(trial)+' '+str(walker)+' '+str(ensemble[walker][trial+1][0])+' '+str(ensemble[walker][trial+1][1])+' '+str(ensemble[walker][trial+1][2])\
                       +' '+str(ensemble[walker][trial+1][3])+' '+str(ensemble[walker][trial+1][4])+' '+str(ensemble[walker][trial+1][5])+' '+\
                       str(ensemble[walker][trial+1][6])+' '+str(ensemble[walker][trial+1][7])+' '+str(ensemble[walker][trial+1][8])+'\n')
print "All done!"
if (len(sys.argv) == 2):
    save.close()
    print "Data saved to", sys.argv[1], '.'
