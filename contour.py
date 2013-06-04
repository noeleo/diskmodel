import math
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from pylab import *
import sys

'''
This code plots sigma contours of chi^2 in the space of 2 different parameters.

Ex:  python contour.py MCMC_Chains/mcmc_ensemble1113 'ensemble' 'beta'
'''

if len(sys.argv) < 4:
    print 'Error:  I need more arguments.  Give me something like this:'
    print "python contour.py MCMC_Chains/mcmc_ensemble1113 'ensemble' 'beta'"
    exit()

if sys.argv[2] != 'ensemble' and sys.argv[2] != 'mh':
    print 'Error:  Give me a mode as sys.argv[2].  Either "ensemble" or "mh".'
    exit()

if sys.argv[2] == 'ensemble':
    f = open(sys.argv[1],'r')
    print 'Reading in data from', sys.argv[1], 'as an ensemble run.'
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
        innerRadsteps.append(float(info[2]))
        outerRadsteps.append(float(info[3]))
        grainSizesteps.append(float(info[4])) #Note that this is still log.
        diskMasssteps.append(float(info[5]))  #This one too.
        powerLawsteps.append(float(info[6]))
        grainEfficiencysteps.append(float(info[7]))
        beltMasssteps.append(float(info[8]))
        chisteps.append(float(info[9]))
        acceptance.append(float(info[10])) #There's no more chi2; it wasn't very interesting.
        line = f.readline()
    f.close()
    print 'Data read!  Beginning to make boxes...'

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
        acceptance.append(float(info[10])) #Position 9 contains chi2.
        line = f.readline()
    f.close()
    print 'Data read!  Beginning to make boxes...'

#Make Boxes
xaxis = grainSizesteps
xminmax = [0.0,0.7]

if sys.argv[3] == 'beta':
    yaxis = grainEfficiencysteps
    yminmax = [0.2,0.7] 
if sys.argv[3] == 'M_D':
    yaxis = diskMasssteps
    yminmax = [-3.6,-2.3] 
if sys.argv[3] == 'R_in':
    yaxis = innerRadsteps
    yminmax = [55,80] #Inner Radius

density = 30.
dx = (xminmax[1]-xminmax[0])/density
dy = (yminmax[1]-yminmax[0])/density

xarray = numpy.arange(xminmax[0],xminmax[1],dx)
xarray = xarray[0:int(density)]
yarray = numpy.arange(yminmax[0],yminmax[1],dy)
yarray = yarray[0:int(density)]

length = len(xaxis)
chop = int(math.ceil(length*0.10)) #Ignore the first 10% of the chain.

#Each entry constitutes 100/(length-chop) %

z = []
for numy in range(0,int(density)):
    z.append([])
    for numx in range(0,int(density)):
        z[numy].append(0)
print 'Boxes created, attempting to fill.'
        
#Fill Boxes
def Binner(xvalue, yvalue):
    xbin = numpy.digitize(xvalue, xarray)
    ybin = numpy.digitize(yvalue, yarray)
    bin = [xbin,ybin]
    return bin

for i in range(length - chop):
    x_0 = xaxis[chop + i]
    y_0 = yaxis[chop + i]
    if x_0 >= xminmax[0] and x_0 < xminmax[1]:
        if y_0 >= yminmax[0] and y_0 < yminmax[1]:
            bin = Binner([x_0], [y_0])
            z[bin[1] - 1][bin[0] - 1] += 1.

print "Rescaling."
everything = []
for i in range(0,len(z)):
    everything.extend(z[i])
#Rescale so that the sum of all elements is 1
grand_total = sum(everything)
zsort = [i/grand_total for i in everything if i != 0]
for numy in range(0,int(density)):
    for numx in range(0,int(density)):
        z[numy][numx] /= grand_total
zsort = sorted(zsort)
zsort.reverse()

print "Determining contour levels."
running_sum = 0.
for i in range(0,len(zsort)):
    if running_sum < 0.341*2:
        running_sum += zsort[i]
    else:
        sigma1 = zsort[i-1]
        start = i
        break
for j in range(start,len(zsort)):
    if running_sum < 0.477*2:
        running_sum += zsort[j]
    else:
        sigma2 = zsort[j-1]
        start = j
        break
for k in range(start,len(zsort)):
    if running_sum < 0.498*2:
        running_sum += zsort[k]
    else:
        sigma3 = zsort[k-1] 
        start = k
        break

contours = [sigma3,sigma2,sigma1,zsort[0]]

print 'Boxes filled!  Plotting now.'

#Locate Best Fit to overplot
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
    print 'This value is seen at step(s):', bestmodels
top = bestmodels[0] 

fig = plt.figure(figsize=(6,6))

CS = plt.contour(xarray,yarray,z,contours, colors=('k',), linewidths = (2,))
CF = plt.contourf(xarray,yarray,z,contours, colors=('b','r','y'))
plt.plot(xaxis[top],yaxis[top], 'o', color='w', markersize=8)
gcf().subplots_adjust(bottom=0.125)
gcf().subplots_adjust(left=0.175)
plt.tick_params(labelsize=16)
plt.xlabel(r'log($a$ $[\mu m]$)', fontsize=18)
if sys.argv[3] == 'beta':
    plt.ylabel(r'$\beta$', fontsize=18)
if sys.argv[3] == 'M_D':
    plt.ylabel(r'log($M_D$ $[M_\oplus]$)', fontsize=18)
if sys.argv[3] == 'R_in':
    plt.ylabel(r'$R_{in}$ [AU]', fontsize=18)
plt.xticks(numpy.arange(xminmax[0],xminmax[1],0.1))


#Buenzli's radius
if sys.argv[3] == 'R_in':
    ax = fig.add_subplot(111)
    rectangle = Rectangle((-5,60.4), 10, 1.7, alpha=0.5, facecolor="grey")
    ax.add_patch(rectangle)
    
    plt.axhline(y=62.1, xmin=-5, xmax=5, color='grey', linewidth=2)
    plt.axhline(y=60.4, xmin=-5, xmax=5, color='grey', linewidth=2)


#plt.savefig('ContourGSIR.eps')
plt.show()

        