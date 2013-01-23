import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from scipy.special import gammainc

'''
Interprets text files from streamergrid.py in order to map chi^2 (or delta chi^2) vs. F_t
save.write(str(SC)+' '+str(F_t)+' '+str(chi_t)+'\n')
'''

#For SG_SCFT_0724, x = SC and y = F_t

xvalues = []
yvalues = []
chivalues = []

f = open(sys.argv[1],'r')
line = f.readline()
while line != '':
    info = line.split()
    xvalues.append(float(info[0]))
    yvalues.append(float(info[1]))
    chivalues.append(float(info[2]))
    line = f.readline()
f.close()

size_squared = len(chivalues)
size = int(math.sqrt(size_squared)) #This only makes sense if you have a square grid

#Here, I make an array of the best chi^2 values for each value of x.
chibest = []
bestmodel = [0,0,1e99] #[SC,F_t,chisq]
for i in range(0,size):
    chibest.append(1e99)
for x in range(0,size):
    for y in range(0,size):
        if chivalues[size*x+y] < chibest[y]:
            chibest[y] = chivalues[size*x+y]
        if chivalues[size*x+y] < bestmodel[2]:
            bestmodel[0] = x
            bestmodel[1] = y
            bestmodel[2] = chivalues[size*x+y]
#Now, I obtain the probability of each state
DoF = 246244 - 2 #2*number of visibilities - 2 parameters to fit in the grid

def chisq_pdf(chi2, DoF):
    arg = (chi2 - bestmodel[2])#/math.sqrt(2*DoF) #This quantity converges to the normal distribution as DoF --> Infinity, and DoF is pretty huge.
    #prob = 1/math.sqrt(2*math.pi)*math.e**(-0.5*arg**2)
    prob = 1.0 - gammainc(float(2)/2, arg/2)
    return prob

def p_value(chi2, DoF):
    return 1.0 - gammainc(DoF/2, chi2/2)
    
chi2_best = chisq_pdf(bestmodel[2], DoF)
probabilities = [chisq_pdf(i,DoF) for i in chibest]
#n_sigma = [n_sigma(i) for i in probabilities]

print xvalues[0:size]
print chibest

#Plot raw chi^2 values
plt.figure(1) 
plt.plot(xvalues[0:size], [i-bestmodel[2] for i in chibest], '.-', linewidth=2, color='b')
plt.axvline(x=0.585, ymin=0, ymax=1e9, color='r', linewidth=2)
plt.title(r'Best $\chi ^2$ Values', fontsize=18)
plt.xlabel('Fraction of Flux in Streamers', fontsize=16)
plt.ylabel(r'$\Delta \chi ^2$', fontsize=16)

#Plot probabilities
plt.figure(2)
plt.plot(xvalues[0:size], [p_value(i-bestmodel[2],2) for i in chibest], '.-', linewidth=2)
plt.axvline(x=0.585, ymin=0, ymax=1e9, color='r', linewidth=2)
plt.title(r'Probability', fontsize=18)
plt.xlabel('Fraction of Flux in Streamers', fontsize=16)
plt.ylabel(r'$P$', fontsize=16)
plt.show()

'''
#Plot n_sigma
plt.figure(2)
plt.plot(xvalues[0:size], n_sigma, '.-', linewidth=2)
plt.title(r'Number of Standard Deviations', fontsize=18)
plt.xlabel('Fraction of Flux in Streamers', fontsize=16)
plt.ylabel(r'$N_\sigma$', fontsize=16)
plt.show()
'''

