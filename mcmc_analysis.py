import math
from numpy import *
import matplotlib.pyplot as plt
import sys
from disk_as import Disk
from matplotlib.patches import Rectangle
import pylab

'''
A more elegantly coded version of mcmcreader, which I got sick of.
This is the presentation part of the mcmc_ensemble code using data saved into a text file.

Ex:  python mcmcreader.py MCMC_Chains/mcmc_ensemble1113
'''

class Ensemble:
    
    def __init__(self, filename, chop):
        #parameters = [R_in,R_out,a,M_D,p,beta,M_B,chi^2,accept]
        #dimensions are walker, trial, parameter
        self.parameters = ['R_in','R_out','log(a)','log(M_D)','p','beta','log(M_B)']
        self.latex_params = [r'$R_{in}$ (AU)', r'$R_{out}$ (AU)', r'log(a) ($\mu m$)', r'log($M_D$) ($M_{\oplus}$)', r'$p$', r'$\beta$', r'log($M_B$) ($M_{\oplus}$)', r'$\chi ^2$']
        self.n_param = len(self.parameters)
        self.n_walkers = 100
        #Recreate the Ensemble
        ensemble = [[] for i in range(self.n_walkers)]
        f = open(filename,'r')
        print 'Reading in data from', filename
        line = f.readline()
        while line != '':
            info = [float(x) for x in line.split()]
            ensemble[int(info[1])].append(array(info[2:self.n_param+4]))
            line = f.readline()
        f.close()
        self.ensemble = array(ensemble)
        self.chop = int(chop) #Ignore the first part of the ensemble.
        
    def Quick_Stats(self):
        for param in range(self.n_param):
            print 'For', self.parameters[param], ':  Mean =', average(self.ensemble[:,self.chop:,param]), \
            ', Median =', median(self.ensemble[:,self.chop:,param]), ', STD =', std(self.ensemble[:,self.chop:,param])
        print 'Acceptance fraction =', average(self.ensemble[:,self.chop:,self.n_param+1])

    def Find_Best_Fit(self):
        chi_min = self.ensemble[:,:,self.n_param].min()
        mindex = where(self.ensemble[:,:,self.n_param] == chi_min)
        print 'The minimum chi^2 was', chi_min, 'which occurred at walker', mindex[0], 'trial(s)', mindex[1]
        self.bestmodel = self.ensemble[mindex[0][0],mindex[1][0],0:self.n_param]
    
    def Sigma_Calc(self):
        self.sigmas = []
        for param in arange(self.n_param):
            deltas = [abs(x - self.bestmodel[param]) for x in self.ensemble[:,self.chop:,param].ravel()]
            deltas.sort()
            self.sigmas.append(deltas[int(0.6827*len(deltas))])
        
    def Binner(self, n_bins):
        bins = []
        bar_heights = []
        for param in arange(self.n_param):
            range = [self.bestmodel[param] - 5*self.sigmas[param], self.bestmodel[param] + 5*self.sigmas[param]]
            histoheights, histobins = histogram(self.ensemble[:,self.chop:,param], bins=n_bins, normed=False, range=range)
            bins.append(histobins)
            bar_heights.append(histoheights)
        return [bins,bar_heights]
    
    def acor(self):
        #This one takes a long time...
        def C_f(array,T):
            avg = average(array)
            M = len(array)
            return sum([(array[T+m]-avg)*(array[m]-avg) for m in range(M-T)])/(M-T)
        
        def tau(array):
            return 1 + 2*sum([C_f(array,T) for T in [i+1 for i in range(len(array)-1)]])/C_f(array,0)
        
        print 'Calculating autocorrelation times.'
        taus = []
        for param in range(self.n_param):
            taus.append(tau(self.ensemble[:,:,param].transpose().ravel()))
        print 'Autocorrelation times for R_in, a, M_D, beta, and M_B:', taus
        exit()
        
    def Chain_Plot(self):
        plt.figure(1, figsize=(9,8))
        f_index = 0
        for param in [7,0,2,3,5,6]: #chi^2, R_in, a, M_D, beta, M_B
            plt.subplot(320+f_index)
            plt.plot(range(len(self.ensemble[0,:,param])), average(self.ensemble[:,:,param],axis=0))
            plt.xlabel('Steps', fontsize=16)
            plt.ylabel(self.latex_params[param], fontsize=16)
            f_index += 1
        plt.subplots_adjust(wspace=0.6, hspace=0.4)
        
    def Histo_Plot(self, n_bins):
        #I highly doubt this works...
        histo = self.Binner(n_bins)
        plt.figure(2, figsize=(9,7))
        f_index = 1
        for param in [0,2,3,5,6]:
            bins = histo[0][param]
            heights = histo[1][param]
            plt.subplot(320+f_index)
            pylab.plot(.5*(bins[1:]+bins[:-1]), [float(x)/len(self.ensemble[:,self.chop:,param].ravel()) for x in heights], linewidth=2)
            plt.axvline(x=self.bestmodel[param], ymin=0, ymax=100, color='k', linewidth=3)
            if param==0:
                plt.axvline(x=62.1, ymin=0, ymax=100, color='c', linewidth=3)
                plt.axvline(x=60.4, ymin=0, ymax=100, color='c', linewidth=3)
            plt.xlabel(self.latex_params[param], fontsize=16)
            plt.ylabel('Fraction', fontsize=16)
            plt.yticks(arange(0,0.17,0.04))
            f_index += 1
        plt.subplots_adjust(wspace=0.4, hspace=0.7)

    def SED_Plot(self):
        plt.figure(3)
        disk = Disk(self.bestmodel[0], self.bestmodel[1], 10**self.bestmodel[2], 10**self.bestmodel[3], self.bestmodel[4], self.bestmodel[5], 10**self.bestmodel[6])
        disk.plotSED()
        print 'SED chi^2 =', disk.computeChiSquared()
        print 'T_avg =', disk.calculateGrainTemperature(0.5*(self.bestmodel[0]+self.bestmodel[1])*1.496e11)
    
    def Modes(self):
        mode_finder = self.Binner(30)
        self.modes = []
        for param in arange(self.n_param):
            bins = mode_finder[0][param]
            heights = mode_finder[1][param]
            self.modes.append(bins[where(array(heights)==array(heights).max())])
    
    def Analyze(self):
        print '-------------------------------------------------'
        self.Quick_Stats()
        print '-------------------------------------------------'
        self.Find_Best_Fit()
        print '-------------------------------------------------'
        self.Sigma_Calc()
        self.Modes()
        for param in range(self.n_param):
            print self.parameters[param], '=', self.bestmodel[param], '+/-', self.sigmas[param]
        print 'And for the modes, we have', [float(x) for x in self.modes]
        self.Chain_Plot()
        self.Histo_Plot(40)
        self.SED_Plot()
        plt.show()