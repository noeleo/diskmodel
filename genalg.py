import random
from disk import Disk
from visgen import VisibilityGenerator
import sys
import os
import time

# need to ssh with -X parameter because of chi plotting
class Genetic:
    
    # TODO: Gaussian instead of linear for some parameters
    # TEMP: width is now just 1.05*inner
    
    def __init__(self):
        # these 2 numbers cannot be changed without messing up the size of
        # each successive generation because of the way things are bred
        # set the number of top phenotypes to breed from
        self.numTop = 10
        # set the population size by combining each top phenotype
        # with another top twice, wihout breeding with itself
        # = n*(n-1) combinations and we keep the top n so = n^2
        self.popSize = self.numTop**2
        
        # set the ranges in the form [lower, upper, vary, distribution]
        # distribution can be linear/log
        ranges = []
        ranges.append([1.301, 2.301, 0.1, 'log']) # inner Rad 20 to 200 [52.99,53.01,.001,'linear'])
        # TEMP: width is now just 1.05*inner
        ranges.append([1.0, 200.0, 25.0, 'linear']) # width
        ranges.append([-4.0, 3.0, 0.25, 'log']) # grain size [.029,.031,.001,'linear'])
        ranges.append([-5.0, -2.0, 0.25, 'log']) # disk mass [1e-5, 1e-2, 1e-4, 'linear']
        ranges.append([0.1, 1.0, 0.1, 'linear']) # power law [.2, 1.99, .1, 'linear']
        ranges.append([0.2, 3.5, 0.2, 'linear']) # grain efficiency
        self.ranges = ranges
        self.vars = len(ranges)
        # first construct random genotype population
        # random linear scale selection
        self.phenotypes = {}
        population = []
        for i in range(self.popSize):
            phenotype = []
            for g in range(self.vars):
                phenotype.append(random.uniform(ranges[g][0], ranges[g][1]))
            population.append(tuple(phenotype))
        self.population = population
        
        # set visibility data
        self.image_width = 512
        self.inclination_angle = 84.3
        self.position_angle = 70.3
    
    """
    returns a tuple with first element list of phenotypes sorted by
    chi values, second element dictionary mapping phenotypes to chi
    """
    def evaluatePopulation(self):
        # initialize dummy Disk
        disk = Disk(5, 10, .1, .1, .1, .1)
        vis = VisibilityGenerator(self.image_width, self.inclination_angle, self.position_angle, 'genalg.fits')
        chiToPhen = {}
        chis = []
        argChis = []
        # find chi squared for all phenotypes
        for phenotype in self.population:
            # adjust for log scales
            args = list(phenotype)
            for i in range(self.vars):
                if self.ranges[i][3] == 'log':
                    args[i] = 10**args[i]
            # outer radius = inner radius + width
            #args[1] = args[0]+args[1]
            # TEMP: outer = inner*1.05
            args[1] = 1.05*args[0]
            disk.changeParameters(*args)
            # compute and store the chi-squared values
            # find the SED chi-squared
            sed_chi = disk.computeChiSquared()
            # find the visibility chi-squared
            vis_chi = vis.computeChiSquared(disk)
            # combine and store the overall chi-squared
            chi = sed_chi + vis_chi
            chiToPhen[chi] = phenotype
            chis.append(chi)
            # used for output and determining minimum
            self.phenotypes[tuple(args)] = chi
            # for status file
            argChis.append(tuple(args) + (sed_chi, vis_chi, chi))
        # sort the chi values
        chis.sort()
        # create a list of phenotypes from best to worst
        phens = []
        phenToChi = {}
        for chi in chis:
            phen = chiToPhen[chi]
            phens.append(phen)
            phenToChi[phen] = chi
        
        phens = tuple(phens)
        return (phens, argChis)
    
    """
    run the genetic algorithm
    """
    def runAlgorithm(self, numIterations, outputFile, statusFile):
        # NOT GUARANTEED TO PRODUCE SAME POPULATION SIZE EACH TIME
        print "generations completed:", "0/"+str(numIterations),"\r",
        sys.stdout.flush()
        status = open('genetic/' + statusFile, 'w')
        for i in xrange(numIterations):
            start = time.time()
            # evaluate the current population
            phenRanked, argChi = self.evaluatePopulation()
            # write to the status file
            for valz in argChi:
                outStr = ''
                for par in valz:
                    outStr += str(par) + ' '
                status.write(outStr + '\n')
            status.flush()
            # produce offspring only from the top ranked phenotypes
            bestPop = phenRanked[:self.numTop]
            newPopulation = []
            for a in range(self.numTop):
                for b in range(self.numTop):
                    # do not breed with yourself
                    if a==b: continue
                    newPhen = []
                    for y in range(self.vars):
                        rand = random.random()
                        if rand < 0.5:
                            newPhen.append(bestPop[a][y])
                        else:
                            newPhen.append(bestPop[b][y])
                        # add some randomness
                        rand = random.random()
                        if rand < 0.25:
                            newPhen[y] = random.uniform(max(self.ranges[y][0], newPhen[y]-self.ranges[y][2]),
                                                        min(self.ranges[y][1], newPhen[y]+self.ranges[y][2]))
                    newPopulation.append(tuple(newPhen))
            for phen in bestPop:
                newPhen = list(phen)
                for y in range(self.vars):
                    rand = random.random()
                    if rand < 0.25:
                        newPhen[y] = random.uniform(max(self.ranges[y][0], newPhen[y]-self.ranges[y][2]),
                                                    min(self.ranges[y][1], newPhen[y]+self.ranges[y][2]))
                newPopulation.append(tuple(newPhen))
            self.population = newPopulation
            # display some information about how far along we are
            end = time.time()
            curGenTime = (end-start)/60
            timeLeft = curGenTime*(numIterations-i-1)
            print "generations completed:", str(i+1)+"/"+str(numIterations), 'with %.2f minutes left' % timeLeft, "\r",
            sys.stdout.flush()
        status.close()
        print ""
        
        print "Finding minimum..."
        chip = {}
        chis = []
        # just temporarily finding the minimum
        for phen, chi in self.phenotypes.items():
            chis.append(chi)
            chip[chi] = phen
        chis.sort()
        
        print "Writing to file..."
        # write everything to a file
        fname = 'genetic/' + outputFile
        if os.path.exists(fname):
            #print 'overwriting ' + fname
            pass
        f = open(fname, 'w')
        for chi in chis:
            phen = chip[chi]
            st = ''
            for p in phen:
                st += str(p) + ' '
            f.write(st + str(chi) + '\n')
        f.close()
        
        print "Done!"        
        # return the top fit as a tuple
        return (chis[0], chip[chis[0]])

iterations = int(sys.argv[1])
out_file = sys.argv[2]
status_file = sys.argv[3]
g = Genetic()
g.runAlgorithm(iterations, out_file, status_file)
