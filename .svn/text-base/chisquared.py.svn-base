from disk import Disk
import numpy
import time

# DEPRECATED

class ChiComputer:
    
    def __init__(self):
        self.innerRads = [100]#numpy.arange(5, 250, 25)#numpy.arange(5, 200, 5)
        self.outerRads = numpy.arange(100, 140, 3)#numpy.arange(5, 200, 5)
        temp = numpy.arange(-2, 0, .25)#numpy.arange(-2, 1.8, .5)
        self.grainSizes = [10**x for x in temp]
        temp = numpy.arange(-6, -4, .2)
        self.diskMasses = [10**x for x in temp]#numpy.arange(.001, .02, .001)
        self.powerLaws = numpy.arange(.2, 1.99, .3)#numpy.arange(.2, 1.99, .4)
        self.grainEfficiencies = numpy.arange(.2, 1.99, .3)#numpy.arange(.2, 1.99, .4)
    
    def computeChiSquared(self, disk):
        disk.generateInterpolation()
        lamma = disk.sample_lambda
        model_flux = disk.interpol_flux
        actual_flux = disk.sample_flux
        error = disk.sample_error
        chi_squared = 0
        for i in range(len(lamma)):
            # the others are very close and we don't want them to mess up the chi-squared
            if lamma[i] > 1e1:
                chi_squared += ((model_flux[i]-actual_flux[i])/error[i])**2
        return chi_squared
    
    def findBest(self):
        startTime = time.clock()
        ranker = Ranker(20)
        numModels = 0
        for i in self.innerRads:
            for o in self.outerRads:
                # doesn't make sense to have larger inner radius
                if i >= o:
                    continue
                for g in self.grainSizes:
                    print [i, o, g]
                    for m in self.diskMasses:
                        for p in self.powerLaws:
                            for e in self.grainEfficiencies:
                                numModels += 1
                                disk = Disk(i, o, g, m, p, e)
                                chi_squared = self.computeChiSquared(disk)
                                ranker.update((i, o, g, m, p, e), chi_squared)
        endTime = time.clock()
        compTimeHours = (endTime-startTime)/(60**2)
        rankings = ranker.getRankings()
        for score in sorted(rankings.iteritems()):
            print score[0], ':', score[1]
        print 'computation took: ', compTimeHours, 'hours'
        print 'calculated ', numModels, 'models'
        
class Ranker:
    
    """
    This class stores the rankings of the lowest minimums.
    """
    
    def __init__(self, numToTrack):
        self.best = {}
        self.numToTrack = numToTrack
        self.benchmark = float('inf')
    
    def update(self, item, score):
        if score <= self.benchmark:
            # add this one to the ranking
            self.best[score] = item
            if len(self.best) > self.numToTrack:
                # we need to replace one of the rankings
                self.best.pop(self.benchmark)
                self.benchmark = max(self.best.keys())
            elif len(self.best) == self.numToTrack:
                # this only happens once when we need to
                # to set a benchmark
                self.benchmark = max(self.best.keys())
    
    def getRankings(self):
        return self.best.copy()

chi = ChiComputer()
chi.findBest()