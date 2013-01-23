from disk import Disk
from visgen import VisibilityGenerator
import numpy
import time
import sys
import os

"""
Minimizes the chi-squared value of a disk model
You can run this script directly to compute the chi-squared values
for a grid of parameters as specified in the ChiMinimizer constructor.
The only argument from the command line is the file name you want
saved in the "chi" subdirectory. Text files will be created numbered
from 0 (to limit the size of each file) until all values are created.
"""
class ChiMinimizer:
    
    # TEMP: using spawn method due to memory issues
    # TEMP: outer radius fixed to 1.05*inner
    
    def __init__(self, chi_file):
        self.chi_file = 'chi/' + chi_file
        self.fits_file = 'chimin.fits'
        # inner radius
        self.innerRads = numpy.arange(50,70,2)#(40, 80, 5)
        # outer radius currently fixed to 1.05*inner
        self.outerRads = numpy.arange(135, 350, 25)
        # grain size
        temp = numpy.arange(-1.5,0,.1)#(-1, 1.5, .25)
        self.grainSizes = [10**x for x in temp]
        # disk mass
        temp = numpy.arange(-4.5,-3,.1)#(-4, -3, .25)
        self.diskMasses = [10**x for x in temp]
        # power law
        # TEMP: constant
        self.powerLaws = [1]
        # grain efficiency
        self.grainEfficiencies = numpy.arange(.2, 1.99, .1)#same
        # set visibility data
        self.image_width = 512
        self.inclination_angle = 84.3
        self.position_angle = 70.3
    
    def findBest(self):
        f = open(self.chi_file, 'w')
        def writeChi(*args):
            towrite = ''
            for arg in args[:-1]:
                towrite += str(arg) + ' '
            towrite += str(args[-1]) + '\n'
            f.write(towrite)
        
        startTime = time.time()
        numModels = 0
        # initialize our models
        disk = Disk(5, 10, 10, .002, .5, .5)
        vis = VisibilityGenerator(512, self.inclination_angle, self.position_angle, self.fits_file)
        for i in self.innerRads:
            """
            for o in self.outerRads:
                # doesn't make sense to have larger inner radius
                if i >= o:
                    continue
            """
            o = 1.05*i
            for g in self.grainSizes:
                print 'model:', str(numModels), 'radii:', str(i) + '-' + str(o), 'grainsize:', g
                for m in self.diskMasses:
                    for p in self.powerLaws:
                        # explicitly flush the data to file at this point
                        f.flush()
                        for e in self.grainEfficiencies:
                            numModels += 1
                            args = (i, o, g, m, p, e)
                            # compute chi-squared for the SED
                            disk.changeParameters(*args)
                            sed_chi = disk.computeChiSquared()
                            # compute chi-squared for the visibilities
                            vis_chi = vis.computeChiSquared(disk)
                            # combine the two separate chi values
                            chi_squared = sed_chi + vis_chi
                            # write to file
                            writeChi(i, o, g, m, p, e, sed_chi, vis_chi, chi_squared)
        f.close()
        # performance calculations
        endTime = time.time()
        compTimeSeconds = endTime-startTime
        # performance output
        print 'computation took', compTimeSeconds/3600, 'hours'
        print 'calculated', numModels, 'models'
        print 'took about', compTimeSeconds/numModels, 'seconds/model'

file_name = sys.argv[1]
if os.path.exists(file_name):
    ok = raw_input('OK to overwrite' + file_name + '? (y/n): ')
    if ok != 'y':
        exit()
    os.system('rm ' + file_name)
chi = ChiMinimizer(file_name)
chi.findBest()