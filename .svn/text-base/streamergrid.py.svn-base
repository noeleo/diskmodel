from visgen_streamers import VisibilityGenerator
from disk import Disk
import numpy as np
import sys
import math
from scipy.special import gamma

#A 2D grid that fits the visibilities with 2 parameters:  Streamer Contribution and Total Flux. (IR used to be tested too.)

#Make File to Save
if (len(sys.argv) < 3):
    print 'You are not saving, fyi.  If you want to save, the file name is sys.argv[2].'
if (len(sys.argv) == 3):
    save = open(sys.argv[2],'w')

wid = 1024
ang = 84.3
rot = 70.3
disk_scale = 1
streamer_scale = 1
fits = sys.argv[1]

disk = Disk(28.75, 93.75, 0.00358, 7.865e-6, 0.5, 0.323, 2e-8)
visgen = VisibilityGenerator(disk, wid, ang, rot, disk_scale, streamer_scale, fits)

density = 50
'''
R_min = 40.0
R_max = 80.0
R_step = (R_max - R_min)/density
'''
F_min = 0.005 #These limits set because our mcmc ensemble suggests the total flux should be 7.1 \pm 0.3.  This gives a broad range around that.
F_max = 0.009
F_step = (F_max - F_min)/density

SC_min = 0.0  #SC for streamer contribution
SC_max = 1.0
SC_step = (SC_max - SC_min)/density

Bestmodel = [0,0,1e100] #Total Flux, Streamer Contribution, Chi^2
count = 0

for F_t in np.arange(F_min,F_max,F_step):
    #disk.changeParameters(R_t, R_t*1.05, 0.00358, 7.865e-6, 0.5, 0.323, 2e-8)
    streamer_flux = visgen.getStreamerFlux()
    disk_flux = disk.calculateFlux(0.00132949)/225.5e9 #Our wavelength in meters, divided by 225.5 because it's nu*F_nu 
    for SC in np.arange(SC_min,SC_max,SC_step):
        print "Beginning: F_t =", F_t, "and Streamer Contribution =", SC
        visgen.changeParameters(disk, 1, 1)
        streamer_scalar = SC*F_t/streamer_flux
        disk_scalar = (1 - SC)*F_t/disk_flux
        #streamer_scalar = F_t/(streamer_flux + disk_flux)
        visgen.changeParameters(disk, disk_scalar, streamer_scalar)
        chi_t = visgen.computeChiSquared(disk)
        print "Chi^2 =", chi_t
        
        if chi_t < Bestmodel[2]:
            Bestmodel[0] = F_t
            Bestmodel[1] = SC
            Bestmodel[2] = chi_t
        
        count += 1
        #print "Completed: F_t =", F_t, "and Streamer Contribution =", SC
        #print visgen.getTotalFlux()
        print "Finished", count, "of", density**2
        if (len(sys.argv) == 3):
            save.write(str(SC)+' '+str(F_t)+' '+str(chi_t)+'\n')

print "Bestmodel:", Bestmodel
if (len(sys.argv) == 3):
    save.close()

'''
Used continuity, varied Inner Radius and Total Flux

Bestmodel: [50.0, 0.0054000000000000003, 456089.59] when
density = 10

R_min = 50.0
R_max = 80.0
R_step = (R_max - R_min)/density

F_min = 0.004
F_max = 0.011
F_step = (F_max - F_min)/density

Second grid...
Bestmodel: [51.0, 0.0048999999999999998, 456088.62]
But when I use visfit, I get chi^2 = 456135.0...

Third grid
Bestmodel: [48.0, 0.0074999999999999963, 456088.81]

'''

'''
With continuity and fixed radius at 61.25 AU, just varying total flux,
Bestmodel: [0.0078499999999999855, 0, 456103.56]
'''

'''
Broke continuity, varied Total Flux and Streamer Contribution, fixed Inner Radius to 61.25 AU
Bestmodel: [0.0074999999999999963, 0.30000000000000004, 456084.69]
'''