from visgen_streamers import VisibilityGenerator
from disk import Disk
import sys

wid = 1024
ang = 84.3
rot = 70.3
disk_scale = 1
streamer_scale = 1
fits = sys.argv[1]
#disk_params = map(float, sys.argv[2:])

R_t = 61.25
#F_t = 0.00785
#SC = 0.3

disk = Disk(R_t, R_t*1.05, 0.00358, 7.865e-6, 0.5, 0.323, 2e-8)
visgen = VisibilityGenerator(disk, wid, ang, rot, disk_scale, streamer_scale, fits)
#streamer_flux = visgen.getStreamerFlux()
#disk_flux = disk.calculateFlux(0.00132949)/225.5e9*disk_scale #Our wavelength in meters, divided by 225.5 because it's nu*F_nu
#streamer_scalar2 = SC*F_t/streamer_flux
#disk_scalar2 = (1 - SC)*F_t/disk_flux
#flux_scalar = F_t/(streamer_flux + disk_flux)
#visgen.changeParameters(disk, flux_scalar, flux_scalar)


# this will call generateFits() and generateVisibility()
# which will in turn produce the files as an effect
print 'chi-squared:', visgen.computeChiSquared(disk)
print 'Actual Total Flux, Disk Flux =', visgen.getTotalFlux()
print 'Theoretical Disk Flux =', disk.calculateFlux(0.00132949)/225.5e9#*flux_scalar
print 'Theoretical Streamer Flux =', visgen.getStreamerFlux()

'''
To make the total flux with streamers as if there were no streamers,
disk_scale = streamer_scale = 0.427495539186 = 0.00658735982805/0.0154091896271, no streamers/unnormalized with streamers

To normalize,
Disk Flux = 1 if disk_scale = 151.805886744
Streamer Flux = 1 if streamer_scale = 113.35516812

These values apply to the parameters disk = Disk(70.5, 74.0, 0.00358, 7.865e-6, 0.5, 0.323, 2e-8)
'''