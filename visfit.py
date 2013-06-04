from visgen import VisibilityGenerator
from disk_as import Disk
import sys

wid = 1024
ang = 84.3
rot = 70.3
fits = sys.argv[1]
disk_params = map(float, sys.argv[2:])

disk = Disk(disk_params[0], disk_params[1], 10**disk_params[2], 10**disk_params[3], disk_params[4], disk_params[5], 10**disk_params[6])
visgen = VisibilityGenerator(wid, ang, rot, fits)
sed_chi = disk.computeChiSquared()
vis_chi = visgen.computeChiSquared(disk)
print 'Visibility chi-squared:', vis_chi
print 'SED chi-squared:', sed_chi
print 'Sum:', vis_chi + sed_chi