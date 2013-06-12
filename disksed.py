from disk_as import Disk
import sys

if (len(sys.argv) == 1):
    innerRad = raw_input('inner radius (AU): ')
    outerRad = raw_input('outer radius (AU): ')
    grainSize = raw_input('grain size (microns): ')
    diskMass = raw_input('disk mass (earth masses): ')
    powerLaw = raw_input('power law paramter (0-1): ')
    grainEfficiency = raw_input('grain efficiency parameter (0-1): ')
    beltMass = raw_input('belt mass (earth masses): ')
else:
    params = sys.argv[1:8]
disk = Disk(*params)
#print 'Calculating chi-squared...'
print 'Chi-squared:', disk.computeChiSquared()
print 'Now plotting...'
disk.plotSED()