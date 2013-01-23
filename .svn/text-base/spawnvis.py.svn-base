from disk import Disk
import visgen
from visgen import VisibilityGenerator
import sys

width = int(sys.argv[1])
angle = float(sys.argv[2])
rot = float(sys.argv[3])
filename = sys.argv[4]
disk_params = sys.argv[5:]
disk_params = map(float, disk_params)

disk = Disk(*disk_params)
f = open('spawn.txt', 'w')
vis = VisibilityGenerator(width, angle, rot, filename)
chi = vis.computeChiSquared(disk)
f.write(str(chi) + '\n')
f.close()