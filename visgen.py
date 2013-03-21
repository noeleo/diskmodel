import math
import numpy
import pyfits
import os
from disk_as import Disk

# TEMPO: NONREDUCED CHI

"""
For the generator to work, have folders named "fits" and "images" in
the directory from which you are calling this script. There should
also be a folder named "hd61005" in which is stored the visibility
file "hd61005_all.vis". From there, just run with command line arguments:
number of pixels across, inclination angle, rotation angle, and name of
the fits file you would like created in the fits directory. Everything
in the image directory is temporary and will be deleted on the next run.
The "images" directory holds the generated visibilities for the model
and the data, both with file extension "uvf".
"""
class VisibilityGenerator:
    
    """
    creates a visibility generator with width pixels,
    angle inclination and rotation rotation (in degrees)
    """
    def __init__(self, width, angle, rotation, fits_file):
        self.width = width
        # convert to radians
        self.incline = angle*math.pi/180.0
        self.rotation = (rotation+90)*math.pi/180.0
        
        # name files we need
        self.fits_name = fits_file
        self.fits_model_image = 'images/' + fits_file
        self.mir_model_image = 'images/model.mp'
        self.mir_model_vis = 'images/model.vis'
        self.fits_model_vis = 'images/model.uvf'
        # get the data files we need
        self.mir_data_vis = 'hd61005/hd61005_all.ppm.vis'
        self.fits_data_vis = 'hd61005/hd61005_all.ppm.uvf'
        
        # get the data from the observed visibility
        if not os.path.exists(self.fits_data_vis):
            command('fits in=' + self.mir_data_vis + ' op=uvout out=' + self.fits_data_vis)
        # store data for computing the chi squared later on
        observed = pyfits.getdata(self.fits_data_vis, 'PRIMARY').data[:,0,0,0,0]
        self.data_vis_real = observed[:,0]
        self.data_vis_imag = observed[:,1]
        self.data_vis_weight = observed[:,2]
    
    """
    generate a FITS file with name fits_file
    """
    def generateFITS(self, disk):
        # remove all current files which we will be working with
        os.system('rm -Rf images/*')
        # remove if the fits file already exists
        if os.path.exists('fits/' + self.fits_name):
            os.system('rm fits/' + self.fits_name)
            #print 'fits/' + self.fits_name, 'was overwritten.'
        
        # square images
        center = self.width/2.0
        
        # in units of AU
        outerRadius = disk.getOuterRadius()
        
        # now find the scaling factor
        pxwid = 3.0*outerRadius/self.width
        
        total_flux = 0
        
        # x and y are now transformed into new coordinate plane
        # center is at (0,0) with points sampled at center of each pixel
        
        # initialize an empty 2-D "image" of 0's
        image = [[0]*self.width for x in range(self.width)]
        # exploit the diagonal symmetry for a sub-2x speedup
        for j in xrange(self.width/2):
            y_face = (j-center+0.5)*pxwid

            for i in xrange(self.width):
                x_face = (i-center+0.5)*pxwid
                
                # first rotate, then transform for inclination angle
                x = x_face*math.cos(self.rotation) + y_face*math.sin(self.rotation)
                y = (-x_face*math.sin(self.rotation) + y_face*math.cos(self.rotation)) / math.cos(self.incline)
                
                radius = math.sqrt(x**2+y**2)
                # have flux in Jy/m^2, now need in Jy/pixel
                flux = disk.calculatePointFlux(radius, 225.5)/math.cos(self.incline)*((pxwid*1.496e11)**2)
                if flux:
                    # only assign if the flux is NOT 0
                    image[j][i] = flux
                    image[-j-1][-i-1] = flux
                    total_flux += flux*2
        # convert the Python image array into a numpy 4-D array for MIRIAD compatibility
        #print "Total Flux =", total_flux
        image = numpy.array([[image]])
        
        # calculate data necessary for header
        stardist = disk.getStarDistance()
        pix_parsec = pxwid/stardist/3600.0 # AU to degrees
        
        # create hdu object to encapsulate data with the header
        hdu = pyfits.PrimaryHDU(image)
        head = hdu.header
        
        head.update('OBJECT', 'HD-61005')
        head.update('SIMPLE', 'T')
        head.update('NAXIS', 4)
        head.update('NAXIS1', self.width)
        head.update('NAXIS2', self.width)
        head.update('NAXIS3', 1)
        head.update('NAXIS4', 1)
        head.update('CDELT1', -1.0*pix_parsec)
        head.update('CRPIX1', center+0.5)
        #head.update('CRVAL1', (7 + 35/60.0 + 47.4407/3600.0)*15)  #Old.  Not sure why this is different from others.
        #head.update('CRVAL1', (7 + 35/60.0 + 47.46192/3600.0)*15)  #J2000
        head.update('CRVAL1', (7 + 35/60.0 + 47.420/3600.0)*15) #Corrected for PM.  I'm lying to MIRIAD because the header in the ppm corrected data file is wrong.
        head.update('CTYPE1', 'RA---SIN')
        head.update('CDELT2', pix_parsec)
        head.update('CRPIX2', center+0.5)
        #head.update('CRVAL2', -1.0*(32 + 12/60.0 + 13.30/3600.0))  #Old
        #head.update('CRVAL2', -1.0*(32 + 12/60.0 + 14.0431/3600.0))  #J2000
        head.update('CRVAL2', -1.0*(32 + 12/60.0 + 13.32/3600.0)) #Corrected for PM
        head.update('CTYPE2', 'DEC--SIN')
        head.update('EPOCH', 2000)
        head.update('BSCALE', 1.0)
        head.update('BZERO', 0.0)
        head.update('CDELT3', 8.0e11)
        head.update('CRVAL3', -1.0)
        head.update('CRPIX3', 1.0)
        head.update('CTYPE3', 'STOKES')
        head.update('CDELT4', 4.0e9)
        head.update('CRPIX4', 1.0)
        head.update('CRVAL4', 225.538e9)
        head.update('CTYPE4', 'FREQ-LSR')
        head.update('BUNIT', 'JY/PIXEL')
        head.update('BTYPE', 'INTENSITY')

        # write to file
        hdu.writeto(self.fits_model_image)
        
        """
        print 'AU^2/pixel:', str(round(pxwid**2,4))
        print 'total flux (Jy):', total_flux
        print 'total flux (Jy*MHz):', str(round(total_flux*225.5e9/1e6,2))
        print 'actual SED flux at 225.5GHz:', str(round(disk.calculateFlux(3.0e8/225.5e9)/1e6,2))
        """
        
    """
    generate the visibilities of the model using MIRIAD routines
    """
    def generateVisibility(self):
        # make sure we have the files we need
        if not os.path.exists(self.mir_data_vis):
            print self.mir_data_vis, "visibility file does not exist."
            exit()
        
        # convert FITS model image into a Miriad image
        command('fits in=' + self.fits_model_image + ' op=xyin out=' + self.mir_model_image)
        # also copy this into the fits directory just for records
        command('cp -R ' + self.fits_model_image + ' fits/' + self.fits_name)
        # make the Miriad image into a visibility
        command('uvmodel model=' + self.mir_model_image + ' vis=' + self.mir_data_vis
                  + ' options=replace out=' + self.mir_model_vis)
        # convert the Miriad visibilities into FITS visibilities
        command('fits in=' + self.mir_model_vis + ' op=uvout out=' + self.fits_model_vis)
    
    def computeChiSquared(self, disk):
        # do all the prereqs
        self.generateFITS(disk)
        self.generateVisibility()
        # get the data from created files
        # model = pyfits.getdata(self.fits_model_vis, 'PRIMARY').data
        model_hdu = pyfits.open(self.fits_model_vis)
        model = model_hdu[0].data.field(5)
        # compute the chi squared
        # = sum of weight*[(model-data)**2] for real and imaginary
        numVis = len(model)
        # use numpy vector operations
        real = model[:,0,0,0,0,0]
        imag = model[:,0,0,0,0,1]
        chi = float((((real-self.data_vis_real)**2)*self.data_vis_weight).sum())
        chi += float((((imag-self.data_vis_imag)**2)*self.data_vis_weight).sum())
        # degrees of freedom = 2*[number of visibilities] - [number parameters to fit]
        # TEMPO: 4 parameters instead of 6 to fit
        #dof = 2*numVis-4
        chi = chi#/dof
        model_hdu.close(closed=True)
        return chi

"""
silences output to the command line
"""
def command(cmd):
    os.system(cmd + ' > /dev/null 2>&1')

"""
vis = VisibilityGenerator(512, 84.3, 70.3, 'whatev.fits')
chi = vis.computeChiSquared(Disk(66.167, 116.439, 2.803, 0.00108, 0.9917, 3.4478))
print chi
"""
"""
wid = int(sys.argv[1])
ang = float(sys.argv[2])
rot = float(sys.argv[3])
fits = sys.argv[4]

start = time.time()
visgen = VisibilityGenerator(wid, ang, rot, fits)
visgen.generateFITS(Disk(60, 100, 10, .002, .5, .5))
visgen.generateVisibility()
end = time.time()
print 'total time:', str(round(end-start,3)), 'sec'
"""