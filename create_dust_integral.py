import os
import re
import pyfits
import numpy
import math
from scipy import integrate
from scipy.interpolate import interp1d

'''
Here, we tabulate some integrals that we need for the disk model with astrosilicate grain opacities.
'''

# define constants
h_const = 6.6260695729e-34
c_const = 2.99792458e8
k_const = 1.380648813e-23
grain_density = 2.7e3
star_luminosity = 3.839e26 * 0.583

print 'Reading in SED data...'
# read and store data from the Kurucz-Lejeune model, which is used to compute absorbed energy
f = open('Model Spectrum.txt','r')
R_star = 0.84*6.955*1e8
dist = 35*3.09e16
KL_lambda = []
KL_flux = []
line = f.readline()
while line != '':
  # convert nanometers to meters
  KL_lambda.append(float(line[4:12])*1e-9)
  KL_flux.append(float(line[29:38]))*1e23
  line = f.readline()
f.close()
KL_lambda.append(1e4*1e-6)
KL_flux.append(2.914e-13*1e23/1e-2)
KL_flux_function = interp1d(KL_lambda, KL_flux, bounds_error=False, fill_value=0)

def calculateIncomingFlux(radius, lamma):
  return KL_flux_function(lamma)*4*math.pi*(R_star/radius)**2

def calculatePlanckFunction(temperature, lamma):
  numerator = 2*h_const*(c_const**2)
  exponent = h_const*c_const/(lamma*k_const*temperature)
  denominator = (lamma**5)*(math.e**exponent-1)
  return numerator/denominator

# create temperature fits file
temp = numpy.arange(0,1000,1)
hdu = pyfits.PrimaryHDU(temp)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_temperature.fits')

#create radius fits file
rad = [x*1.496e11 for x in numpy.arange(0,300,1)] #Sample with a step size of 1 AU
hdu = pyfits.PrimaryHDU(rad)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_radii.fits')

# start grain size fits file
grains = []

# start QB matrix, Q matrix, and QF matrix
QB_integral = []
bigQ = []
QF_integral = []

for dirname, dirnames, filenames in os.walk('./dust'):
  regex = re.compile('./dust/astrosil_(\d*.\d*)mic')
  if re.match(regex, dirname):
    grain_size = float(re.match(regex, dirname).groups()[0]) * 1e-6
    grains.append(grain_size)
    lambder = [float(x)*1e-6 for x in pyfits.open(dirname + '/lambda.fits')[0].data]
    kapper = [float(x)/10 for x in pyfits.open(dirname + '/kappa.fits')[0].data]
    albeder = [float(x) for x in pyfits.open(dirname + '/albedo.fits')[0].data]
    q = []
    for i in range(len(lambder)):
      q.append(float(4/3) * kapper[i] * grain_density * grain_size * (1-albeder[i]))
    bigQ.append(q)
    qFunction = interp1d(lambder, q, bounds_error=False, fill_value=0) 
    lammax = max(lambder)
    lammin = min(lambder)
    QB_integral_values = []
    for i in temp:
      QBval = integrate.quad(lambda l: qFunction(l)*calculatePlanckFunction(i, l), lammin, lammax)[0]
      QB_integral_values.append(QBval)
    QB_integral.append(integral_values)
    QF_integral_values = []
    for i in rad:
      QFval = integrate.quad(lambda l: qFunction(l)*calculateIncomingFlux(i, l), lammin, lammax)[0]
      QF_integral_values.append(QFval)
    QF_integral.append(QFval)
		
#Save data into separate FITS files.
np_lambda = numpy.array(lambder)
hdu = pyfits.PrimaryHDU(np_lambda)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_lambda.fits')

np_Q = numpy.array(bigQ)
hdu = pyfits.PrimaryHDU(np_Q)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_Q.fits')

np_result = numpy.array(QB_integral)
hdu = pyfits.PrimaryHDU(np_result)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_QBintegrals.fits')

np_result = numpy.array(QF_integral)
hdu = pyfits.PrimaryHDU(np_result)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_QFintegrals.fits')

np_grains = numpy.array(grains)
hdu = pyfits.PrimaryHDU(np_grains)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_grain_sizes.fits')