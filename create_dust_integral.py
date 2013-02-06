import os
import re
import pyfits
import numpy
import math
from scipy import integrate
from scipy.interpolate import interp1d

# define constants
h_const = 6.6260695729e-34
c_const = 2.99792458e8
k_const = 1.380648813e-23
grain_density = 2.7e3
star_luminosity = 3.839e26 * 0.583

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

# start grain size fits file
grains = []

# start result matrix
result = []

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
			q.append(
				float(4/3) * kapper[i] * grain_density * grain_size * (1-albeder[i])
				)
		qFunction = interp1d(lambder, q, bounds_error=False, fill_value=0) 
		lammax = max(lambder)
		lammin = min(lambder)
		integral_values = []
		for i in temp:
			val = integrate.quad(lambda l: qFunction(l)*calculatePlanckFunction(i, l), lammin, lammax)[0]
			integral_values.append(val)
		result.append(integral_values)

np_result = numpy.array(result)
hdu = pyfits.PrimaryHDU(np_result)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_integrals.fits')

np_grains = numpy.array(grains)
hdu = pyfits.PrimaryHDU(np_grains)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('./dust/compiled_grain_sizes.fits')