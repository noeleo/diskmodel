import math
import numpy
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
from scipy.special import gamma
import warnings
import pyfits
import time


class Disk:

    stefbolt_const = 5.67e-8
    h_const = 6.6260695729e-34
    k_const = 1.380648813e-23
    c_const = 2.99792458e8
    
    """
    inner and outer radius in units of AU
    grain size in units of microns
    disk mass in units of earth masses
    """
    def __init__(self, innerRadius, outerRadius, grainSize, diskMass, powerLaw, grainEfficiency, beltMass):
        # set the grain density to 2.7 grams/cm^3 in kilograms/m^3
        self.grainDensity = 2.7e3
        # set the star temperature to solar temperature in Kelvin
        self.starTemperature = 5500.0
        # set the star radius to solar radius in meters
        self.starRadius = 6.955e8*.84
        # set the star luminosity to solar luminosity in Watts
        self.starLuminosity = 3.839e26*.583
        # set the star distance to 35 parsecs in meters
        self.starDistance = 35*3.08568e16
        # disable overflow warnings
        warnings.filterwarnings("ignore", "overflow encountered in double_scalars")
        # convert the radii to meters
        self.innerRadius = float(innerRadius)*1.496e11
        self.outerRadius = float(outerRadius)*1.496e11
        # convert grain size to meters
        self.grainSize = float(grainSize)/1e6
        # convert to kilograms
        self.diskMass = float(diskMass)*5.9742e24
        self.powerLaw = float(powerLaw)
        self.grainEfficiency = float(grainEfficiency)
        self.beltMass = float(beltMass)*5.9742e24
        # calculate zeta(grainEfficency + 4)
        argument = self.grainEfficiency + 4
        self.zeta = 0
        for n in range(1, 100):
            self.zeta += 1.0/(n**argument)
        # calculate the surface density at 100 AU (sigma_100)
        j = 2 - self.powerLaw
        m = (self.outerRadius**j) - (self.innerRadius**j)
        self.surfaceSigma = self.diskMass*j/(2*math.pi*((100*1.496e11)**self.powerLaw)*m)
        self.grainMass = self.grainDensity*(4.0/3.0*math.pi*(self.grainSize**3))
        
        # read and store data from the SED model
        f = open('Model Spectrum.txt','r')
        radius = 0.84*6.955*1e8
        dist = 35*3.09e16
        self.data_lambda = []
        self.data_flux = []
        line = f.readline()
        while line != '':
            # convert nanometers to meters
            lamma = float(line[4:12])*1e-9
            flux = float(line[29:38])*4*3.14159*(radius/dist)**2*1e23
            # multiply by nu
            flux *= self.c_const/lamma
            self.data_lambda.append(lamma)
            self.data_flux.append(flux)
            line = f.readline()
        f.close()
        self.data_lambda.append(1e4*1e-6)
        self.data_flux.append(2.914e-13*4*3.14159*(radius/dist)**2*1e23*self.c_const/1e-2)
        self.data_lambda = self.convertToMicrons(self.data_lambda)
        
        # read and store observed fluxes
        f = open('Observed Fluxes.txt', 'r')
        self.sample_lambda = []
        self.sample_flux = []
        self.sample_error = []
        line = f.readline()
        while line != '':
            # read wavelength, flux, and error
            info = line.split()
            # convert microns to meters, but use microns
            lam = float(info[0])*1e-6
            self.sample_lambda.append(float(info[0]))
            # make flux in Janksys into Jansky*Hz by muliplying with nu
            flu = float(info[1])*self.c_const/lam
            self.sample_flux.append(flu)
            # do the same thing with error
            err = float(info[2])*self.c_const/lam
            self.sample_error.append(err)
            line = f.readline()
        f.close()
        
        #read and store IRAS fluxes, which will only be used for display purposes, not fitting
        f = open('IRS Spectrum.txt', 'r')
        self.IRS_lambda = []
        self.IRS_flux = []
        self.IRS_error = []
        line = f.readline()
        while line != '':
            info = line.split()
            lam = float(info[0])
            self.IRS_lambda.append(lam)
            flu = float(info[1])*self.c_const/(lam*1e-6) #Want nu*F_nu
            self.IRS_flux.append(flu)
            err = float(info[2])/flu*self.c_const/(lam*1e-6)*100
            self.IRS_error.append(err)
            line = f.readline()
        self.ast_avg_err = numpy.mean(self.IRS_error)
        #print "Average % error in IRS Spectrum =",self.ast_avg_err,"%" #Turns out it's 5.89%
        f.close()
        
        #Read in Q*B table and associated arrays.
        compiled_temp = [float(x) for x in pyfits.open('./dust/compiled_temperature.fits')[0].data]
        compiled_grain_sizes = [float(x) for x in pyfits.open('./dust/compiled_grain_sizes.fits')[0].data]
        compiled_integrals = pyfits.open('./dust/compiled_integrals.fits')[0].data
        compiled_q = pyfits.open('./dust/compiled_Q.fits')[0].data
        self.sorted_q = numpy.array([y for (x,y) in sorted(zip(compiled_grain_sizes,compiled_q))])
        self.sorted_lambda = numpy.array([float(x) for x in pyfits.open('./dust/compiled_lambda.fits')[0].data])
        self.sorted_integrals = numpy.array([y for (x,y) in sorted(zip(compiled_grain_sizes,compiled_integrals))])
        self.sorted_grain_sizes = numpy.array(sorted(compiled_grain_sizes))
        self.sorted_temp = numpy.array(compiled_temp)
        #The following need to be put in calculateGrainTemperature when we do a distribution
        self.grain_close = min(self.sorted_grain_sizes, key=lambda y: math.fabs(y-self.grainSize))
        self.grain_index = numpy.where(self.sorted_grain_sizes==self.grain_close)[0][0]
        self.integral_list = [self.sorted_integrals[self.grain_index][x] for x in range(len(self.sorted_temp))]

        # generate interpolation function
        loglamb = map (math.log10, self.data_lambda)
        logflux = map(math.log10, self.data_flux)
        # if out of bounds, interpolates to 0
        logFunct = interp1d(loglamb, logflux, bounds_error=False, fill_value=0)
        self.interpol_funct = lambda x: 10**(logFunct(math.log10(x)))
    
        #Add in an asteroid belt with a fixed temperature and mass
        self.asteroid_radius = 1.0e-6 #arbitrary
        self.M_aster = 4/3*math.pi*self.asteroid_radius**3*self.grainDensity 
        self.n_asteroids = self.beltMass/self.M_aster
        self.Temp_a = 100

        self.grain_temps = []
        for i in range(1000):
            lhs = self.starLuminosity/(16*(math.pi**2)*(radius**2))
            integral_close = min(self.integral_list, key=lambda y: math.fabs(y-lhs))
            integral_index = numpy.where(self.integral_list==integral_close)[0][0]
            temperature = self.sorted_temp[integral_index]
            self.grain_temps.append(temperature)
        
    """
    changes the paramters to the disk
    """
    def changeParameters(self, innerRadius, outerRadius, grainSize, diskMass, powerLaw, grainEfficiency, beltMass):        
        # convert the radii to meters
        self.innerRadius = float(innerRadius)*1.496e11
        self.outerRadius = float(outerRadius)*1.496e11
        # convert grain size to meters
        self.grainSize = float(grainSize)/1e6
        # convert to kilograms
        self.diskMass = float(diskMass)*5.9742e24
        self.powerLaw = float(powerLaw)
        self.grainEfficiency = float(grainEfficiency)
        self.beltMass = float(beltMass)*5.9742e24
        self.n_asteroids = self.beltMass/self.M_aster
        # calculate zeta(grainEfficency + 4)
        argument = self.grainEfficiency + 4
        self.zeta = 0
        for n in range(1, 100):
            self.zeta += 1.0/(n**argument)
        # calculate the surface density at 100 AU (sigma_100)
        j = 2 - self.powerLaw
        m = (self.outerRadius**j) - (self.innerRadius**j)
        self.surfaceSigma = self.diskMass*j/(2*math.pi*((100*1.496e11)**self.powerLaw)*m)
        self.grainMass = self.grainDensity*(4.0/3.0*math.pi*(self.grainSize**3))
    
    """
    gets the inner and outer radii in AU
    """
    def getOuterRadius(self):
        return self.outerRadius/1.496e11
    
    def getInnerRadius(self):
        return self.innerRadius/1.496e11
    
    """
    gets the star distance in parsecs
    """
    def getStarDistance(self):
        return self.starDistance/3.08568e16
    
    """
    converts a list of meters into microns
    """
    def convertToMicrons(self, lst):
        return [x*1e6 for x in lst]
    
    # TODO: make sure this is the same as calculate flux
    """
    Takes a radius (in AU) and frequency (in GHz)
    Returns the point flux at that radius and frequency in Jansky's
    """
    def calculatePointFlux(self, radius, frequency):
        radius = float(radius)*1.496e11
        # make sure that we are inside of the ring
        if radius > self.outerRadius or radius < self.innerRadius:
            return 0
        # convert frequency to GHz
        lamma = self.c_const/(frequency*1e9)
        Q = self.qFunction(lamma)
        flux = 1e26*Q*(self.grainSize**2)*self.calculateGrainDistribution(radius)*self.calculateGrainBlackbody(radius, lamma)/(self.starDistance**2)
        return flux
    
    """
    input lambda wavelength in meters
    return nu*B_nu(lambda) in Jansky*Hz
    """
    def calculateFlux(self, lamma):
        # integrate returns a list of integral value and error, we only want value
        fluxIntegral = integrate.quad(lambda radius: radius*self.calculateGrainDistribution(radius)
                                                 *self.calculateGrainBlackbody(radius, lamma), self.innerRadius, self.outerRadius,
                                                 limit=5)[0]
        # scale by nu
        nu = self.c_const/lamma
        Q = self.qFunction(lamma)
        flux = Q*nu*2*math.pi*fluxIntegral*1e26*(self.grainSize**2)/(self.starDistance**2)
        return flux
    
    """
    returns the surface density at a given radius
    """
    def calculateGrainDistribution(self, radius):
        surfaceDensity = self.surfaceSigma*((radius/(100*1.496e11))**(-self.powerLaw))
        return surfaceDensity/self.grainMass
    
    """
    returns B_nu(T) in Janskys ( = 10^-26 * Watts / m^2 / Hz )
    """
    def calculateGrainBlackbody(self, radius, lamma):
        try:
            nu = self.c_const/lamma
            exponent = self.h_const*nu/(self.k_const*self.calculateGrainTemperature(radius))
            numerator = 2*self.h_const*(nu**3)*math.pi 
            denominator = (self.c_const**2)*(math.e**exponent - 1)
            grainBlackbody = numerator/denominator
        except OverflowError:
            return 0
        return grainBlackbody
    
    def calculatePlanckFunction(self, temperature, lamma):
        numerator = 2*self.h_const*(self.c_const**2)
        exponent = self.h_const*self.c_const/(lamma*self.k_const*temperature)
        denominator = (lamma**5)*(math.e**exponent-1)
        return numerator/denominator

    """
    Approximates the temperature of a grain using a precalculated table of integrals.
    """
    def calculateGrainTemperature(self, radius):
        return self.grain_temps[int(round(radius))]
        '''
        lhs = self.starLuminosity/(16*(math.pi**2)*(radius**2))
        
        integral_close = min(self.integral_list, key=lambda y: math.fabs(y-lhs))
        integral_index = numpy.where(self.integral_list==integral_close)[0][0]
        
        temperature = self.sorted_temp[integral_index]
        '''
        '''
        plt.plot(self.sorted_temp,integral_list)
        print "grain size =", self.grainSize
        print "approx grain size =", grain_close
        print "T =", temperature
        print "lhs =", lhs
        print "Integral =", integral_close
        '''
        return temperature
    
    """
    Returns the emissivity of a grain at a given wavelength from a lookup table.
    """
    def qFunction(self, lamma):
        lamma_close = min(self.sorted_lambda, key=lambda y: math.fabs(y-lamma))
        lamma_index = numpy.where(self.sorted_lambda==lamma_close)[0][0]
        #print lamma_close, lamma, lamma_index
        return self.sorted_q[self.grain_index][lamma_index]
    
    '''
    def calculateGrainTemperature(self, radius):
        #x = time.time()
        lhs = self.starLuminosity/(16*(math.pi**2)*(radius**2))
        # initialize temperature for binary search
        start = 0.0
        end = 1000.0
        # make threshold 0.1% of desired value
        threshold = 0.001*lhs
        while True:
            temp = (end+start)/2 #I don't like this.  Let's implement Newton's method later.
            print "Temperature =", temp
            #rhs = integrate.quad(lambda l: self.qFunction(l)*self.calculatePlanckFunction(temp, l), self.lammin, self.lammax)[0]
            rhs = RectBivariateSpline(self.sorted_grain_sizes,self.sorted_temp,self.sorted_integrals)
            print "Right =", rhs, "Left =", lhs
            diff = rhs-lhs
            if math.fabs(diff) < threshold:
                #y = time.time()
                print 'took ' + str(y-x) + ' seconds'
                return temp
            elif diff > 0:
                # decrease t
                end = temp
            else:
                # increase t
                start = temp
        '''
    
    """
    returns nu*B_nu(lambda) in Jansky*Hz of the host star
    """
    def calculateStarBlackbody(self, lamma):
        nu = self.c_const/lamma
        exponent = self.h_const*nu/(self.k_const*self.starTemperature)
        numerator = 2*self.h_const*(nu**3)
        denominator = (self.c_const**2)*(math.exp(exponent)-1)
        nu = self.c_const/lamma
        # not sure what to scale by (4pi?)
        return numerator/denominator*nu*1e26*(self.starRadius**2)/(self.starDistance**2)
    
    """
    generates lambda and nu*B_nu values in meters and Janksy*Hz, respectively
    """
    def generateModel(self):
        self.generateAsteroids()
        # sample from .1 microns to 10^4 microns
        x = numpy.arange(-7, -2, .01)
        x = [10**power for power in x]
        y = [self.calculateFlux(lamma) for lamma in x]
        self.disk_lambda = self.convertToMicrons(x)
        self.model_lambda = self.convertToMicrons(x)
        self.disk_flux = y
        z = self.asteroid_flux
        self.model_flux = map(lambda lamma,flux1,flux2: flux1+flux2+self.interpol_funct(lamma*1e6), x, y, z)
    
    def calculateAsteroidBelt(self, lamma):
        try:
            nu = self.c_const/lamma
            exponent = self.h_const*nu/(self.k_const*self.Temp_a)
            numerator = 2*self.h_const*(nu**3)*math.pi*self.n_asteroids #Multiply by number of asteroids
            denominator = (self.c_const**2)*(math.e**exponent - 1)
            asteroidBlackbody = numerator/denominator*nu*1e26*self.asteroid_radius**2/self.starDistance**2
        except OverflowError:
            return 0
        return asteroidBlackbody
    
    def generateAsteroids(self):
        x = numpy.arange(-7, -2, 0.01)
        x = [10**power for power in x]
        y = [self.calculateAsteroidBelt(lamma) for lamma in x]
        self.asteroid_lambda = self.convertToMicrons(x)
        self.asteroid_flux = y
        #print 'Warm Belt Flux at 1.3mm = ', self.calculateAsteroidBelt(1.3e-3)*1.3e-3/self.c_const, 'Jy'

    def generateInterpolation(self):
        # generate interpolated data for the actual data we have
        def calculateInterpol(lam):
            # for less than 10^-6, ignore disk because so faint
            if lam*1e-6 > 1e-6:
                return self.calculateFlux(lam*1e-6)+self.interpol_funct(lam)+self.calculateAsteroidBelt(lam*1e-6)
            else:
                return self.interpol_funct(lam)
        self.interpol_flux = map(calculateInterpol, self.sample_lambda)
    
    """
    plots the SED of the star and disk
    """
    def plotSED(self):
        self.generateModel()
        self.generateInterpolation()
        
        # plot the observed data
        plt.errorbar(self.sample_lambda, self.Lsun(self.sample_flux), yerr=self.Lsun(self.sample_error), fmt='o', label = 'Observed Data', color='k')
        plt.loglog(self.IRS_lambda, self.Lsun(self.IRS_flux), label = 'IRS Spectrum', linewidth=2, color='y')

        # plot the disk model
        plt.loglog(self.model_lambda, self.Lsun(self.model_flux), '-', label="Best Fit Model", linewidth=2, color='b')
        plt.loglog(self.disk_lambda, self.Lsun(self.disk_flux), '--', label="Disk Model", linewidth=2, color='g')
        #plt.loglog(self.sample_lambda, self.interpol_flux, "o", label = 'Interpolated Model Data')
        plt.loglog(self.asteroid_lambda, self.Lsun(self.asteroid_flux), ':', label='Warm Dust Belt', linewidth=2, color='r')
        plt.errorbar([1.3e-3*1e6], self.Lsun([7.1e-3*self.c_const/1.3e-3]), self.Lsun([1.5e-3*self.c_const/1.3e-3]), fmt='s', color='m') #This is our data point, with 20% systematic uncertainty added in quadrature.
        
        # format and display the plot 
        plt.tick_params(labelsize=16)
        plt.xlabel('$\lambda$ $(\mu m)$', fontsize=20)
        plt.ylabel(r'$L_{\nu}$ $(L_{\odot})$', fontsize=20)
        plt.legend()
        plt.text(0.12, 0.60, 'HD 61005', fontsize=24)
        plt.xlim(1e-1,3e3)
        plt.ylim(1e-7,2e0)
        #plt.savefig('SED_mcmcbelt0704.eps')
        plt.show()
    
    """
    computes the chi-squared value of the disk and model data
    """
    def computeChiSquared(self):
        self.generateInterpolation()
        lamma = self.sample_lambda
        model_flux = self.interpol_flux
        actual_flux = self.sample_flux
        error = self.sample_error
        chi_squared = 0
        numToFit = 0
        for i in range(len(lamma)):
            # the others are very close and we don't want them to mess up the chi-squared
            if lamma[i] > 1e1 and lamma[i] < 1e3:
                numToFit += 1
                #print 'model:', model_flux[i]
                #print 'actual:', actual_flux[i]
                chi_squared += ((model_flux[i]-actual_flux[i])/error[i])**2
        # degrees of freedom = [number of samples] - [number parameters to fit]
        # TEMPO: parameters 6 normally
        dof = numToFit - 5
        chi_squared = chi_squared#/dof
        return chi_squared
    
    '''
    Converts fluxes that are in nu*F_nu into Solar Luminosity units
    '''
    def Lsun(self, fluxes):
        SolLum = [4*math.pi*self.starDistance**2*i/3.839e26/1e26 for i in fluxes]
        return SolLum
    
    '''
    Use Beckwith's formula to get a realistic disk mass
    '''
    def diskMassestimate(self):
        nu = 235.5e9
        kappa_nu = 0.1*(nu/1e12)**self.grainEfficiency
        Mass = 0.1*(0.02/kappa_nu) #Solar Masses
        return Mass/3.00266478958e-06 #Conversion to Earth Masses
    
    '''
    Gets the approximate temperature of the disk
    '''
    def disktemp(self):
        x = numpy.arange(-7, -2, .01)
        x = [10**power for power in x]
        y = [self.calculateFlux(lamma) for lamma in x]
        max = [0,0]
        for i in range(0,len(y)):
            if y[i] > max[1]:
                max[0] = x[i]
                max[1] = y[i]
        lambda_peak = max[0]
        T = 2.898e-3/lambda_peak
        return T
                
        
    '''
    Need to debug?  Look no further!  Well...there's nothing further anyway, but yeah, use this.
    '''
    def testplot(self):
        x = numpy.arange(-10, 5, .001)
        x = [10**power for power in x]
        y = [self.calculateGrainBlackbody(70.,lamma) for lamma in x]
        plt.loglog(x, y)
        plt.show()
