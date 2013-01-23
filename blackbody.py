import math
import numpy
from decimal import *
import matplotlib.pyplot as plt

global h, k, c, b
h = Decimal(6.6260695729e-34)
k = Decimal(1.380648813e-23)
c = Decimal(2.99792458e8)
b = Decimal(2.8977685e-3)

def bbfunction(temp):
    def f(lam):
        numer = 2*h*(c**2)
        expo = h*c/(lam*k*temp)
        denom = (lam**4)*((Decimal(math.e)**expo)-1)
        return numer/denom
    return f

def blackbody(temp):
    maxlam = b/temp
    center = round(math.log(maxlam, 10))
    x = numpy.arange(center-2, center+5, .01)
    x = map(lambda u: Decimal(10**u), x)
    y = map(bbfunction(temp), x)
    
    x = map(lambda u: u*Decimal(1e6), x)
    plt.loglog(x, y)
    plt.xlabel('lambda (microns)')
    plt.ylabel('F*lambda')
    plt.ylim(ymin=1e-16)
    plt.show()

while True:
    temp = raw_input("Enter temperature: ")
    if temp == "exit":
        break
    temp = Decimal(temp)
    blackbody(temp)
    