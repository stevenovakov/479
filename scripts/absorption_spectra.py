#
#
# Author:
#   Steve Novakov - Sept 2014
#

import scipy.constants as scicon
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.integrate import quad

def MaxwellBoltzmannOnedir(v, alpha):
  return math.e**(-1 * v**2 / (2 * alpha**2))
# remember to multiply by (2/scicon.pi)**0.5 / alpha**3

def LorentzAbsorb(dnu, nu0, vel, gamma):
  return 1 / ((dnu - nu0*vel/scicon.c)**2 + gamma**2/4.0)
# remember to multiply by gamma/(2*scicon.pi)

def IntegralFunction(vel, dnu, nu0, alpha, gamma):
  return LorentzAbsorb(dnu, nu0, vel, gamma)*MaxwellBoltzmannOnedir(vel, alpha)

rc('text', usetex=True)

dnu = 2.0 # GHz
N = 2000
nuspace = np.linspace(-1*dnu, dnu, N)

rb_nu0_D1 = 377107.463380 # GHz
rb_gamma_D1 = 0.0057500 # GHz
rb_nu0_D2 = 384230.484468 # GHz
rb_gamma_D2 = 0.0060666 # GHz

rb_isat = 2e-3 * 10000.0 # W/m**2

mass_rb87 = 86.909180520 * scicon.physical_constants["atomic mass unit-kilogram relationship"][0] # KG

T = 295.0 # Kelvin
alpha = (scicon.k * T / mass_rb87)**0.5
mean_velocity = 2 * alpha * (2/scicon.pi)**0.5
vmax = 10 * mean_velocity

print str(vmax) + "(m/s)?"

#NV = 1000
#vspace = np.linspace(0, vmax, NV)

absorption_response = []

constant_maxwell = 1/(alpha * (2 * scicon.pi)**0.5)
constant_lorentz = rb_gamma_D1/(2*scicon.pi)


#
# TODO Steve - this is a kluge - use more accurate integration methods
#
for n in nuspace:
  integral = 0 #LorentzAbsorb(n, rb_nu0_D1, rb_gamma_D1)
  integral = quad(IntegralFunction, -1*vmax, vmax, args=(n, rb_nu0_D1, alpha, rb_gamma_D1))[0]

  #absorption_response.append(math.e**(-1*integral*constant_factor))
  absorption_response.append(math.e**(-1*integral*constant_maxwell*constant_lorentz))

plt.plot(nuspace, absorption_response, linewidth=3, color="red")

plt.title("Doppler Broadened Absorption Feature for Rb87 D1 Transition")
plt.ylabel(r"Transmission Ratio : $ e^{- \tau(\nu)} $")
plt.xlabel(r"$ \nu - \nu_0 $ (GHz)")

plt.show()