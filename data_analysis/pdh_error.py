import scipy.constants as scicon
import math
import numpy as np
import matplotlib.pyplot as pt
from matplotlib import rc
from scipy.integrate import quad
import sys

rc('text', usetex=True)
pt.ion()

# F=3 -> F'=3 D2 transition (MHz)
rb_85_d2_fp3 = 384.230e6 - 1.265e3 + 100.205
# FWHM bandwidth around center (MHz)
gamma = 6.066

rb_molar = 85.46e-3 #kg
nrb = 7.4e9 * 1e4
# N/m**3 , density at 20C as a vapour

L = 1e-2

def n(w, w0, gamma):
  # r = (w0**2 - w**2) / ((w0**2-w**2)**2 - (gamma**2 * w**2))

  r = (w0 - w) / ((2*w)*((w0-w)**2 + (gamma/2)**2))

  r *= nrb * scicon.e**2 / (2 * scicon.m_e * scicon.epsilon_0)

  return 1 + r

def a(w, w0, gamma):
  #r = w**2 / ((w0**2 -w**2)**2 + (gamma**2 * w**2))

  r = (gamma/2)**2/((w0-w)**2 + (gamma/2)**2)

  #r *= scicon.N_A * scicon.e**2 / (2 * scicon.m_e * scicon.epsilon_0 * scicon.c)
  r *= nrb * scicon.e**2 / (gamma * scicon.m_e * scicon.epsilon_0 * scicon.c)

  return r

def trans(m, m0, g):
  w = 1e6 * m * 2 * math.pi
  w0 = 1e6 * m0 * 2 * math.pi
  gamma = 1e6 * g * 2 * math.pi
  return math.e**(-1*a(w, w0, gamma)*L)

def phase(m, m0, g):
  w = 1e6 * m * 2 * math.pi
  w0 = 1e6 * m0 * 2 * math.pi
  gamma = 1e6 * g * 2 * math.pi
  return (1 - n(w, w0, gamma))*(w/scicon.c)*L

#
# Generate
#     T
#     phase T
#     PDH in-phase
#     PDH quadrature
#

span = 6*gamma #MHz
fmod = 20 #MHz
N = 1e3

w = np.linspace(-1*(span+fmod), span+fmod, N)

t = []
p = []

for f in w:
  t.append(trans(f+rb_85_d2_fp3, rb_85_d2_fp3, gamma))
  p.append(phase(f+rb_85_d2_fp3, rb_85_d2_fp3, gamma))

tf = []

for am,ph in zip(t,p):
  tf.append(am*math.e**(1j*ph))

window = 0

select = w[0]

while(select - w[0] < fmod):
  window += 1
  select = w[window]

inphase = []
quadrature = []
wsub = w[window:-1*window]

for s in xrange(0, len(tf)-2*window):

  center = s + window
  plus = s + 2*window
  minus = s

  tfr = tf[center].real
  tfi = tf[center].imag
  tfpr = tf[plus].real
  tfpi = tf[plus].imag
  tfmr = tf[minus].real
  tfmi = tf[minus].imag

  inphase.append(tfr*(tfpr-tfmr) + tfi*(tfpi-tfmi))
  quadrature.append(tfi*(tfpr + tfmr) - tfr*(tfpi + tfmi))

inphase = np.array(inphase)
quadrature = np.array(quadrature)
wsub = np.array(wsub)
p = np.array(p)*1000
t = np.array(t)
#

# PLOTTING
#

fig = pt.figure()

ax1 = fig.add_subplot(311)
ax1.plot(w, t, linewidth=3, color="#000066")

ax2 = fig.add_subplot(312)
ax2.plot(w, p, linewidth=3, color="#003399")

ax3 = fig.add_subplot(313)
ax3.plot(wsub, inphase, linewidth=3, label="in phase", color="#0099FF")
ax3.plot(wsub, quadrature, linewidth=3, label="quadrature", color="#FF3300")

ylim1 = ax1.get_ylim()
ylim2 = ax2.get_ylim()
ylim3 = ax3.get_ylim()
scale3 = 1.1
ylim3 = [ylim3[0]*scale3, ylim3[1]*scale3]

ax1.vlines([0.0], ylim1[0], ylim1[1])
ax2.vlines([0.0], ylim2[0], ylim2[1])
ax3.vlines([0.0], ylim3[0], ylim3[1])

ax1.set_ylim(ylim1)
ax2.set_ylim(ylim2)
ax3.set_ylim(ylim3)

ax1.set_xlim([-1*span, span])
ax2.set_xlim([-1*span, span])
ax3.set_xlim([-1*span, span])

ax1.set_ylabel(r"$\tau(\omega)$ (normalized)")
ax2.set_ylabel(r"$\phi_T(\omega)$ (mRad)")
ax3.set_ylabel(r"$\epsilon$ (arbitrary)")
ax3.set_xlabel(r"$\Delta f = f - f_0$ (MHz)")
ax3.legend(loc="upper left",prop={'size':11})

ax1.grid(True, which='both')
ax1.set_xticklabels([])
ax2.grid(True, which='both')
ax2.set_xticklabels([])
ax3.grid(True, which='both')

fig.show()
raw_input("..")
pt.close(fig)