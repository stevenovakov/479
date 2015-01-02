#
# ENPH 479
# University of British Columbia
#
# Creation Date:
#			2014-12-28
#
# Author(s):
#			Steve Novakov


#
# TODO Steve:
#   - twin axes for absorption spectrum (on same subplot)
#   - peak search fo true resonance feature
#   - centering around true resonance feature
#   - establish GHz/second scale
#   - make the display of the absorption spectrum and modxfer spectrum
#     nicer by having same vertical scale (will have to move up/down)
#

import os
import matplotlib.pyplot as pt
import numpy as np

# this is the data root, DROPBOX/ENPH\ 479/rubidium
# specify power/temp/old data directory later
data_dir = os.environ['DATA479']
root_dir = os.getcwd()

def OSDirAppend(target_dir):

  if 'nt' in os.name:
    return "\\" + target_dir
  else:
    return "/" + target_dir

def FindResonancesFFT():

  # DOES NOT WORK VERY WELL AT THE MOMENT, USE SG VERSION

  from scipy.fftpack import fft, ifft

  dump_dir = root_dir

  source_dir_85 = data_dir + OSDirAppend('power') + OSDirAppend('85')
  source_dir_87 = data_dir + OSDirAppend('power') + OSDirAppend('87')

  for filename in os.listdir(source_dir_85):

    file_path = source_dir_85 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    t_fit = file_data['t1']
    v_fit = file_data['v1']

    nv = len(v_fit)

    pad = nv*10

    v_filt = fft(v_fit, pad)

    cutoff = int(pad/250.)

    v_filt[cutoff:pad-cutoff] = np.zeros(pad-2*cutoff)

    v_filt = ifft(v_filt)

    v_filt = v_filt[0:nv]

    pt.plot(t_fit, v_fit, t_fit, v_filt)
    pt.show()

    if(raw_input("Enter c to close PROGRAM: ") == "c"):
      break

  return


def FindResonancesSG():

  from savitzkygolay import SavitzkyGolay

  def naive_peak_find(vdata, offset):

    peak_list = []

    for i in xrange(1, len(vdata)-1):
      deriv1 = vdata[i]-vdata[i-1]
      deriv2 = vdata[i+1]-vdata[i]

      # only care about "dips"
      if (deriv1 < 0.0) and (deriv2 > 0.0):
        peak_list.append(i+offset)
        continue

    return peak_list

  dump_dir = root_dir

  source_dir_85 = data_dir + OSDirAppend('power') + OSDirAppend('85')
  source_dir_87 = data_dir + OSDirAppend('power') + OSDirAppend('87')

  feature_separation_85 =

  for filename in os.listdir(source_dir_85)[0:1]:

    file_path = source_dir_85 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    t_fit = file_data['t1']
    v_fit = file_data['v1']

    # Piecewise fit using
    #     11 points per window (must be odd)
    #     9th order polynomial (must be > 1, < window size)

    num_passes = 50

    v_sgf = np.copy(v_fit)

    for i in xrange(0, num_passes):
      v_sgf = SavitzkyGolay(v_sgf, 61, 3)

    vn = len(v_fit)
    window = int(vn/4.)

    v_sgf2 = v_sgf[window:vn-window]

    # NOTE:
    #   noticed that sometimes there is a peak near the center, but
    #   it is inconsistent in appearence so only the last two peaks are used

    peaks = naive_peak_find(v_sgf2, window)[-2:]#

    # boolean, true if peaks are in front of center, as they are lower
    # actual frequency than the center of the spectrum
    flipped = (peaks[-1] > int(vn/2.))

    tpeaks = [t_fit[a] for a in peaks]

    # UNCOMMENT FOLLOWING TO VISUALLY INSPECT

    pt.plot(t_fit, v_fit, t_fit, v_sgf, linewidth=3)
    pt.vlines(tpeaks, min(v_fit), max(v_fit))
    pt.show()

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break

  return

def GeneratePowerFigures():
  '''generates images in rootdirectory/figures of power data'''

  dump_dir = root_dir + OSDirAppend('figures')

  source_dir_85 = data_dir + OSDirAppend('power') + OSDirAppend('85')
  source_dir_87 = data_dir + OSDirAppend('power') + OSDirAppend('87')

  for filename in os.listdir(source_dir_85):

    print "Plotting " + filename + "\n"

    file_path = source_dir_85 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1.plot(file_data['t1'], file_data['v1'], linewidth=2, color='#324259')

    ax1b = ax1.twinx()
    ax1b.plot(file_data['t2'], file_data['v2'], linewidth=2, color='#4646EA')

    ax1.get_xaxis().set_visible(False)

    ax2 = fig.add_subplot(212)
    ax2.plot(file_data['t3'], file_data['v3'], linewidth=2, color='#ff2b63')

    fig.show()

    if(raw_input("Enter c to close PROGRAM: ") == "c"):
      break

    pt.close(fig)


  return

FindResonancesSG()
#GeneratePowerFigures()






