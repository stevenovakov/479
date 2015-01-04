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
from savitzkygolay import SavitzkyGolay

# this is the data root, DROPBOX/ENPH\ 479/rubidium
# specify power/temp/old data directory later
data_dir = os.environ['DATA479']
root_dir = os.getcwd()

# peak distances calculated using Dan Steck Rb 87, 85 data

co13_co23_87 = 156.947070 / 2.0
co34_co23_85 = 63.40161 / 2.0

tzero = 0.0
mhzps = 1.0
tdelta = 1.0

source_dir_85 = data_dir + OSDirAppend('power') + OSDirAppend('85')
source_dir_87 = data_dir + OSDirAppend('power') + OSDirAppend('87')

write_params_85 = root_dir + OSDirAppend('85params.csv')
write_params_87 = root_dir + OSDirAppend('87params.csv')

def OSDirAppend(target_dir):

  if 'nt' in os.name:
    return "\\" + target_dir
  else:
    return "/" + target_dir

#
# for plotting, each file needs the following data:
#
# tzero = time of second peak (outward most), this is the origin/offset
# mhzps = MHz per second resolution of the plot
# isflip = is the plot flipped (goes down in frequency from L -> R)
#
# These "FindRes" functions find this write it to a file in the data dir
#


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

def FindRes85Power():

  wfile = open(write_params_85, 'wb')

  # RB85 VARIABLE POWER DATA

  for filename in os.listdir(source_dir_85):

    print "Processing " + filename

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

    if flipped:
      tzero = tpeaks[-1]
      tdelta = tpeaks[-1] - tpeaks[-2]
    else:
      tzero = tpeaks[0]
      tdelta = tpeaks[1] - tpeaks[0]

    mhzps = co34_co23_85 / abs(tdelta)

    wfile.write(",".join([filename, str(tzero), str(mhzps), str(flipped)])+"\n")

    # UNCOMMENT FOLLOWING TO VISUALLY INSPECT

    # pt.plot(t_fit, v_fit, t_fit, v_sgf, linewidth=3)
    # pt.vlines(tpeaks, min(v_fit), max(v_fit))
    # pt.show()

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break


  wfile.close()

  return

def FindRes87Power():

  wfile = open(write_params_87, 'wb')

  for filename in os.listdir(source_dir_87):

    print "Processing " + filename

    file_path = source_dir_87 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    t_fit = file_data['t1']
    v_fit = file_data['v1']

    # Piecewise fit using
    #     11 points per window (must be odd)
    #     9th order polynomial (must be > 1, < window size)

    num_passes = 30

    v_sgf = np.copy(v_fit)

    for i in xrange(0, num_passes):
      v_sgf = SavitzkyGolay(v_sgf, 51, 5)

    vn = len(v_fit)
    window = int(vn*0.202)

    v_sgf2 = v_sgf[window:vn-window]

    # NOTE:
    #   noticed that sometimes there is a peak near the center, but
    #   it is inconsistent in appearence so only the last two peaks are used

    peaks = naive_peak_find(v_sgf2, window)[-2:]#

    # boolean, true if peaks are in front of center, as they are lower
    # actual frequency than the center of the spectrum
    flipped = (peaks[-1] > int(vn/3.))

    tpeaks = [t_fit[a] for a in peaks]

    if flipped:
      tzero = tpeaks[-2]
      tdelta = tpeaks[-1] - tpeaks[-2]
    else:
      tzero = tpeaks[1]
      tdelta = tpeaks[1] - tpeaks[0]

    mhzps = co13_co23_87 / abs(tdelta)

    wfile.write(",".join([filename, str(tzero), str(mhzps), str(flipped)])+"\n")

    # UNCOMMENT FOLLOWING TO VISUALLY INSPECT

    # pt.plot(t_fit, v_fit, t_fit, v_sgf, linewidth=3)
    # pt.vlines(tpeaks, min(v_fit), max(v_fit))
    # pt.show()

  wfile.close()

  return


#
#
# The next set of functions deals with plotting figures
#
#

def params_dict(source_file):

  import string

  pdict = {}
  read_file = open(source_file, "rb")

  # key, tzero, mhzps, flipped

  for line in read_file:

    parameters = string.split(line,",")

    pdict[parameters[0]] = {}

    pdict[parameters[0]]["tzero"] = float(parameters[1])

    pdict[parameters[0]]["mhzps"] = float(parameters[2])

    pdict[parameters[0]]["flipped"] = bool(parameters[3])

  return pdict

def PowerFigures85():
  '''generates images in rootdirectory/figures of power data'''

  dump_dir = root_dir + OSDirAppend('figures')

  params_85 = params_dict(root_dir + OSDirAppend('85params.csv'))

  for filename in os.listdir(source_dir_85):

    print "Plotting " + filename + "\n"

    file_path = source_dir_85 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    mhz_scale = (file_data['t1'] - params_85[filename]['tzero']) * \
      params_85[filename]['mhzps']

    master_aom = file_data['v1']
    pdh_absorb = file_data['v2']
    pdh_error = file_data['v3']

    # flip if flipped

    if(params_85[filename]['flipped']):
      mhz_scale = mhz_scale[::-1] * -1.0
      master_aom = master_aom[::-1]
      pdh_absorb = pdh_absorb[::-1]
      pdh_error = pdh_error[::-1]

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1.plot(mhz_scale, master_aom, linewidth=2, color='#324259')

    ax1b = ax1.twinx()
    ax1b.plot(mhz_scale, pdh_absorb, linewidth=2, color='#4646EA')

    ax1.get_xaxis().set_visible(False)

    ax2 = fig.add_subplot(212)
    ax2.plot(mhz_scale, pdh_error, linewidth=2, color='#ff2b63')

    save_file = dump_dir + OSDirAppend(filename.replace(".csv", ".png"))

    fig.savefig(save_file, bbox_inches='tight')

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break

    pt.close(fig)

  return

def PowerFigures87():

  dump_dir = root_dir + OSDirAppend('figures')

  params_87 = params_dict(root_dir + OSDirAppend('87params.csv'))

  # 87 DATA IS NOT VERY EASY TO EXTRACT PEAKS FROM, WE JUST USE THE DATA
  # FROM THE FIRST EXTRACTION FOR ALL FILES AS THEY ARE SIMILAR ENOUGH

  fzero87 = os.listdir(source_dir_87)[0]

  for filename in os.listdir(source_dir_87):

    print "Plotting " + filename + "\n"

    file_path = source_dir_87 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    mhz_scale = (file_data['t1'] - params_87[fzero87]['tzero']) * \
      params_87[fzero87]['mhzps']

    master_aom = file_data['v1']
    pdh_absorb = file_data['v2']
    pdh_error = file_data['v3']

    # flip if flipped

    if(params_87[filename]['flipped']):
      mhz_scale = mhz_scale[::-1] * -1.0
      master_aom = master_aom[::-1]
      pdh_absorb = pdh_absorb[::-1]
      pdh_error = pdh_error[::-1]

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1.plot(mhz_scale, master_aom, linewidth=2, color='#324259')

    ax1b = ax1.twinx()
    ax1b.plot(mhz_scale, pdh_absorb, linewidth=2, color='#4646EA')

    ax1.get_xaxis().set_visible(False)

    ax2 = fig.add_subplot(212)
    ax2.plot(mhz_scale, pdh_error, linewidth=2, color='#ff2b63')

    save_file = dump_dir + OSDirAppend(filename.replace(".csv", ".png"))

    fig.savefig(save_file, bbox_inches='tight')

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break

    pt.close(fig)

  return

#
# In Shell run:
#
# >> FindResXXPower()
# >> PowerFiguresXX()
#