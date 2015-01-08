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
from math import ceil

from savitzkygolay import SavitzkyGolay


# this is the data root, DROPBOX/ENPH\ 479/rubidium
# specify power/temp/old data directory later

data_dir = ""

if 'nt' in os.name:
  data_dir = "E:\\Documents\\Dropbox\\ENPH 479\\rubidium"
else:
  data_dir = os.environ['DATA479']


root_dir = os.getcwd()

# peak distances calculated using Dan Steck Rb 87, 85 data

co13_co23_87 = 156.947070 / 2.0
co34_co23_85 = 63.40161 / 2.0

tzero = 0.0
mhzps = 1.0
tdelta = 1.0

def OSDirAppend(target_dir):

  if 'nt' in os.name:
    return "\\" + target_dir
  else:
    return "/" + target_dir


source_dir_85 = data_dir + OSDirAppend('power') + OSDirAppend('85')
source_dir_87 = data_dir + OSDirAppend('power') + OSDirAppend('87')
source_dir_temps = data_dir + OSDirAppend('temperature_2')

write_params_85 = root_dir + OSDirAppend('85params.csv')
write_params_87 = root_dir + OSDirAppend('87params.csv')
write_params_temps = root_dir + OSDirAppend('tempparams.csv')
write_params_lock = root_dir + OSDirAppend('lockparams.csv')

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

    pt.plot(t_fit, v_fit, t_fit, v_sgf, linewidth=3)
    pt.vlines(tpeaks, min(v_fit), max(v_fit))
    pt.show()

  wfile.close()

  return


def FindResTemps():

  wfile = open(write_params_temps, 'wb')

  for filename in os.listdir(source_dir_temps):

    print "Processing " + filename

    file_path = source_dir_temps + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    t_fit = file_data['t1']
    v_fit = file_data['v1']

    v_sgf = np.copy(v_fit)

    if "87" in filename:

      # Piecewise fit using
      #     11 points per window (must be odd)
      #     9th order polynomial (must be > 1, < window size)

      num_passes = 30

      for i in xrange(0, num_passes):
        v_sgf = SavitzkyGolay(v_sgf, 51, 3)

      vn = len(v_fit)
      window = int(vn*0.202)

    #   v_sgf2 = v_sgf[window:vn-window]

    #   # NOTE:
    #   #   noticed that sometimes there is a peak near the center, but
    #   #   it is inconsistent in appearence so only the last two peaks are used

    #   peaks = naive_peak_find(v_sgf2, window)[-2:]#

    #   # boolean, true if peaks are in front of center, as they are lower
    #   # actual frequency than the center of the spectrum
    #   flipped = (int(sum(peaks)/float(len(peaks))) > int(vn/3.))

    #   tpeaks = [t_fit[a] for a in peaks]

    #   if flipped:
    #     tzero = tpeaks[-2]
    #     tdelta = tpeaks[-1] - tpeaks[-2]
    #   else:
    #     tzero = tpeaks[1]
    #     tdelta = tpeaks[1] - tpeaks[0]

    #   mhzps = co13_co23_87 / abs(tdelta)

    else:

      num_passes = 50

      for i in xrange(0, num_passes):
        v_sgf = SavitzkyGolay(v_sgf, 51, 3)

    #   vn = len(v_fit)
    #   window = int(vn/4.)

    #   v_sgf2 = v_sgf[window:vn-window]

    #   # NOTE:
    #   #   noticed that sometimes there is a peak near the center, but
    #   #   it is inconsistent in appearence so only the last two peaks are used

    #   peaks = naive_peak_find(v_sgf2, window)[-2:]#

    #   # boolean, true if peaks are in front of center, as they are lower
    #   # actual frequency than the center of the spectrum
    #   flipped = (int(sum(peaks)/float(len(peaks))) > int(vn/3.))

    #   tpeaks = [t_fit[a] for a in peaks]

    #   if flipped:
    #     tzero = tpeaks[-1]
    #     tdelta = tpeaks[-1] - tpeaks[-2]
    #   else:
    #     tzero = tpeaks[0]
    #     tdelta = tpeaks[1] - tpeaks[0]

    #   mhzps = co34_co23_85 / abs(tdelta)

    # wfile.write(",".join([filename, str(tzero), str(mhzps), str(flipped)])+"\n")

    # # UNCOMMENT FOLLOWING TO VISUALLY INSPECT

    #
    # Automatic code is not working well here, do manually unless time to fix
    #

    pt.plot(t_fit, v_fit, t_fit, v_sgf, linewidth=3)
    #pt.vlines(tpeaks, min(v_fit), max(v_fit))
    pt.show()

    flipped = raw_input("Flipped?: ")

    inner_peak = float(raw_input("Inner Peak (s): "))

    outer_peak = float(raw_input("Outer Peak (s): "))

    tzero = inner_peak
    mhzps = 1/(abs(inner_peak - outer_peak))

    if "87" in filename:
      mhzps *= co13_co23_87
    else:
      mhzps *= co34_co23_85

    wfile.write(",".join([filename, str(tzero), str(mhzps), str(flipped)])+"\n")

  wfile.close()

  return

#
#
# The next set of functions deals with plotting figures
#
#

def round_up(num, divisor):
  rem = (num%divisor)
  if rem > 0:
    return num + 5 - num%divisor
  else:
    return num

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

def slope_dict(source_file):

  import string

  pdict = {}
  read_file = open(source_file, "rb")

  # key, tzero, mhzps, flipped

  for line in read_file:

    parameters = string.split(line,",")

    pdict[parameters[0]] = {}

    pdict[parameters[0]]["mvps"] = float(parameters[1])

    pdict[parameters[0]]["dmvps"] = float(parameters[2])

    pdict[parameters[0]]["offset"] = float(parameters[3])

    pdict[parameters[0]]["nmvpp"] = float(parameters[4])

    pdict[parameters[0]]["in1"] = int(parameters[5])

    pdict[parameters[0]]["in2"] = float(parameters[6])

  return pdict

def PowerFigures85():
  '''generates images in rootdirectory/figures of power data'''

  dump_dir = root_dir + OSDirAppend('figures')

  params_85 = params_dict(root_dir + OSDirAppend('85params.csv'))
  slope_85 = slope_dict(root_dir + OSDirAppend('lockparams85.csv'))

  for filename in os.listdir(source_dir_85):

    print "Plotting " + filename + "\n"

    file_path = source_dir_85 + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    mhz_scale = (file_data['t1'] - params_85[filename]['tzero']) * \
      params_85[filename]['mhzps']

    master_aom = file_data['v1']*1000
    pdh_absorb = file_data['v2']*1000
    pdh_error = file_data['v3']*1000

    # lock point slope stuff

    lock_slope = slope_85[filename]['offset'] + \
      file_data['t1'] * slope_85[filename]['mvps']

    mvpmhz = abs(slope_85[filename]['mvps']/params_85[filename]['mhzps'])
    dmvpmhz = slope_85[filename]['dmvps']/params_85[filename]['mhzps']

    nmvpp = slope_85[filename]['nmvpp']

    # flip if flipped

    if(params_85[filename]['flipped']):
      mhz_scale = mhz_scale[::-1] * -1.0
      master_aom = master_aom[::-1]
      pdh_absorb = pdh_absorb[::-1]
      pdh_error = pdh_error[::-1]
      lock_slope = lock_slope[::-1]

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1b = ax1.twinx()
    ax2 = fig.add_subplot(212)

    xmin = -200
    xmax = 600

    ax1.set_xlim([xmin, xmax])
    ax2.set_xlim([xmin, xmax])
    ax1b.set_xlim([xmin, xmax])

    ax1.plot(mhz_scale, master_aom, linewidth=2, color='#324259')
    ax1.set_ylabel('Master (mV, Black)')
    ax1.set_ylim([-360, -180])

    ax1b.plot(mhz_scale, pdh_absorb, linewidth=2, color='#90b2fb')
    ax1b.set_ylabel('Mod Transfer (mV, Blue)')
    ymin, ymax = ax1b.get_ylim()
    diff = round_up(ceil((ymax-ymin)/9.0), 5)
    if diff < 5.0:
      diff = 5.0
    print diff
    ax1b.set_ylim([ymax-9*diff, ymax])
    ax1b.set_yticks([ymax-(9-a)*diff for a in reversed(xrange(0, 10))])
    ax1b.set_yticklabels([ymax-(9-a)*diff for a in reversed(xrange(0, 10))])

    ax2.plot(mhz_scale, pdh_error, linewidth=2, color='#ff2b63')
    ax2.set_ylabel('PDH Error Signal (mV)')
    ax2.set_xlabel('Frequency from F\'=3,4 Crossover (MHz)')

    ax2lims = ax2.get_ylim()

    ax2.plot(mhz_scale, lock_slope, linewidth=2, color='black')

    ax2.set_ylim(ax2lims)

    ax2.annotate("Slope (mV/MHz): %.4f +/- %.4f" % (mvpmhz, dmvpmhz),\
    xy=(0.55, 0.9), xycoords='axes fraction', fontsize=10)
    ax2.annotate("Noise @ Lock (mVpp): %.4f" % (nmvpp),\
    xy=(0.55, 0.8), xycoords='axes fraction', fontsize=10)

    ax1.grid(True, which='both')
    ax1.set_xticklabels([])
    ax2.grid(True, which='both')


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

    master_aom = file_data['v1']*1000
    pdh_absorb = file_data['v2']*1000
    pdh_error = file_data['v3']*1000

    # flip if flipped

    if(params_87[filename]['flipped']):
      mhz_scale = mhz_scale[::-1] * -1.0
      master_aom = master_aom[::-1]
      pdh_absorb = pdh_absorb[::-1]
      pdh_error = pdh_error[::-1]

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1b = ax1.twinx()
    ax2 = fig.add_subplot(212)

    xmin = -200
    xmax = 500

    ax1.set_xlim([xmin, xmax])
    ax2.set_xlim([xmin, xmax])
    ax1b.set_xlim([xmin, xmax])

    ax1.plot(mhz_scale, master_aom, linewidth=2, color='#324259')
    ax1.set_ylabel('Master (mV, Black)')
    ax1.set_ylim([-360, -180])

    ax1b.plot(mhz_scale, pdh_absorb, linewidth=2, color='#90b2fb')
    ax1b.set_ylabel('Mod Transfer (mV, Blue)')
    ymin, ymax = ax1b.get_ylim()
    diff = round_up(ceil((ymax-ymin)/9.0), 5)
    if diff < 5.0:
      diff = 5.0
    print diff
    ax1b.set_ylim([ymax-9*diff, ymax])
    ax1b.set_yticks([ymax-(9-a)*diff for a in reversed(xrange(0, 10))])
    ax1b.set_yticklabels([ymax-(9-a)*diff for a in reversed(xrange(0, 10))])

    ax2.plot(mhz_scale, pdh_error, linewidth=2, color='#ff2b63')
    ax2.set_ylabel('PDH Error Signal (mV)')
    ax2.set_xlabel('Frequency from F\'=2,3 Crossover (MHz)')

    ax1.grid(True, which='both')
    ax1.set_xticklabels([])
    ax2.grid(True, which='both')

    save_file = dump_dir + OSDirAppend(filename.replace(".csv", ".png"))

    fig.savefig(save_file, bbox_inches='tight')

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break

    pt.close(fig)

  return

def TempFigures():

  # just plot for now until we figure out how to calculate the paramaters more
  # reliably

  params_temp = params_dict(root_dir + OSDirAppend('tempparams.csv'))

  dump_dir = root_dir + OSDirAppend('figures')

  for filename in os.listdir(source_dir_temps):

    print "Plotting " + filename + "\n"

    file_path = source_dir_temps + OSDirAppend(filename)
    file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
      names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

    mhz_scale = (file_data['t1']- params_temp[filename]['tzero']) * \
      params_temp[filename]['mhzps']

    master_aom = file_data['v1']*1000
    pdh_absorb = file_data['v2']*1000
    pdh_error = file_data['v3']*1000

    # flip if flipped

    # BUG:
    # figures don't flip even if "flipped" is True in CSV

    if(params_temp[filename]['flipped']):
      mhz_scale = mhz_scale[::-1] * -1.0
      master_aom = master_aom[::-1]
      pdh_absorb = pdh_absorb[::-1]
      pdh_error = pdh_error[::-1]

    fig = pt.figure()

    ax1 = fig.add_subplot(211)
    ax1.plot(mhz_scale, master_aom, linewidth=2, color='#324259')
    ax1.set_ylabel('Master (mV, Black)')

    ax1b = ax1.twinx()
    ax1b.plot(mhz_scale, pdh_absorb, linewidth=2, color='#90b2fb')
    ax1b.set_ylabel('Mod Transfer (mV, Blue)')

    ax1.get_xaxis().set_visible(False)

    ax2 = fig.add_subplot(212)
    ax2.plot(mhz_scale, pdh_error, linewidth=2, color='#ff2b63')
    ax2.set_ylabel('PDH Error Signal (mV)')
    ax2.set_xlabel('Frequency (MHz)')

    save_file = dump_dir + OSDirAppend(filename.replace(".csv", ".png"))

    fig.savefig(save_file, bbox_inches='tight')

    # if(raw_input("Enter c to close PROGRAM: ") == "c"):
    #   break

    pt.close(fig)

  return


def LockPointSlopeNoise():
  'fairly interactive. Will store slope data in mV/s, rely on peak finder \
  to convert to mV/MHz. Will store noise value in mV'

  def LineFit(axislimits, time, error, givenguess=[1.0, 1.0]):

    def func(x, a, b):
      return a*x + b

    import scipy.optimize as opt

    index1 = 0
    index2 = 0

    for i in xrange(0, len(time)):

      if (not index1):
        if (time[i] > axislimits[0]):
          index1 = i

      if (not index2):
        if (time[i] > axislimits[1]):
          index2 = i
          break

    print "Index Bounds: ", index1, index2

    fittime = time[index1:index2]

    fiterror = error[index1:index2]

    guess = givenguess

    #fits = np.polyfit(fittime, fiterror, 1)
    fits, fcov =  opt.curve_fit(func, fittime, fiterror, guess)

    print "Fitting Result:"
    print fits[0]

    mvps = fits[0]
    dmvps = np.sqrt(fcov[0,0])
    offset = fits[1]

    flatline = fiterror - fittime*mvps # removes slope, flat line w/ noise
    mvppn = (max(flatline) - min(flatline)) # vpp, not vrms

    fitline = offset + mvps * fittime

    return fittime, fitline, mvps, dmvps, offset, mvppn, index1, index2, fiterror

  wfile = open(write_params_lock, 'wb')

  pt.ion()
  fig = pt.figure()

  dump_dir = root_dir

  for directory in [source_dir_87]:

    for filename in os.listdir(directory):

      fig.clf()

      file_path = directory + OSDirAppend(filename)
      file_data = np.genfromtxt(file_path, dtype=float, delimiter=',', \
        names=['t1', 'v1', 't2', 'v2', 't3', 'v3'])

      # in seconds

      sec_scale = file_data['t1']

      # in mV

      master_aom = file_data['v1']*1000
      pdh_absorb = file_data['v2']*1000
      pdh_error = file_data['v3']*1000

      ax1 = fig.add_subplot(111)
      ax1.plot(sec_scale, pdh_absorb, linewidth=2, color='#324259')
      ax1.set_ylabel('Mod Transfer (mV, Blue)')

      ax1b = ax1.twinx()
      ax1b.plot(sec_scale, pdh_error, linewidth=2, color='#ff2b63')
      ax1b.set_ylabel('PDH Error Signal (mV)')

      ax1.set_xlabel('Time(s)')

      fig.show()

      raw_input("Reposition Error Signal Around Lock Point,"+ \
        "then press any key to continue...")

      ax1x = ax1.get_xlim()

      print "Time bounds: ", ax1x

      fit_time, fit_line, mvps, dmvps, offset, vppn, in1, in2, fiterror = \
        LineFit(ax1x, sec_scale, pdh_error)

      ax1b.plot(fit_time, fit_line, linewidth=2, color='black')
      ax1b.plot(fit_time, fiterror, linewidth=2, color='green')

      print "MVPS, MVPPN: ", mvps, dmvps, vppn

      prompt = raw_input("good? (y/n/c): ")

      while(prompt != "y" and prompt != "c"):
        raw_input("Reposition Error Signal Around Lock Point,"+ \
          "then press any key to continue...")

        suggestion = raw_input("Specify guessses (slope/offset)?: ")

        guesses = [1.0, 1.0]

        if suggestion != "":
          guesses = np.array(suggestion.split(",")).astype(float)

        ax1x = ax1.get_xlim()

        print "Time bounds: ", ax1x

        fit_time, fit_line, mvps, dmvps, offset, vppn, in1, in2, fiterror = \
          LineFit(ax1x, sec_scale, pdh_error, guesses)

        ax1b.plot(fit_time, fit_line, linewidth=2, color='black')
        ax1b.plot(fit_time, fiterror, linewidth=2, color='green')

        print "MVPS, MVPPN: ", mvps, vppn

        prompt = raw_input("good? (y/n): ")

      if prompt == "c":
        return

      wfile.write(",".join([filename, str(mvps), str(dmvps), str(offset), str(vppn), \
        str(in1), str(in2)])+"\n")

    pt.close(fig)

    wfile.close()

    return


#
# In Shell run:
#
# >> FindResXXPower()
# >> PowerFiguresXX()
#