__author__ = 'Mehrnoosh Tahani; also using original code by Ray Furuya'

import os, sys, subprocess, shlex, glob, json
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

# ~~~~~~~~~ File Paths ~~~~~~~~~
currentDir = os.path.abspath(os.getcwd())
dataDir = os.path.join(currentDir, '..', 'Data/')
outDir = os.path.join(currentDir, '..', 'Output/')
saveFilePath = os.path.join(outDir, 'AngleHeatMap.pdf')
iext_name = 'iextN633MO808YYNN.fit'
Qext_name = 'qextN633MO808YYNN.fit'
Uext_name = 'uextN633MO808YYNN.fit'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Cutoff values ~~~~~~~~~
I_cutoff = 10
PI_cutoff = 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Stokes HDU~~~~~~~~~~~~~~~~~~~~~~~
# Stokes I
ifits_name = os.path.join(dataDir, iext_name)
ihdu = fits.open(ifits_name)
isdt = ihdu['PRIMARY', 1].data[0]
iedt = ihdu['VARIANCE', 1].data[0]
# Stokes Q
Qfits_name = os.path.join(dataDir, Qext_name)
Qhdu = fits.open(Qfits_name)
Qsdt = Qhdu['PRIMARY', 1].data[0]
QMap = Qsdt.copy()
# Stokes U
Ufits_name = os.path.join(dataDir, Uext_name)
Uhdu = fits.open(Ufits_name)
Usdt = Uhdu['PRIMARY', 1].data[0]
UMap = Usdt.copy()

# ~~~~~~~~~~ Header and WCS ~~~~~~~
ifits_name = os.path.join(dataDir, iext_name)
ihdu = fits.open(ifits_name)
ifhd = ihdu[0].header
ifbd = ihdu[1].data
wcs = WCS(ifhd).celestial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating Fig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wx0, wy0 = wcs.wcs_pix2world(100., 100., 1)
fig = plt.figure(facecolor='w', figsize=(12, 10), dpi=100)

gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
gs.update(left=0.16, right=0.85, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
plt.rcParams['axes.titlesize'] = 25
ax = plt.subplot(gs[0,0], projection=wcs, facecolor='w')
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('Right Ascention (J2000)', minpad=0.75, fontsize=20)
lat.set_axislabel('Declination (J2000)', minpad=-0.3, fontsize=20)
lon.set_ticklabel(size=18)
lat.set_ticklabel(size=18)
# lon.set_major_formatter('hh:mm:ss.s')
# lat.set_major_formatter('dd:mm:ss.s')
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_separator(('h', "m", 's'))
# Tick/label spacing and properties
lon.set_ticks(spacing=60 * u.arcsec, color='black', exclude_overlapping=True)
lat.set_ticks(spacing=60 * u.arcsec, color='black', exclude_overlapping=True)
# Minor ticks
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lat.set_minor_frequency(10)
# Tick, tick label, and axis label positions
lon.set_ticks_position('b') # bottom or top, if both "bt"
lon.set_ticklabel_position('b')
lon.set_axislabel_position('b')
lat.set_ticks_position('l') # left or right, if both 'lr'
lat.set_ticklabel_position('l')
lat.set_axislabel_position('l')

im = ax.imshow(isdt, cmap=plt.cm.gist_yarg, aspect='equal', vmin= 0, vmax = 3500)#gist_yarg
# Color bar
cbax = plt.subplot(gs[0,1]) # Place it where it should be defined by GRIDSPEC.
cb = plt.colorbar(cax=cbax, mappable=im, orientation='vertical', ticklocation='right', norm=mpl.colors.Normalize(vmin=0, vmax=5000))
cb.ax.tick_params(labelsize=20) #  Fontsize of colorbar values
cb.set_label(r'Intensity [mJy beam$^{-1}$]', fontsize=18, labelpad=20)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating Fig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~~~~~~~~~ Ploting Field Lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Ploting the B lines
ConvertCoordAngle = []
rotate = 90 # rotating e vectors to get b vectors
step = 10 # show lines for every 10 steps
scale = 1.5 #length of the lines
p_lenx =  QMap.shape[1] #length of the x axis in pixels
p_leny =  QMap.shape[0] #length of the y axis in pixels
data_p = [ [1]*len(QMap[0]) for i in range(len(QMap)) ] #magnitude of the vectors = choosing length of 1, you can change this to make it equal to the polarization fraction for each point
data_a = np.degrees( 0.5 * np.arctan2( UMap , QMap  ) )  # angle of the vectors in degrees
# Loop to add line by line
for y in range(0, p_leny, step):
    for x in range(0, p_lenx, step):
        r = data_p[y][x] * scale
        a = np.radians(data_a[y][x] + rotate)
        x1 = x + r * np.sin(a)
        y1 = y - r * np.cos(a)
        x2 = x - r * np.sin(a)
        y2 = y + r * np.cos(a)
        # Add the lines:
        l = mlines.Line2D([x1, x2], [y1, y2], color='red',linewidth=2, zorder=2, alpha = 0.6)
        ax.add_line(l)
# ~~~~~~~~~~~~~~~~~~~~~~~~~ Ploting Field Lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.show()
