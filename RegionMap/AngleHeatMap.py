__author__ = 'Mehrnoosh Tahani'

import os, sys, subprocess, shlex, glob, json
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec # GRIDSPEC !

# ~~~~~~~~~ File Paths ~~~~~~~~~
currentDir = os.path.abspath(os.getcwd())
dataDir = os.path.join(currentDir, 'Data/')
outDir = os.path.join(currentDir, 'Output/')
saveFilePath = os.path.join(outDir, 'AngleHeatMap.pdf')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Cutoff values ~~~~~~~~~
I_cutoff = 10
PI_cutoff = 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~File Names~~~~~~~~~~~~~~~~~~~~~~~
iext_name = 'iextN633MO808YYNN.fit'
iextInsideAST_name = 'iextN633MO808YYNNinsideAST.fit'
Qext_name = 'qextN633MO808YYNN.fit'
Uext_name = 'uextN633MO808YYNN.fit'
Piext_name = 'piextN633MO808YYNN.fit'
PFext_name = 'pfextN633MO808YYNN.fit'
# ~~~~~~~~~~~~~~~End of File Names~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~ Header and WCS ~~~~~~~
ifits_name = os.path.join(dataDir, iext_name)
ihdu = fits.open(ifits_name)
ifhd = ihdu[0].header
ifbd = ihdu[1].data
wcs = WCS(ifhd).celestial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Stokes HDU~~~~~~~~~~~~~~~~~~~~~~~
ifitsInsideAST = os.path.join(dataDir, iextInsideAST_name)
ihduInsideAST = fits.open(ifitsInsideAST)
isdtInsideAST = ihduInsideAST['PRIMARY', 1].data[0]
iedtInsideAST = ihduInsideAST['VARIANCE', 1].data[0]
# Stokes I
isdt = ihdu['PRIMARY', 1].data[0]
iedt = ihdu['VARIANCE', 1].data[0]
# Stokes Q
Qfits_name = os.path.join(dataDir, Qext_name)
Qhdu = fits.open(Qfits_name)
Qsdt = Qhdu['PRIMARY', 1].data[0]
Qedt = Qhdu['VARIANCE', 1].data[0]
# Stokes U
Ufits_nameCl = os.path.join(dataDir, Uext_name)
Uhdu = fits.open(Ufits_nameCl)
Usdt = Uhdu['PRIMARY', 1].data[0]
Uedt = Uhdu['VARIANCE', 1].data[0]
# Polarized I
Pifits_name = os.path.join(dataDir, Piext_name)
Pihdu = fits.open(Pifits_name)
Pisdt = Pihdu['PRIMARY', 1].data[0]
Piedt = Pihdu['VARIANCE', 1].data[0]
# Polarized Fraction
PFfits_name = os.path.join(dataDir, PFext_name)
PFhdu = fits.open(PFfits_name)
PFsdt = PFhdu['PRIMARY', 1].data[0]
PFedt = PFhdu['VARIANCE', 1].data[0]
# B Angle
copiedU = Usdt.copy()
copiedQ = Qsdt.copy()
Anglesdt = np.degrees(0.5* np.arctan2(copiedU, copiedQ))
BAngleRotated = Anglesdt+90
# Apply the selection criteria
BAngleRotated[isdt/iedt<=I_cutoff]=np.nan
BAngleRotated[Pisdt/Piedt<=PI_cutoff]=np.nan
# ~~~~~~~~~~~~~~~~~~~~~~~~ All Stokes HDU~~~~~~~~~~~~~~~~~~~~~~~.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating Fig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = plt.figure(facecolor='w', figsize=(12, 10), dpi=100)
gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
gs.update(left=0.16, right=0.76, bottom=0.08, top=0.93, wspace=0.0, hspace=0.0)
plt.rcParams['axes.titlesize'] = 25
ax = plt.subplot(gs[0,0], projection=wcs, facecolor='w')

lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('Right Ascention (J2000)', minpad=0.8, fontsize=20)
lat.set_axislabel('Declination (J2000)', minpad=-0.3, fontsize=20)
lon.set_ticklabel(size=18)
lat.set_ticklabel(size=18)
# lon.set_major_formatter('hh:mm:ss.s')
# lat.set_major_formatter('dd:mm:ss.s')
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_separator(('h', "m", 's'))
ax.set_xlim(50,345)
ax.set_ylim(10, 375)
# Tick/label spacing and properties
lon.set_ticks(spacing=60 * u.arcsec, color='black', exclude_overlapping=True)
lat.set_ticks(spacing=2 * u.arcmin, color='black', exclude_overlapping=True)
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

im = ax.imshow(BAngleRotated, cmap=plt.cm.seismic, aspect='equal')#gist_yarg
# Color bar
cbax = plt.subplot(gs[0,1])
cb = plt.colorbar(cax=cbax, mappable=im, orientation='vertical', ticklocation='right')
cb.ax.tick_params(labelsize=18) #  Fontsize of colorbar values
cb.set_label(r'Magnetic field angle [$^{\circ}$]', fontsize=18, labelpad=18)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating Fig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

# plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.show()
