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
dataDir = os.path.join(currentDir, 'Data/')
outDir = os.path.join(currentDir, 'Output/')
saveFilePath = os.path.join(outDir, 'AngleHeatMap.pdf')
iextInsideAST_name = 'iextN633MO808YYNNinsideAST.fit'
iext_name = 'iextN633MO808YYNN.fit'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Cutoff values ~~~~~~~~~
I_cutoff = 10
PI_cutoff = 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Stokes HDU~~~~~~~~~~~~~~~~~~~~~~~
ifitsInsideAST = os.path.join(dataDir, iextInsideAST_name)
ihduInsideAST = fits.open(ifitsInsideAST)
isdtInsideAST = ihduInsideAST['PRIMARY', 1].data[0]
iedtInsideAST = ihduInsideAST['VARIANCE', 1].data[0]
# Stokes I
ifits_name = os.path.join(dataDir, iext_name)
ihdu = fits.open(ifits_name)
isdt = ihdu['PRIMARY', 1].data[0]
iedt = ihdu['VARIANCE', 1].data[0]

#AST mask
isdt[np.isnan(isdtInsideAST) ] = 0

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
catfit_name = 'polvecatN633MO808YYNN.FIT'
catfits = os.path.join(dataDir, catfit_name)
hdu_list = fits.open(catfits, memmap=True)
#print Table(hdu_list[1].data)
polvec_data = Table(hdu_list[1].data) # to astropy Table
# MT added
mdf = polvec_data.to_pandas()
mcdf = mdf.dropna()

polvec_dataSelected2 = mcdf[ (mcdf.I / mcdf.DI > I_cutoff) &
             (mcdf.PI / mcdf.DPI > PI_cutoff) &
             (mcdf.P >0.) &
            (mcdf.P <100.)]
polvec_dataSelected2[['P', 'DP']].describe()

polvec_dataSelected = Table.from_pandas(polvec_dataSelected2)

RA0 = polvec_dataSelected['RA'] *180./np.pi + 360. #
DEC0 = polvec_dataSelected['DEC']*180/np.pi
'''
Becase polvec_data['P'] is given in percent, e.g., 20.0,
 first, multiply 1/100 to obtain fraction, yielding 1 degree = 1 %.
 second, multiply 1/3600 to convert 1 arcsecond = 1%
 third, multiply fhd['CDELT2']*3600 = 4 arcsec, so that 1 pixel = 1%
'''
pf2pixsize = 1/100. * 1/3600. * ifhd['CDELT2']*3600.  #  P=1.0% corresponds to 4 arcsec # the 1/100 reduces the size of the vectors, was originally 1/100
# # STARTING POINTS
rotang = 90. #for the magnetic field, instead of dust polarization
RA1  = RA0  + pf2pixsize*polvec_dataSelected['P'] * np.sin((polvec_dataSelected['ANG']+rotang)/180.*np.pi)*180/np.pi
DEC1 = DEC0 + pf2pixsize*polvec_dataSelected['P'] * np.cos((polvec_dataSelected['ANG']+rotang)/180.*np.pi)*180/np.pi
# ENDING POINTS
RA2  = RA0  - pf2pixsize*polvec_dataSelected['P'] * np.sin((polvec_dataSelected['ANG']+rotang)/180.*np.pi)*180/np.pi
DEC2 = DEC0 - pf2pixsize*polvec_dataSelected['P'] * np.cos((polvec_dataSelected['ANG']+rotang)/180.*np.pi)*180/np.pi
# Plotting vectors
step = 1
for x in range(0, RA0.shape[0], step):
    l = mlines.Line2D([RA1[x],RA2[x]], [DEC1[x],DEC2[x]], linewidth=1, color='blue', transform=ax.get_transform(ifhd['RADESYS'].lower()), alpha=0.8, clip_on=True)
    ax.add_line(l)
# ~~~~~~~~~~~~~~~~~~~~~~~~~ Ploting Field Lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.show()
