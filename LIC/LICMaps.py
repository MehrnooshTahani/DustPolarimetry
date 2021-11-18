__author__ = 'Mehrnoosh Tahani; based on a code from Kate Pattle and a function from Susan Clark'

'''IMPORTNAT NOTE: LIC first needs to be installed: pip install licpy'''

import numpy as np
from reproject import reproject_interp
from licpy.lic import runlic
import colorsys
import seaborn as sns
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

#~~~~~~~~~~~~~LIC Function ~~~~~~~~~~~~~~~~~~~~
def lic_plot(lic_data, background_data=None, F_m=0.0, F_M=0.0, cmap="YlOrRd"):
    """
    Code to visualize an LIC plot. By Susan Clark.

    lic_data        :: output of LIC code
    background_data :: background color map, e.g. density or vector magnitude
    F_m             :: contrast enhancement parameter - see below
    F_M             :: contrast enhancement parameter - see below
    cmap            :: matplotlib recognized colormap

    Contrast Enhancement from http://www.paraview.org/Wiki/ParaView/Line_Integral_Convolution#Image_LIC_CE_stages
    L_ij = (L_ij - m) / (M - m)
    L = HSL lightness. m = lightness to map to 0. M = lightness to map to 1.
    m = min(L) + F_m * (max(L) - min(L))
    M = max(L) - F_M * (max(L) - min(L))
    F_m and F_M take values between 0 and 1. Increase F_m -> darker colors. Increase F_M -> brighter colors.

    """

    # 1. Compute nhi values
    # 2. Interpolate these values onto cmap to find corresponding RGBA value
    # 3. Convert RGB value to HSV, use these Hues + Saturations
    # 4. Assign Lightness as lic amplitude
    # 5. Display HLS map.

    # Normalize background data
    # if background_data == None:
    #    background_data = np.ones(lic_data.shape)
    hues = background_data / np.nanmax(background_data)

    sats = np.ones(lic_data.shape)
    licmax = np.nanmax(lic_data)
    licmin = np.nanmin(lic_data)

    # Contrast enhancement
    m = licmin + (F_m * (licmax - licmin))
    N = licmax - (F_M * (licmax - licmin))
    vals = (lic_data - m) / (N - m)

    y, x = hues.shape

    # Map background data onto RGB colormap
    clmap = ListedColormap(sns.color_palette(cmap, 256))
    background_data_rgb = clmap(background_data / np.nanmax(background_data))

    # Only need RGB, not RGBA
    background_data_rgb = background_data_rgb[:, :, 0:3]

    # Map to Hue - Saturation - Value
    hsv = mpl.colors.rgb_to_hsv(background_data_rgb)

    # to work in hls instead of hsv
    hs = hsv[:, :, 0].flatten()
    ls = vals.flatten()
    ss = hsv[:, :, 1].flatten()
    r = np.zeros(len(hues.flatten()))
    g = np.zeros(len(hues.flatten()))
    b = np.zeros(len(hues.flatten()))

    maxls = np.nanmax(ls)
    minls = np.nanmin(ls)

    # Translate HLS to RGB
    for i in range(len(hues.flatten())):
        r[i], g[i], b[i] = colorsys.hls_to_rgb(hs[i], ls[i], ss[i])

    r = r.reshape(lic_data.shape)
    g = g.reshape(lic_data.shape)
    b = b.reshape(lic_data.shape)

    rgb = np.zeros((y, x, 3), np.float_)
    rgb[:, :, 0] = r
    rgb[:, :, 1] = g
    rgb[:, :, 2] = b

    return rgb
#~~~~~~~~~~~~~LIC Function ~~~~~~~~~~~~~~~~~~~~.

# ~~~~~~~~~~~~~~ Paths ~~~~~~~~~~
'''If you use the code as it is, you just need to put your fits file in the Data directory.
The next 3 lines read the path, and you can run it without any changes'''
parentDir = os.path.abspath(os.path.dirname(os.getcwd()))
dataDir = os.path.join(parentDir, 'Data/')
outDir = os.path.join(parentDir, 'Output/')
# ~~~~~~~~~~~~~ Paths ~~~~~~~~~~~.
# ~~~~~~~ Selection Criteria ~~~~
I_cutoff = 10
PI_cutoff = 3
# ~~~~~~~ Selection Criteria ~~~~.
# ~~~~~~~~~~ Read Stokes Param ~~~~~
iext_name = 'iextN633MO808YYNN.fit'
Qext_name = 'qextN633MO808YYNN.fit'
Uext_name = 'uextN633MO808YYNN.fit'
# Stokes I
iext_hdulist = fits.open(dataDir + iext_name)
isdt = iext_hdulist['PRIMARY', 1].data[0]
iedt = iext_hdulist['VARIANCE',1].data[0]
IMap = isdt.copy()
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
# WCS
Forwcs = iext_hdulist[0].header
wcs = WCS(Forwcs).celestial
# ~~~~~~~~~~ Read Stokes Param ~~~~~.

# ~~~~~~~~~ Apply Selection Criteria ~~~~~~~~
QMap[(isdt / np.sqrt(iedt)) <= I_cutoff] = np.nan
UMap[(isdt / np.sqrt(iedt)) <= I_cutoff] = np.nan
# ~~~~~~~~~ Apply Selection Criteria ~~~~~~~~.
ang = np.degrees(0.5 * np.arctan2(UMap, QMap))

# ~~~~~~~~~~~~ Creat & Run LIC ~~~~~~~~~~~~~
tex = runlic(np.sin(np.radians(ang)), np.cos(np.radians(ang)), 5)
clrmap = 'brg' # can try other color maps
LICMap = lic_plot(tex, isdt, cmap=clrmap, F_M=1, F_m=1.3)
# ~~~~~~~~~~~~ Creat & Run LIC ~~~~~~~~~~~~~.
# ~~~~~~~~~~~~ Make the Final Plot ~~~~~~~~~
import matplotlib.gridspec as gridspec
fig = plt.figure(facecolor='w', figsize=(8, 6), dpi=100)
gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[1,0.05])
ax = plt.subplot(gs[0,0], projection=wcs, facecolor='w')
im = ax.imshow(LICMap.data, origin='lower',alpha=0.6)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('Right Ascention (J2000)', minpad=0.75, fontsize=14)
lat.set_axislabel('Declination (J2000)', minpad=-0.3, fontsize=14)
lon.set_ticklabel(size=12)
lat.set_ticklabel(size=12)
# lon.set_major_formatter('hh:mm:ss.s')
# lat.set_major_formatter('dd:mm:ss.s')
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_separator(('h', "m", 's'))
# Tick/label spacing and properties
lon.set_ticks(spacing=4 * u.arcmin, color='black', exclude_overlapping=True)
lat.set_ticks(spacing=4 * u.arcmin, color='black', exclude_overlapping=True)
# Minor ticks
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lat.set_minor_frequency(5)
# Tick, tick label, and axis label positions
lon.set_ticks_position('b') # bottom or top, if both "bt"
lon.set_ticklabel_position('b')
lon.set_axislabel_position('b')
lat.set_ticks_position('l') # left or right, if both 'lr'
lat.set_ticklabel_position('l')
lat.set_axislabel_position('l')
# saveFilePath = outDir + 'NGC6334LIC.pdf'
# plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none')
plt.show()
# ~~~~~~~~~~~~ Make the Final Plot ~~~~~~~~~.
