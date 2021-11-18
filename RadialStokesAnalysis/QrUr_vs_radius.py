__author__ = 'Mehrnoosh Tahani'
# Created at 2020-06-26

import os, sys, subprocess, shlex, glob, json
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
import astropy as ast
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import cm
import statistics as st
import matplotlib.ticker as ticker
from scipy.interpolate import UnivariateSpline
from Classes.RadialPatternFinding import RadialParam
from scipy.interpolate import splev, splrep

# ~~~~~~~~~ File Paths ~~~~~~~~~~~~
parentDir = os.path.abspath(os.path.dirname(os.getcwd()))
dataDir = os.path.join(parentDir, 'Data/')
outDir = os.path.join(parentDir, 'Output/')
SaveFileDirectory = os.path.join(outDir, '3D/')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Cutoff values ~~~~~~~~~
StokesI_cutoff = 10.
PolarizationI_cutoff = 3.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~Plot Q/Ur vs radius~~~~~~~~~
def plotandSave(qrFunc, urFunc, radius, qrFuncStdev, urFuncStdev, nCount, bubbleNumber, cx, cy, r, tx, directory):
    # ~~~~~~~~~~~Xmin and Xmax~~~~~~ (non-efficient coding for now)
    # ~~~~~~~~~~end of xmin and xmax ~~~~~~~~~
    saveFilePath = directory + 'Bubble' + str(bubbleNumber) + 'QURadialNoI.pdf'
    fig, ax1 = plt.subplots(figsize=(8,6))
    ax2 = ax1.twinx()
    ax1.set_ylabel(r'<Q$_r$> or <U$_r$> (mJy/beam)', fontsize=18, labelpad=-6)
    ax1.set_xlabel('radius (arc-second)', fontsize=18, labelpad=7)
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(tickLocator))
    # plt.title('SNR(I)>'+str(StokesI_cutoff) +', SNR(PI)> '+str(PolarizationI_cutoff)+ 'Bubble ' + str(bubbleNumber) , fontsize=14)
    # For NO PI selection Criteria
    # plt.title('SNR(I)>'+str(StokesI_cutoff) +', No SNR(PI), Bubble ' + str(bubbleNumber) +'\n', linespacing =1 , fontsize=14)
    plt.title('Bubble ' + str(bubbleNumber)  , fontsize=18)
    ax2.set_ylabel('Number of pixels', fontsize=14, color='black', labelpad=-3)
    # ax1.plot(radius, qrFunc, color='green', marker='o', label= r'<Q$_r$>', markersize=4)
    forLeg1 = ax1.errorbar(radius, qrFunc, yerr=qrFuncStdev, color='blue', marker='o', ls='none',label=r'<Q$_r$>', markersize = 6)
    # ax1.plot(radius, urFunc, color='red', marker='s', label= r'<U$_r$>', markersize=4)
    forLeg2 = ax1.errorbar(radius, urFunc, yerr=urFuncStdev, color='red', marker='o', ls='none',label=r'<Q$_r$>', markersize = 6, alpha = 0.5)
    ax1.set_ylim(-25,25)
    ax1.set_xlim(xmin,xmax)
    ax1.tick_params(labelsize=12)
    plt.yticks(fontsize=12)
    ax1.axvline(x=r*3.999999, label='bubble radius', ls ='--', c = 'black')
    # For Simpson bubbles we have thickness as well:
    forLeg3 = ax2.scatter(radius, nCount,color='black', marker='x', label = 'Pixels', s = 14)
    # ax1.legend(bbox_to_anchor=(-0.34, .4, 0.6, 0.6))
    # plt.savefig(saveFilePath, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.legend((forLeg1, forLeg2, forLeg3), (r'<Q$_r$>', r'<U$_r$>', 'Total Pixels'), loc=[0.02, 0.81])
    spl = UnivariateSpline(radius, qrFunc)
    spl2 = UnivariateSpline(radius, urFunc)
    xs = np.linspace(xmin2, xmax, xmax-xmin2)
    yu = np.zeros(len(radius))
    ax1.plot(xs, spl(xs), lw=1.2, linestyle='--', color='blue', alpha = 0.5)
    # ax1.plot(xs, yu, lw=1.2, linestyle = '-.', color='red', alpha = 0.4)
    ax1.plot(radius, yu, lw=1.2, linestyle='-.', color='red', alpha=0.4)
    plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none')
    # plt.show()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~ Header and WCS ~~~~~~~
ifits_name = os.path.join(dataDir, 'iextN633MO808YYNN.fit')
ihdu = fits.open(ifits_name)
ifhd = ihdu[0].header
ifbd = ihdu[1].data
wcs = WCS(ifhd).celestial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~Find the bubbles info~~~~~~~~~.
bubbles = BubbleShortenedPicked(wcs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
bubbleIndex = 0
bubbleThickness = 0
# ~~Make Results for each Bubble~~~~.
for bubbleIndex in range(len(bubbles.RegionRA)):
    bubbleNumber = bubbleIndex + 1
    if bubbleNumber !=11 and bubbleNumber !=12 and bubbleNumber !=13 and bubbleNumber !=2:
        maskedRegion = RadialParam(dataDir, bubbles.RegionX[bubbleIndex], bubbles.RegionY[bubbleIndex], bubbles.RegionRadius[bubbleIndex], bubbles.RegionThickinArcSec[bubbleIndex], StokesI_cutoff, PolarizationI_cutoff)
        plot = plotandSave(maskedRegion.pixelinEachAngleBin, maskedRegion.angleBin, maskedRegion.radius, bubbleNumber, SaveFileDirectory)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

