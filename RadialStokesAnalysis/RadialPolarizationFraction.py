__author__ = 'Mehrnoosh Tahani'
# Created at 2020-06-26

import os, sys, subprocess, shlex, glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from Classes.RadialPatternFinding import RadialParam
from scipy.interpolate import UnivariateSpline
import statistics as st
import matplotlib.ticker as ticker

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
def plotandSave(PfFunc, PIFunc, radius, PfFuncStdev, PiFuncStdev, nCount, bubbleNumber, cx, cy, r, tx, directory):
    saveFilePath = directory + 'Bubble' + str(bubbleNumber) + 'RadialPF.pdf'
    fig, ax1 = plt.subplots(figsize=(8,6))
    ax2 = ax1.twinx()
    ax1.set_ylabel(r'Polarization Fraction (%)', fontsize=18, labelpad=-3)
    ax1.set_xlabel('radius (arc-second)', fontsize=18, labelpad=7)
    ax1.tick_params(labelsize=12)
    plt.yticks(fontsize=12)
    plt.title('Bubble ' + str(bubbleNumber)  , fontsize=18)
    ax2.set_ylabel('Number of pixels', fontsize=14, color='black', labelpad=-3)
    forLeg1 = ax1.errorbar(radius,  PfFunc, yerr=PfFuncStdev, color='blue', marker='o', ls='none',label=r'<Q$_r$>', markersize = 6)
    ax1.axvline(x=r*3.999999, label='bubble radius', ls ='--', c = 'black')
    # For Simpson bubbles we have thickness as well:
    forLeg3 = ax2.scatter(radius, nCount,color='black', marker='x', label = 'Pixels', s = 14)
    plt.legend((forLeg1, forLeg3), ('<PF>', 'Total Pixels'), loc=[0.02, 0.81])
    spl = UnivariateSpline(radius, PfFunc)
    xs = np.linspace(xmin2, xmax, xmax-xmin2)
    ax1.plot(xs, spl(xs), lw=1.2, linestyle='--', color='blue', alpha = 0.5)
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

