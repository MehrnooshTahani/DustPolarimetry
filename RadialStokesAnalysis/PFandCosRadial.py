__author__ = 'Mehrnoosh Tahani'
# Created at 2020-06-26

import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import UnivariateSpline
from Classes.RadialPatternFinding import RadialParam
from scipy.interpolate import splev, splrep
from Classes.BubblesCircleShortened import BubbleShortenedPicked

# ~~~~~~~~~ File Paths ~~~~~~~~~~~~
parentDir = os.path.abspath(os.path.dirname(os.getcwd()))
dataDir = os.path.join(parentDir, 'Data/')
outDir = os.path.join(parentDir, 'Output/')
SaveFileDirectory = os.path.join(outDir, 'FinalCos/')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Cutoff values ~~~~~~~~~
StokesI_cutoff = 10.
PolarizationI_cutoff = 3.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~Plot Q/Ur vs radius~~~~~~~~~
# Plot Specs
xmin,xmax, xmin2 = 0, 100, 0
tickLocator = 20 #you can define these as dynamic variables
def plotandSave(angFunc, PFFunc, radius, angFuncStdev, PFFuncStdev, nCount, bubbleNumber, cx, cy, r, tx, directory, name):
    # ~~~~~~~~~~~Xmin and Xmax~~~~~~ (non-efficient coding for now)
    # ~~~~~~~~~~end of xmin and xmax ~~~~~~~~~
    saveFilePath = directory + 'Bubble' + str(bubbleNumber) + 'PFCosRadialAngle.pdf'
    fig, ax1 = plt.subplots(figsize=(8,6))
    plt.gcf().subplots_adjust(right=0.8, left=0.15)
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    ax1.set_ylabel(r'$\cos(\theta_r)$', fontsize=19, labelpad=0)
    ax1.set_xlabel('radius [arc-second]', fontsize=19, labelpad=7)
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(tickLocator))
    ax1.tick_params(axis='both', which='major', labelsize=14)
    plt.title('Bubble ' + str(bubbleNumber)  , fontsize=18)
    ax2.set_ylabel('Number of pixels', fontsize=14, color='black', labelpad=-3)
    ax2.set_ylim(0, 205)
    ax3.set_ylim(0, 18)
    ax3.set_ylabel(r'<Polarization Fraction> (%)', fontsize=15, color='black', labelpad=0.5)
    ax2.spines["right"].set_position(("axes", 1.1))
    forLeg1 = ax1.errorbar(radius, angFunc, yerr=angFuncStdev, color='blue', marker='o', ls='none',label=r'<Q$_r$>', markersize = 6)
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(-0.05,1.05)
    ax1.axhline(y=0.95, label='bubble radius', ls ='--', c = 'black', alpha = 0.6)
    ax1.axvline(x=r*3.999999, label='bubble radius', ls ='--', c = 'black')
    # For Simpson bubbles we have thickness as well:
    forLeg3 = ax2.scatter(radius, nCount,color='black', marker='x', label = 'Pixels', s = 14)
    forLeg4 = ax3.errorbar(radius, PFFunc, yerr=PFFuncStdev, color='red', marker='o', ls='none',label=r'<Q$_r$>', markersize = 6)
    leg= plt.legend((forLeg1, forLeg4, forLeg3), (r'$\cos{(\theta_r)}$', r'<PF>', r'total pixels'), title = r"  $\bf{Bubble~" + name + "}$", loc=[0.02, 0.02], fontsize = 'x-large')
    leg._legend_box.align = "left"
    # leg.get_frame().set_alpha(0.5)
    leg.get_frame().set_facecolor((1, 1, 1, 0.5))
    leg.get_frame().set_edgecolor((0,0,0, 0.8))
    leg.get_title().set_fontsize('14')
    plt.tight_layout()
    # plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.show()
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
        plot = plotandSave(maskedRegion.cosAvgAfterAveg, maskedRegion.Pf_rAvg, maskedRegion.radius, maskedRegion.cosAvgAfterAveg_Stdev, maskedRegion.Pf_rStdev, maskedRegion.N_r, bubbleNumber, bubbles.RegionX[bubbleIndex], bubbles.RegionY[bubbleIndex], bubbles.RegionRadius[bubbleIndex], bubbles.RegionThickinArcSec[bubbleIndex], SaveFileDirectory)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

