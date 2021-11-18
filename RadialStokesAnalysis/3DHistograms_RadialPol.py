__author__ = 'Mehrnoosh Tahani'

'''
This file plots 3D histograms for radial polarization angles
'''

import os, sys, subprocess, shlex, glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from Classes.RadialPatternFinding import RadialParam
import statistics as st
import matplotlib.ticker as ticker

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Directory Paths ~~~~~~~~~~~~~~~~~
parentDir = os.path.abspath(os.path.dirname(os.getcwd()))
dataDir = os.path.join(parentDir, 'Data/')
outDir = os.path.join(parentDir, 'Output/')
SaveFileDirectory = os.path.join(outDir, '3D/')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~ Selection Criteria ~~~~~~~~~~~
StokesI_cutoff = 10.
PolarizationI_cutoff = 3.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
# ~~~~~~~~~~~~~~~~~~ Reading Fits Files and wcs ~~~~~~~~~~~~~
# Stokes I
iext_name = 'iextN633MO808YYNN.fit'
ifits_name = os.path.join(dataDir, iext_name)
ihdu = fits.open(ifits_name)
ifhd = ihdu[0].header
isdt = ihdu['PRIMARY',1].data[0] # Science data in mJy/beam, this has I,
iVardt = ihdu['VARIANCE',1].data[0] # Error data in mJy/beam
wcs = WCS(ifhd).celestial # coordinate system of the file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

# ~~~~~~Plot Q/Ur vs radius~~~~~~~~~
def plotandSave(pixels, angleBin, radius, bubbleNumber, directory):
    # Path
    saveFilePath = directory + 'Bubble' + str(bubbleNumber) + '3DRarialAngle.pdf'
    # Figure shape
    fig = plt.figure(num = None,  figsize=(15,10), dpi=100 ,facecolor= 'w', edgecolor='k')
    ax = fig.add_subplot(111, projection = '3d')
    for i, r in enumerate(radius):
        if i<plt.cm.tab20c.N:
            colors = plt.cm.tab20c(i)
        elif i-plt.cm.tab20c.N< plt.cm.Pastel1.N:
            colors= plt.cm.Pastel1(i-plt.cm.tab20c.N)
        else:
            colors= plt.cm.Pastel2(i-plt.cm.tab20c.N - plt.cm.Pastel1.N)
        ax.bar(angleBin, pixels[:, i], zs=r, zdir='y', color=colors, width= 15, alpha=0.8, edgecolor='b')
    ax.set_zlabel('\n' + 'Number of pixel', fontsize=18, linespacing=1)
    ax.set_xlabel('\n' + 'Angle (degree)', fontsize=18, linespacing=1)
    ax.set_ylabel('\n' + 'radius (arc-sec)', fontsize=18, linespacing=1)
    ax.set_title('Bubble {}'.format(bubbleNumber), fontsize=19, y=1)
    ax.tick_params(labelsize=12)
    ax.set_xlim(-90, 90)
    plt.xticks(angleBin, angleBin)
    # plt.show()
    plt.savefig(saveFilePath, dpi=90, facecolor=fig.get_facecolor(), edgecolor='none')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

# ~~~~Class to mask regions and find rings for each bubble~~~~

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

