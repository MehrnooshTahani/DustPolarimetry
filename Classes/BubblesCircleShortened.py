__author__ = 'Mehrnoosh Tahani'

import os, sys, subprocess, shlex, glob, json
import numpy as np
# import scipy as sc
import astropy.units as u
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy.coordinates import FK5
from astropy.table import Table

long1 = 351.479;
lat1 = +0.643;
radius1 = 68.;

# ~~~~~~~~~ File Paths ~~~~~~~~~~~~
parentDir = os.path.abspath(os.path.dirname(os.getcwd()))
dataDir = os.path.join(parentDir, 'Data/')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

class BubbleShortenedPicked():
    def __init__(self, wcs):
        # ~~~~~~~~ location of the regions ~~~~~~~~
        long1 = 351.479; lat1= +0.643; radius1 = 68.; name1 = 'G' + str(long1) + '+' + str(lat1); dist1 = 0; categ1 = 'Q'
        long2= 351.462; lat2= +0.556; radius2 = 68.; name2 = 'G' + str(long2) + '+' + str(lat2); dist2 = 0; categ2 = 'G'
        long3 = 351.424; lat3 = +0.650; radius3 = 99.; name3 = 'G' + str(long3) + '+' + str(lat3); dist3 = 0; categ3 = 'G'
        long4 = 351.420; lat4 = +0.637; radius4 = 52.; name4 = 'G' + str(long4) + '+' + str(lat4); dist4 = 0; categ4 = 'Q'
        long5 = 351.383; lat5 = +0.737; radius5 = 303.; name5 = 'G' + str(long5) + '+' + str(lat5); dist5 = 1.3; categ5 = 'K'
        long6 = 351.367; lat6= +0.640; radius6 = 107.; name6 = 'G' + str(long6) + '+' + str(lat6); dist6 = 1.3; categ6 = 'K'
        long7= 351.311; lat7= +0.663; radius7 = 119.; name7 = 'G' + str(long7) + '+' + str(lat7); dist7 = 1.3; categ7 = 'K'
        long8= 351.246; lat8 = +0.673; radius8 = 131.; name8 = 'G' + str(long8) + '+' + str(lat8); dist8 = 1.3; categ8 = 'K'
        long9 = 351.170; lat9 = +0.704; radius9 = 171.; name9 = 'G' + str(long9) + '+' + str(lat9); dist9 = 0; categ9 = 'K'
        long10 = 351.153; lat10 = +0.623; radius10 = 68.; name10 = 'G' + str(long10) + '+' + str(lat10); dist10 = 0; categ10 = 'Q'
        long11 = 351.348; lat11 = +0.593; radius11 = 125.; name11 = 'G' + str(long11) + '+' + str(lat11); dist11 = 0; categ11 = 'G'
        long12 = 351.130; lat12 = +0.449; radius12 = 549.; name12 = 'G' + str(long12) + '+' + str(lat12); dist12 = 0; categ12 = 'K'
        long13 = 350.995; lat13 = +0.654; radius13 = 428.; name13 = 'G' + str(long13) + '+' + str(lat13); dist13 = 0; categ13 = 'K'
        RegionLong = []
        RegionLat = []
        self.RegionRA = []
        self.RegionDec = []
        self.RegionDec = []
        self.RegionX = []
        self.RegionY = []
        self.RegionRadius = []
        self.RegionThickinArcSec = []
        self.RegionLable = []
        self.RegionDist = []
        self.RegionCategory = []
        self.RegionRadiusinRADeg = []
        self.RegionRadiusinRAarcSec=[]
        # ~~~~~~~~~~~~~~~~ Rest will be after reading the fits file, can change to a class instead ~~~~~~~~~

        #~~~~~~~~~~~~~~~~ Simpson et al 2012 bubbles ~~~~~~~~~~~~~~~~
        SimpsonFilename = os.path.join(dataDir, 'SimpsonBubblesFinal.txt')
        Long, Lat, InXDiam, InYDiam, OutXDiam, Reff, Thickness, Eccentricity, ElipsPositionAnlge, Hit, PositionDispersion, Flag, Ra, Dec = np.loadtxt(
            SimpsonFilename, usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), skiprows=2, unpack=True)

        # ~~~~~~~~~~~~~~~~ For the HII regions and etc. ~~~~~~~~~~~~~~~
        for index in range(13):
            RegionLong.append(locals()['long{}'.format(index+1)])
            RegionLat.append(locals()['lat{}'.format(index+1)])
            self.RegionRadiusinRADeg.append(locals()['radius{}'.format(index + 1)]/3600)
            self.RegionRadius.append(locals()['radius{}'.format(index + 1)]/3.999999)# 10 pixels = 39.9999999999996 arcsec
            gc = SkyCoord(l=RegionLong[index] * u.degree, b=RegionLat[index] * u.degree, frame='galactic')
            tempRa = (gc.fk5.ra * u.degree).value
            tempDec = (gc.fk5.dec * u.degree).value
            self.RegionRA.append(tempRa)
            self.RegionDec.append(tempDec)
            tempX, tempY = wcs.wcs_world2pix(tempRa, tempDec, 0) #Python indeces start from 0
            self.RegionX.append(tempX)
            self.RegionY.append(tempY)
            self.RegionLable.append(locals()['name{}'.format(index+1)])
            self.RegionDist.append(locals()['dist{}'.format(index+1)])
            # self.RegionCategory.append(locals()['categ{}'.format(index+1)])
            self.RegionThickinArcSec.append(0.)

        # ~~~~~ For the bubbles Simpson et al. 2012 ~~~~~~
        for index in range(len(Ra)):
            self.RegionRA.append(Ra[index])
            self.RegionDec.append(Dec[index])
            self.RegionRadiusinRADeg.append(Reff[index]/60)
            self.RegionRadiusinRAarcSec.append(Reff[index]*60.)
            self.RegionRadius.append((Reff[index]*60.)/3.999999)#This is in pixel
            self.RegionThickinArcSec.append(Thickness[index]*60.)
            tempX, tempY = wcs.wcs_world2pix(Ra[index], Dec[index], 0)
            self.RegionX.append(tempX)
            self.RegionY.append(tempY)
            self.RegionLable.append(index+1)
            self.RegionDist.append(0)
