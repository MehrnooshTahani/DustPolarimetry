__author__ = 'Mehrnoosh Tahani'

import os, sys, subprocess, shlex, glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import statistics as st

class RadialParam:
    def __init__(self, Dir, cx, cy, rx, tx, i_CutOff, pi_cutOff):
        binWidth = 20
        # ~~~~~ I, Q, and U File names ~~~~
        iext_nameCl = 'iextN633MO808YYNN.fit'
        iextInsideAST_name = 'iextN633MO808YYNNinsideAST.fit'
        Qext_nameCl = 'qextN633MO808YYNN.fit'
        Uext_nameCl = 'uextN633MO808YYNN.fit'
        Piext_nameCl = 'piextN633MO808YYNN.fit'
        PFext_nameCl = 'pfextN633MO808YYNN.fit'
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
        # ~~~~~~~~~~ Header and WCS ~~~~~~~
        ifits_nameCl = os.path.join(Dir, iext_nameCl)
        ihduCl = fits.open(ifits_nameCl)
        ifhdCl = ihduCl[0].header
        ifbdCl = ihduCl[1].data
        wcs = WCS(ifhdCl).celestial
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
        # ~~~~~~~~~ Cutoff values ~~~~~~~~~
        I_cutoff = i_CutOff
        PI_cutoff = pi_cutOff
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
        # ~~~ I, Q, U arrays from files ~~~ (TURN THESE TO A CLASS LATER)
        # Stokes I Inside AST (the units are different)
        ifitsInsideAST = os.path.join(Dir, iextInsideAST_name)
        ihduInsideAST = fits.open(ifitsInsideAST)
        isdtInsideAST = ihduInsideAST['PRIMARY',1].data[0]
        iedtInsideAST = ihduInsideAST['VARIANCE',1].data[0]
        # Stokes I
        isdtCl = ihduCl['PRIMARY',1].data[0]
        iedtCl = ihduCl['VARIANCE',1].data[0]
        # Stokes Q
        Qfits_nameCl = os.path.join(Dir, Qext_nameCl)
        QhduCl = fits.open(Qfits_nameCl)
        QsdtCl = QhduCl['PRIMARY', 1].data[0]
        QedtCl = QhduCl['VARIANCE', 1].data[0]
        # Stokes U
        Ufits_nameCl = os.path.join(Dir, Uext_nameCl)
        UhduCl = fits.open(Ufits_nameCl)
        UsdtCl = UhduCl['PRIMARY', 1].data[0]
        UedtCl = UhduCl['VARIANCE', 1].data[0]
        # Polarized I
        Pifits_nameCl = os.path.join(Dir, Piext_nameCl)
        PihduCl = fits.open(Pifits_nameCl)
        PisdtCl = PihduCl['PRIMARY',1].data[0]
        PiedtCl = PihduCl['VARIANCE',1].data[0]
        # Polarized Fraction
        PFfits_nameCl = os.path.join(Dir, PFext_nameCl)
        PFhduCl = fits.open(PFfits_nameCl)
        PFsdtCl = PFhduCl['PRIMARY',1].data[0]
        PFedtCl = PFhduCl['VARIANCE',1].data[0]
        # Stokes I again for Making (so that it remains untouched while isdt is masked)
        ifitsForMasking = os.path.join(Dir, iext_nameCl)
        ihduForMasking = fits.open(ifitsForMasking)
        isdtForMasking = ihduForMasking['PRIMARY',1].data[0]
        iedtForMasking = ihduForMasking['VARIANCE',1].data[0]
        # Polarized I again for Making (so that it remains untouched while isdt is masked)
        PifitsForMasking_name = os.path.join(Dir, Piext_nameCl)
        PihduForMasking = fits.open(PifitsForMasking_name)
        PisdtForMasking = PihduForMasking['PRIMARY',1].data[0]
        PiedtForMasking = PihduForMasking['VARIANCE',1].data[0]
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
        # ~~Defining Class Parameters ~~~~~.
        self.PFsdtCl = PFsdtCl
        self.QsdtCl = QsdtCl
        self.UsdtCl = UsdtCl
        self.isdtCl = isdtCl
        PI_r = []
        self.Pf_rAvg = []
        self.Q_rAVG = []
        self.U_rAVG = []
        I_r =[]
        self.N_r=[]
        self.radius =[]
        self.Q_rStdev = []
        self.U_rStdev = []
        self.AngleAvg = []
        self.AngleStdev = []
        self.angleBin = np.arange(-80, 100, binWidth)
        self.pixelinEachAngleBin = np.zeros([self.angleBin.shape[0], int(rx)])
        self.cosAVG = []
        self.cosStdev = []
        self.cosAvgAfterAveg = []
        self.cosAvgAfterAveg_Stdev = []
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.
        # ~~loop for pixels in rings around bubble ~~~.
        for radiusIndex in range(int(rx)):
            Pf=0.
            Q_rNormalized = 0 ; Q_rTemp = 0
            U_rNormalized = 0 ; U_rTemp = 0
            locals()['PolAngInRing{}'.format(radiusIndex + 1)] = []
            locals()['PFValuesInRing{}'.format(radiusIndex+1)] = []
            locals()['PIValuesInRing{}'.format(radiusIndex+1)] = []
            locals()['Q_rValuesInRing{}'.format(radiusIndex+1)] = []
            locals()['U_rValuesInRing{}'.format(radiusIndex+1)] = []
            N = 0.
            # ~~ Go in loop  for each pixel to find the rings ~~~
            for xIndex in range(len(self.isdtCl[0,:])):
                for yIndex in range(len(self.isdtCl[:,0])):
                # This goes from 0 to 2xradius
                    if (xIndex - cx)**2  + (yIndex - cy)**2 <= (2*(radiusIndex+1))**2 and (xIndex - cx)**2  + (yIndex - cy)**2 > (2*(radiusIndex))**2:
                        if PFsdtCl[yIndex, xIndex] <= 100 and PFsdtCl[yIndex, xIndex] >= 0:
                            if (PisdtCl[yIndex, xIndex] / np.sqrt(PiedtCl)[yIndex, xIndex]) > PI_cutoff and (isdtCl[yIndex, xIndex] / np.sqrt(iedtCl)[yIndex, xIndex]) > I_cutoff:
                                print(self.QsdtCl[yIndex, xIndex])
                                '''IMPORTANT NOTE: the x is opposite of Ra, so it should be reverse to keep the 
                                convention also matches with Schmidt since x in their formula is ra '''
                                pfi = np.arctan2((cx - xIndex), (yIndex - cy))  # I checked and arctan and arctan2
                                # give exactly the same results
                                Q_rTemp = QsdtCl[yIndex, xIndex]*np.cos(2*pfi) + UsdtCl[yIndex, xIndex]*np.sin(2*pfi)
                                U_rTemp = -QsdtCl[yIndex, xIndex]*np.sin(2*pfi) + UsdtCl[yIndex, xIndex]*np.cos(2*pfi)
                                Angel_temp = np.degrees(0.5*np.arctan2(U_rTemp, Q_rTemp))
                                locals()['PolAngInRing{}'.format(radiusIndex + 1)].append(np.degrees(0.5*np.arctan2(U_rTemp, Q_rTemp)))
                                locals()['PFValuesInRing{}'.format(radiusIndex+1)].append(PFsdtCl[yIndex, xIndex])
                                locals()['PIValuesInRing{}'.format(radiusIndex+1)].append(PisdtCl[yIndex, xIndex])
                                locals()['Q_rValuesInRing{}'.format(radiusIndex + 1)].append(Q_rTemp)
                                locals()['U_rValuesInRing{}'.format(radiusIndex + 1)].append(U_rTemp)
                                N = N + 1
            # Finally before going to next ring, store information and find averages.
            if N > 0:
                self.N_r.append(N)
                self.radius.append((2*radiusIndex+1)*3.999999)#3.999999 is to convert it to arcsec for 8'' maps
                self.Pf_rAvg.append((st.mean(locals()['PFValuesInRing{}'.format(radiusIndex+1)])))
                self.PI_rAVG.append((st.mean(locals()['PIValuesInRing{}'.format(radiusIndex+1)])))
                self.Pf_rStdev.append((st.stdev(locals()['PFValuesInRing{}'.format(radiusIndex+1)]))/np.sqrt(N-1))
                self.PI_rStdev.append((st.stdev(locals()['PIValuesInRing{}'.format(radiusIndex+1)]))/np.sqrt(N-1))
                self.Q_rAVG.append((st.mean(locals()['Q_rValuesInRing{}'.format(radiusIndex+1)])))
                self.U_rAVG.append((st.mean(locals()['U_rValuesInRing{}'.format(radiusIndex+1)])))
                self.Q_rStdev.append((st.stdev(locals()['Q_rValuesInRing{}'.format(radiusIndex+1)]))/np.sqrt(N-1))
                self.U_rStdev.append((st.stdev(locals()['U_rValuesInRing{}'.format(radiusIndex+1)]))/np.sqrt(N-1))
                self.AngleAvg.append(st.mean(locals()['PolAngInRing{}'.format(radiusIndex + 1)]))
                self.AngleStdev.append(st.stdev(locals()['PolAngInRing{}'.format(radiusIndex + 1)]))
                self.cosAVG.append(st.mean(locals()['PolCosInRing{}'.format(radiusIndex + 1)]))
                self.cosStdev.append((st.stdev(locals()['PolCosInRing{}'.format(radiusIndex + 1)]))/(np.sqrt(N-1)))
                AngleAfterAvged = 0.5 * np.arctan2((st.mean(locals()['U_rValuesInRing{}'.format(radiusIndex+1)])), (st.mean(locals()['Q_rValuesInRing{}'.format(radiusIndex+1)])))
                self.cosAvgAfterAveg.append(np.cos(AngleAfterAvged))
                Q_rAvgTemp = st.mean(locals()['Q_rValuesInRing{}'.format(radiusIndex + 1)])
                U_rAvgTemp = st.mean(locals()['U_rValuesInRing{}'.format(radiusIndex + 1)])
                Q_rStdevTemp = (st.stdev(locals()['Q_rValuesInRing{}'.format(radiusIndex + 1)])) / np.sqrt(N - 1)
                U_rStdevTemp = (st.stdev(locals()['U_rValuesInRing{}'.format(radiusIndex + 1)])) / np.sqrt(N - 1)
                # ~~~~~~~~~~~~ For error calc of Cos ~~~~~~
                Func = np.cos(AngleAfterAvged)
                A = (1/(1+(U_rAvgTemp/Q_rAvgTemp)**2) )* (1/Q_rAvgTemp)
                dfunc_dU = -0.5*np.sin(AngleAfterAvged)*A
                dfunc_dQ = 0.5*np.sin(AngleAfterAvged)*A * (U_rAvgTemp/Q_rAvgTemp)
                # delta_Fun = Func * (np.sqrt(((dfunc_dQ)*(Q_rStdevTemp))**2 + ((dfunc_dU)*(U_rStdevTemp))**2))
                delta_Fun = (np.sqrt(((dfunc_dQ)*(Q_rStdevTemp))**2 + ((dfunc_dU)*(U_rStdevTemp))**2))
                # ~~~~~~~~~~~~ For error calc of Cos ~~~~~~.
                self.cosAvgAfterAveg_Stdev.append(delta_Fun)
            for idx, binVal in enumerate(self.angleBin):
                for angelIdx, angelVal in enumerate(locals()['PolAngInRing{}'.format(radiusIndex + 1)]):
                    if (binVal - (binWidth/2)) <= angelVal < (binVal + (binWidth/2)):
                        self.pixelinEachAngleBin[idx, radiusIndex] = self.pixelinEachAngleBin[idx, radiusIndex]+1






