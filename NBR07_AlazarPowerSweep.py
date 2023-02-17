# -*- coding: utf-8 -*-

import os
import time
import Labber
import subprocess
import numpy as np
import fitTools.quasiparticleFunctions as qp
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from scipy.optimize import curve_fit, minimize, Bounds, dual_annealing
from fitTools.utilities import dBm2Watt, Watt2dBm
from time import perf_counter, sleep
from skopt import gp_minimize
import pickle
import logging

SPATH = r"G:\Shared drives\Quasiparticles_2\NBR07_Sept2022\Test_DA_Power_Sweep_No_TWPA"
figpath = SPATH + r"\Figures\\"
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
#logging.basicConfig(filename=SPATH+'MEASUREMENTLOG.log')

pathToExe = r'C:/Users/LFL/lflPython/AlazarDrivers/CS_Average/x64/Release/ATS9371_CS_Average.exe'

client = Labber.connectToServer(timeout=None)

LO = client.connectToInstrument('Rohde&Schwarz RF Source',
                                dict(interface='TCPIP',address='192.168.1.128',startup='Get config'))
SMU = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='23',startup='Get config'))
DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(interface='USB',address='24679',startup='Get config'))

#################
# Measurement parameters
#################
numberTraces = 1
acquisitionLength_sec = 5
origRateMHz = 300
sampleRateMHz = 10 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.

phi = 0.47
T = 30

fd = 4.27256
f0 = 4.2727

I = SMU.getValue('Source current')
LO.setValue('Frequency',fd*1e9)
LO.setValue('Output',True)

now = perf_counter()
for ds in np.arange(0,40,2):
    DAsetting = ds
    DA.setValue('Attenuation',ds)
    StringForFlux = r'{:3.0f}flux_DA{:2.0f}_SR{:d}MHz'.format(phi*1000,DAsetting,sampleRateMHz)
    path = os.path.join(SPATH,"{}".format(StringForFlux))
    

    if not os.path.exists(path):
        os.makedirs(path)

    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    savefile = os.path.join(path,'NBR07_{}.bin'.format(timestamp))
    
    samplesPerPoint = int(max(origRateMHz/sampleRateMHz,1))
    actualSampleRateMHz = origRateMHz/samplesPerPoint
    
    # write metadata to corresponding .txt file
    with open(savefile[0:-4] + ".txt",'w') as f:
        from time import strftime
        f.write(strftime("%c")+'\n')
        f.write("Channels: " + 'AB' + '\n')
        f.write("Acquisition duration: " + str(acquisitionLength_sec) + " seconds." + '\n')
        f.write("Sample Rate MHz: " + str(actualSampleRateMHz) + '\n')
        f.write("LO frequency: "+str(fd) + " GHz\n")
        f.write("flux bias: "+str(I) + " mA\n")
        f.write("DA setting: "+str(DAsetting) + " dB\n")
        f.write("Temperature: "+str(T)+' mK\n')
        f.write("PHI: "+str(phi))
    
    
    Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))
    
    logging.info(Creturn)
    
    data = qp.loadAlazarData(savefile)
    data = qp.uint16_to_mV(data)
    ax = qp.plotComplexHist(data[0],data[1])
    ax.set_title(f'PHI = {phi:.4f}')
    plt.savefig(figpath+f'\\{phi*1000}.png')
    plt.show()
    plt.close()

timeCycle = perf_counter()-now
print(f'This step took {timeCycle:.6f} seconds')
