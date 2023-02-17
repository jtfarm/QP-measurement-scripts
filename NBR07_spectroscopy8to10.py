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

T = 30
ph = 0.415
estimated_freq = 4.2802968*1e9
fd = 4.280087
shift = estimated_freq - fd*1e9

SPATH = r"E:\\NBR07_spectroscopy_{:3.0f}\\".format(ph*1000)
# SPATH = r"E:\\driveNextMode_{:3.0f}\\".format(ph*1000)
figpath = SPATH + r"\Figures\\"
if not os.path.exists(figpath):
    os.makedirs(figpath)
    

logging.basicConfig(filename=SPATH+f'MEASUREMENTLOG_{time.strftime("%Y%m%d_%H%M%S")}.log',filemode='w',level=logging.INFO)

pathToExe = r'C:/Users/LFL/lflPython/AlazarDrivers/CS_Average/x64/Release/ATS9371_CS_Average.exe'

client = Labber.connectToServer(timeout=None)
LO = client.connectToInstrument('Rohde&Schwarz RF Source',
                                dict(interface='TCPIP',address='192.168.1.128',startup='Get config'))
SMU = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='23',startup='Get config'))
VNA = client.connectToInstrument('Agilent Network Analyzer E5071B',dict(interface='GPIB',address='17',startup='Get config'))
SA = client.connectToInstrument('HP Spectrum Analyzer',dict(interface='GPIB',address='30',startup='Get config'))
DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(interface='USB',address='24679',startup='Get config'))
LO2 = client.connectToInstrument('SignalCore SC5511A Signal Generator',dict(name='10002A05',startup='Get config'))

LO.setValue('Output',False)

################
## THESE ARE ONLY VALID FOR THIS RUN!!! 6/14/22
###############
def curFunc(phi):
    return phi*13.1703 - 3.5704 # returns the current in mA for NBR07. ONLY VALID THIS RUN.
def freqFunc(fl):
    w0 = 4.300379
    q0 = 0.0118226
    return w0*(1+q0*np.sin(np.pi*fl/2)*np.arctanh(np.sin(np.pi*fl/2))/(1-np.sin(np.pi*fl/2)*np.arctanh(np.sin(np.pi*fl/2))))**(-0.5)

###############
# fit multiple modes in VNA data
###############
def sumLor(f,A1,A2,A3,f1,shift,Gamma):
    return 1 - A1/(1+(2*(f-f1)/Gamma)**2) - A2/(1+(2*(f-(f1-shift))/Gamma)**2) - A3/(1+(2*(f-(f1-2*shift))/Gamma)**2)
def Lor(f,A1,f1,Gamma):
    return 1 - A1/(1+(2*(f-f1)/Gamma)**2)


#nHours = 12
#nMinutesDelay = 30
#numberTraces = nHours*60//nMinutesDelay
numberTraces = 1
acquisitionLength_sec = 30
origRateMHz = 300    # this is hard coded into the .exe file. Should always be 300
# avgTime = 3e-6
sampleRateMHz = 2 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.




# LO.setValue('Frequency',LOfrequency*1e9)
# PHIS = np.arange(0.5,0.3,-0.001)
powersetting = np.round(np.arange(-25,10,3),decimals=0)
attens = np.round(np.arange(30,0,-3),decimals=0)
powers = np.concatenate((powersetting[0]-attens,powersetting))
# powers = np.round(np.arange(-25,10,0.5),decimals=1)
# PHIS = np.arange(0.348,0.3,-0.001)
freqs = np.round(np.arange(10.2,11,0.025),decimals=3)
NPOINTS = len(powers)*len(freqs)
averageTimeCycle = 0
highI = curFunc(0.22) # mA
lowI = curFunc(0.42) # mA
# lfVNA = Labber.LogFile(r"G:\Shared drives\LFL\Labber\Data\2022\07\Data_0714\NBR07_autoFluxSweep_{}mK".format(int(T)))
lfVNA = Labber.createLogFile_ForData('NBR07_spectroscopy_{:3.0f}'.format(ph*1000),
                                      [{'name':'VNA - S21','complex':True,'vector':True,'x_name':'Frequency','x_unit':'Hz'}],
                                      step_channels=[{'name':'Power','values':powers,'unit':'dBm'},{'name':'Frequency','values':freqs,'unit':'GHz'}])

I = curFunc(ph)

# do I tune the resonator high or low?
LO2.setValue('Output status',False)
if ph <= 0.35:
    SMU.setValue('Source current',lowI*1e-3,rate=0.0002)
else:
    SMU.setValue('Source current',highI*1e-3,rate=0.0002)


# Take VNA trace for background
VNA.setValue('Range type','Center - Span')
VNA.setValue('Center frequency', estimated_freq)
VNA.setValue('Span',3e6)
VNA.setValue('Output enabled',True)
VNA.setValue('Output power',-20)
VNA.setValue('# of averages',400)
sleep(1)
VNA.setValue('Trigger',True)
dBG = VNA.getValue('S21')
xBG = np.arange(dBG['t0'],dBG['t0']+dBG['shape'][0]*dBG['dt'],dBG['dt'])
zBG = dBG['y']

# bring resonator to ph
LO2.setValue('Output status',False)
SMU.setValue('Source current',I*1e-3,rate=0.0001)
        
# set up the DA and sig gen
LO2.setValue('Power',-25)
DA.setValue('Attenuation',30)

for count2,freq in enumerate(freqs):
    
    # set the spectroscopy frequency
    LO2.setValue('Frequency', freq*1e9)
    
    for count,p in enumerate(powers):
        now = perf_counter()
        # set up the drive power
        if p <= -25:
            atten = -1*(p + 25)
            DA.setValue('Attenuation',atten)
        else:
            LO2.setValue('Power',p)
        # LO2.setValue('Power',p)
        
        # # get the current and freq mapping from flux fits
        # estimated_freq = freqFunc(ph)*1e9
        # I = curFunc(ph) # mA
        # logging.info(f'\n\nStarting phi = {ph:.3f}.\nEstimated frequency is {estimated_freq*1e-9:.6f} GHz.\nCurrent on coil is {I:.6f} mA.')
        
        
        # # take a quick VNA trace to center on resonance
        # SMU.setValue('Source current',I*1e-3,rate=0.0002)
        
        # turn on the drive
        LO2.setValue('Output status',True)
        
        # VNA.setValue('Range type','Center - Span')
        # VNA.setValue('Center frequency', estimated_freq)
        # VNA.setValue('Span',10e6)
        # VNA.setValue('Output enabled',True)
        # VNA.setValue('Output power',-15)
        # VNA.setValue('# of averages',100)
        # sleep(1)
        # VNA.setValue('Trigger',True)
        # dF = VNA.getValue('S21')
        # xFind = np.arange(dF['t0'],dF['t0']+dF['shape'][0]*dF['dt'],dF['dt'])
        # zFind = dF['y']
        # fcInd = np.argmin(np.abs(zFind))
        # fcenter = xFind[fcInd]
        
        # do I tune the resonator high or low?
        # LO2.setValue('Output status',False)
        # if ph <= 0.35:
        #     SMU.setValue('Source current',lowI*1e-3,rate=0.0002)
        # else:
        #     SMU.setValue('Source current',highI*1e-3,rate=0.0002)
        
        # with resonator tuned away, take VNA data for background
        # LO2.setValue('Output status',True)
        # VNA.setValue('Range type','Center - Span')
        # VNA.setValue('Center frequency', fcenter)
        # VNA.setValue('Span',3e6)
        # VNA.setValue('Output enabled',True)
        # VNA.setValue('Output power',-20)
        # VNA.setValue('# of averages',100)
        # sleep(1)
        # VNA.setValue('Trigger',True)
        # dBG = VNA.getValue('S21')
        # xBG = np.arange(dBG['t0'],dBG['t0']+dBG['shape'][0]*dBG['dt'],dBG['dt'])
        # zBG = dBG['y']
        
        # bring the resonator to the correct flux
        # LO2.setValue('Output status',False)
        # SMU.setValue('Source current',I*1e-3,rate=0.0001)
        
        # take VNA trace of resonator
        LO2.setValue('Output status',True)
        VNA.setValue('Output enabled',True)
        VNA.setValue('Output power',-35)
        VNA.setValue('# of averages',999)
        sleep(1)
        VNA.setValue('Trigger',True)
        dData = VNA.getValue('S21')
        zData = dData['y']/zBG
        td2 = Labber.getTraceDict(zData,x0=xBG[0],x1=xBG[-1])
        lfVNA.addEntry({'VNA - S21':td2})
    
        
        # turn off VNA, turn on sig gen at fd
        VNA.setValue('Output enabled',False)
        LO.setValue('Frequency',fd*1e9)
        LO.setValue('Output',True)
        
        
        
        # Take Alazar data power sweep
        logging.info(f'\nStarting Alazar acquisition at drive power = {p:.1f}')
        StringForFlux = r'{}GHz\{:3.0f}dBm\SR{:d}MHz'.format(freq,p*10,sampleRateMHz)
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
            f.write("f0: " + str(estimated_freq)+'\n')
            f.write("flux bias: "+str(I) + " mA\n")
            f.write("Temperature: "+str(T)+' mK\n')
            f.write("PHI: "+str(ph)+'\n')
            f.write("clearingPower: " + str(p) +'\n')
            f.write("clearingFreq: " + str(freq) +'\n')
        
        
        Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))
        
        logging.info(Creturn)
        
        # data = qp.loadAlazarData(savefile)
        # data = qp.uint16_to_mV(data)
        # ax = qp.plotComplexHist(data[0],data[1])
        # ax.set_title(f'PHI = {ph:.4f}')
        # plt.savefig(figpath+f'\\{ph*1000}.png')
        # plt.show()
        # plt.close()
        
        # Turrn off JPA pump and LO, reset DA to 20
        LO.setValue('Output',False)
        timeCycle = perf_counter()-now
        print(f'This step took {timeCycle:.6f} seconds')
        averageTimeCycle += (timeCycle-averageTimeCycle)/(count+1)
        Nremain = NPOINTS - (count+1) - (count2)*len(powers)
        print(f'There are {Nremain} steps cycles left, estimated time remaining = {averageTimeCycle*Nremain/3600:.3f} hours')
        logging.info(f'This step took {timeCycle:.6f} seconds.\nThere are {Nremain} steps cycles left, estimated time remaining = {averageTimeCycle*Nremain/3600:.3f} hours')