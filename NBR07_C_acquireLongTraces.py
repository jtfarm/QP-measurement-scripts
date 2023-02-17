# -*- coding: utf-8 -*-

# from VISAdrivers.continuousAlazar import ADC
# from tools.datatools import bin2csv
# import numpy as np
import os
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from scipy.signal import windows, convolve
import time
import Labber
import subprocess
import numpy as np

pathToExe = r'C:/Users/LFL/lflPython/AlazarDrivers/CS_Average/x64/Release/ATS9371_CS_Average.exe'

client = Labber.connectToServer()
LO = client.connectToInstrument('Rohde&Schwarz RF Source',
                                dict(interface='TCPIP',address='192.168.1.128'))


DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator', dict(interface='USB', address='24679'))
DA.startInstrument()

#nHours = 12
#nMinutesDelay = 30
#numberTraces = nHours*60//nMinutesDelay
numberTraces = 1
acquisitionLength_sec = 5
origRateMHz = 300
# avgTime = 3e-6
sampleRateMHz = 10 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.
#DAsetting = 6
# LOfrequency = 4.27244 # GHz

LOfrequency = np.array([4.27267, 4.27256, 4.27244, 4.27233, 4.27222, 4.27210, 4.27199])
# LOfrequency = np.array([4.27267, 4.27256, 4.27244, 4.27222, 4.27210, 4.27199])
# LOfrequency = np.array([4.27233])
# power = np.array([15])
power = np.arange(50,0,-1)
# x = []
# y = []
for f in LOfrequency:
    LO.setValue('Frequency',f*1e9)
    
    freq = f
    for p in power:
        DA.setValue('Attenuation',p)
        
        StringForFlux = r'{}GHz_DA{}'.format(freq,str(p))
        path = r"G:\Shared drives\LFL\Projects\Quasiparticles\NBR07_Mar2022_Temperature\75mK\{}\\".format(StringForFlux)
        # path = r"E:\Quasiparticles\NBR07_Mar2022_Burst\{}\\".format(StringForFlux)

        
        # StringForFlux = r'{}GHz_DA{}'.format(str(f),str(p))
        # path = r"E:\NBR07_Jan25_2022\0p47flux_Freq_Power_Sweep\{}\\".format(StringForFlux)
        # figpath = r"G:\Shared drives\LFL\Projects\Quasiparticles\NBR07_Jan20_2022\PowerSweep\Figures\\"
        
        if not os.path.exists(path):
            os.makedirs(path)
        # if not os.path.exists(figpath):
        #     os.makedirs(figpath)
        
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        savefile = path + 'NBR07_{}.bin'.format(timestamp)
        
        samplesPerPoint = int(max(origRateMHz/sampleRateMHz,1))
        actualSampleRateMHz = origRateMHz/samplesPerPoint
        
        # write metadata to corresponding .txt file
        with open(savefile[0:-4] + ".txt",'w') as f:
            from time import strftime
            f.write(strftime("%c")+'\n')
            f.write("Channels: " + 'AB' + '\n')
            f.write("Acquisition duration: " + str(acquisitionLength_sec) + " seconds." + '\n')
            f.write("Sample Rate MHz: " + str(actualSampleRateMHz) + '\n')
            f.write("LO frequency: "+ str(freq) + " GHz" + '\n')
            f.write("DA setting: "+str(p) + " dB")
        
        
        Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))
        
        print(Creturn)
# adc = ADC()
# adc.configureClock(MS_s = origRateMHz)
# adc.configureTrigger(source='INT')

# for i in range(numberTraces):
#     now = time.perf_counter()
    
#     # acquire data
#     print('Starting acquisition {}'.format(i))
#     timestamp = time.strftime("%Y%m%d_%H%M%S")
#     savefile = path + 'NBR07_{}.bin'.format(timestamp)
#     # write metadata to corresponding .txt file
#     with open(savefile[0:-4] + ".txt",'w') as f:
#         from time import strftime
#         f.write(strftime("%c")+'\n')
#         f.write("Channels: " + 'AB' + '\n')
#         f.write("Acquisition duration: " + str(acquisitionLength_sec) + " seconds." + '\n')
#         f.write("Sample Rate MHz: " + str(actualSampleRateMHz) + '\n')


#     Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))
    
#     time.sleep(nMinutesDelay*60 - (time.perf_counter() - now))
    #sleep(60)