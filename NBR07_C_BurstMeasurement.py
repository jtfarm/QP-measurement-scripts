# -*- coding: utf-8 -*-

# from VISAdrivers.continuousAlazar import ADC
# from tools.datatools import bin2csv
# import numpy as np
import os
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from scipy.signal import windows, convolve
import time
# import Labber
import subprocess

pathToExe = r'C:/Users/LFL/lflPython/AlazarDrivers/CS_Average/x64/Release/ATS9371_CS_Average.exe'

nHours = 6
# nMinutesDelay = 30
#numberTraces = nHours*60//nMinutesDelay
numberTraces = 1
acquisitionLength_sec = nHours*3600
origRateMHz = 300
# avgTime = 3e-6
sampleRateMHz = 1 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.
DAsetting = 20
LOfrequency = 4.2727  #GHz


StringForFlux = r'{}GHz_DA{}_SR{}MHz'.format(LOfrequency,DAsetting,sampleRateMHz)
path = r"F:\NBR07_bursts\NBR07_Mar29_2022_Burst\{}\\".format(StringForFlux)


if not os.path.exists(path):
    os.makedirs(path)

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
    f.write("LO frequency: "+str(LOfrequency) + " GHz")


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