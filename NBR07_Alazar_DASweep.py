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
import fitTools.quasiparticleFunctions as qp
import matplotlib.pyplot as plt
import matplotlib.colors as mplc

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
sampleRateMHz = 1 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.
#DAsetting = 6
# LOfrequency = 4.27244 # GHz
T = 70
Vcurrent = 2.265
phi = 0.44
f0qp = 4.279654
fshift = 0.000187
# LOfrequency = np.arange(f0qp,f0qp-1.5*fshift,-1*fshift)
LOfrequency = np.array([f0qp, f0qp-1*fshift, f0qp-1.5*fshift])
# LOfrequency = np.array([4.27210, 4.27199])

power = np.arange(45,5,-1)
# x = []
# y = []
for f in LOfrequency:
    LO.setValue('Frequency',f*1e9)
    qppeak = np.round((f0qp-f)/fshift,decimals=1)
    
    freq = f
    for p in power:
        DA.setValue('Attenuation',p)
        
        StringForFlux = r'{:.6}GHz_DA{}'.format(freq,str(p))
        path = r"G:\Shared drives\LFL\Projects\Quasiparticles\NBR07_fluxSweep_Apr2022\{}mK\VictorCurrent_{}mA_phi_{}\{}\\".format(T,Vcurrent,phi,StringForFlux)

        
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
            f.write("DA setting: "+str(p) + " dB\n")
            f.write("Temperature: "+str(T)+' mK\n')
            f.write("Victor current: "+str(Vcurrent)+' mA\n')
            f.write('phi: '+str(phi)+'\n')
            f.write('QP peak:' + str(qppeak))
        
        
        Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))
        
        print(Creturn)

    # data = qp.loadAlazarData(savefile)
    # data = qp.GaussianConvolution(data, 14, 10)
    # data = qp.uint16_to_mV(data)
    # # ax = qp.plotComplexHist(data[0],data[1])
    # # ax.set_title(f'{I*1e3:.1f} mA')
    # # plt.savefig(figpath+f'\\{I*1e6}uA.png')
    # # plt.show();
    # # plt.close();
    # x.append(np.mean(data[0]))
    # y.append(np.mean(data[1]))
    
    # fig,ax = plt.subplots()
    # hi = plt.hist2d(data[0],data[1],bins=(80,80),cmap=plt.get_cmap('Greys'))
    # # hi = np.histogram2d(data[0],data[1],bins=(200,200))
    # xc = (hi[1][:-1]+hi[1][1:])/2
    # yc = (hi[2][:-1]+hi[2][1:])/2
    # # guess = [60000,np.mean(data[0]),np.mean(data[1]),1,1,0]
    # # xx,yy,amps,means,varis = qp.fitGaussian(hi,guess)
    # # # f = gaussianMix(heights,widths,means)
    # # qp.make_ellipses2(means,varis,ax,['red'])
    # ax.set_title(f'{p} dB')
    # plt.savefig(figpath+f'\\{p}dB.png')
    # plt.show()
    # plt.close()
    # print(varis)

# fig,ax = plt.subplots(1,1,'none',figsize=[4,3],constrained_layout=True)
# colors = plt.get_cmap('gist_rainbow', len(x))
# norm = mplc.Normalize(vmin=0, vmax=len(x))
# sm = plt.cm.ScalarMappable(cmap=colors, norm=norm)
# sm.set_array([])
# fig.colorbar(sm, aspect=60)
# plt.plot(x,y)
# for i,(xx,yy) in enumerate(zip(x,y)):
#     plt.scatter(xx,yy,c=colors(i))
# plt.savefig(figpath+r'summary.png')
# LO.setValue('Frequency',LOfrequency*1e9)

# stringdesc = f"{int(LOfrequency*1000)}"


# # StringForFlux = r'{}GHz_DA{}_SR{}MHz'.format(LOfrequency,DAsetting,sampleRateMHz)
# path = r"G:\Shared drives\LFL\Projects\Quasiparticles\TestOffsetNoise\\"
# figpath = r"G:\Shared drives\LFL\Projects\Quasiparticles\TestOffsetNoise\figures\\"

# if not os.path.exists(path):
#     os.makedirs(path)
# if not os.path.exists(figpath):
#     os.makedirs(figpath)

# timestamp = time.strftime("%Y%m%d_%H%M%S")
# savefile = path + '{}.bin'.format(stringdesc)

# samplesPerPoint = int(max(origRateMHz/sampleRateMHz,1))
# actualSampleRateMHz = origRateMHz/samplesPerPoint

# # write metadata to corresponding .txt file
# with open(savefile[0:-4] + ".txt",'w') as f:
#     from time import strftime
#     f.write(strftime("%c")+'\n')
#     f.write("Channels: " + 'AB' + '\n')
#     f.write("Acquisition duration: " + str(acquisitionLength_sec) + " seconds." + '\n')
#     f.write("Sample Rate MHz: " + str(actualSampleRateMHz) + '\n')
#     f.write("LO frequency: "+str(LOfrequency) + " GHz")

# # savefile = adc.startTriggeredCapture(acquisitionLength_sec,channel='AB',dataFilePath=savefile,returnfname=True,downsamplerate=sampleRateMHz*1e6)
# Creturn = subprocess.getoutput('"{}" {} {} "{}"'.format(pathToExe,int(acquisitionLength_sec),samplesPerPoint,savefile))

# print(Creturn)

# data = qp.loadAlazarData(savefile)
# data = qp.BoxcarDownsample(data,2e-6,sampleRateMHz*1e6)
# data = qp.uint16_to_mV(data)
# ax,hi = qp.plotComplexHist(data[0],data[1],bins=(80,80),returnHistData=True)
# xc = (hi[1][:-1]+hi[1][1:])/2
# yc = (hi[2][:-1]+hi[2][1:])/2
# guess = [1000,0,0,1,1,0]
# xx,yy,amps,means,varis = qp.fitGaussian(hi,guess)
# # f = gaussianMix(heights,widths,means)
# qp.make_ellipses2(means,varis,ax,['red'])
# print(f'\n\n{stringdesc} gives {varis}\n\n')
# plt.title(f'{stringdesc}')
# plt.savefig(figpath+f'{stringdesc}_IQhist.png')
# plt.show()




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