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

SPATH = r"F:\NBR07_AutoTune_FluxSweep\\"
figpath = SPATH + r"\Figures\\"
if not os.path.exists(figpath):
    os.makedirs(figpath)
    

logging.basicConfig(filename=os.path.join(SPATH,f'MEASUREMENTLOG_{time.strftime("%Y%m%d_%H%M%S")}.log'),filemode='w',level=logging.INFO)

pathToExe = r'C:/Users/LFL/lflPython/AlazarDrivers/CS_Average/x64/Release/ATS9371_CS_Average.exe'

client = Labber.connectToServer(timeout=None)
LO = client.connectToInstrument('Rohde&Schwarz RF Source',
                                dict(interface='TCPIP',address='192.168.1.128',startup='Get config'))
SMU = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='23',startup='Get config'))
SMUj = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='18',startup='Get config'))
VNA = client.connectToInstrument('Agilent Network Analyzer E5071B',dict(interface='GPIB',address='17',startup='Get config'))
SA = client.connectToInstrument('HP Spectrum Analyzer',dict(interface='GPIB',address='30',startup='Get config'))
PUMP = client.connectToInstrument('SignalCore SC5511A Signal Generator',dict(name='10002F25',startup='Get config'))
DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(interface='USB',address='24679',startup='Get config'))

LO.setValue('Output',False)
PUMP.setValue('Output status',False)
DA.setValue('Attenuation',20)
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


##############
# JPA tuneup optimization function
############
def tuneup(x):
    # record 10 traces to average
    p = Powers[x[0]]
    fl = Fluxes[x[1]]
    fp = Freqs[x[2]]
    # print(f'setting flux to {fl*1e3} mA')
    SMUj.setValue('Source current',fl,rate=0.0005)
    PUMP.setValue('Power',p)
    PUMP.setValue('Frequency',fp*1e9)
    
    sleep(0.5)
    dSA = SA.getValue('Signal')
    pSig = np.mean(dBm2Watt(dSA['y']))
    LO.setValue('Output',False)
    
    sleep(0.5)
    dSA = SA.getValue('Signal')
    pNoise = np.mean(dBm2Watt(dSA['y']))
    LO.setValue('Output',True)
    
    snr = Watt2dBm(pSig)-Watt2dBm(pNoise)
    print(f'P = {p:.5f} | f = {fp:.5f} | F = {fl:.8f} | snr = {snr:.4f}')
    logging.info(f'P = {p:.5f} | f = {fp:.5f} | F = {fl:.8f} | snr = {snr:.4f}')
    return -snr

def calb(res):
    if -res.fun - snrOFF > 8:
        return True
    else:
        return False

######################
# bounds for JPA tuneup, Freqs changes dynamically in loop below
#####################
Powers = np.round(np.arange(0,3,0.01),decimals=2)
Fluxes = np.round(np.arange(-26.9,-25.9,0.001)*1e-3,decimals=6)

#nHours = 12
#nMinutesDelay = 30
#numberTraces = nHours*60//nMinutesDelay
numberTraces = 1
acquisitionLength_sec = 5
origRateMHz = 300 # this is hard coded into the .exe file. Should always be 300
# avgTime = 3e-6
sampleRateMHz = 10 # Note this should always be a integer factor of origRateMHz. Such as 15 x 20 = 300.

T = 30


# LO.setValue('Frequency',LOfrequency*1e9)
PHIS = np.arange(0.5,0.3,-0.02)
NPOINTS = len(PHIS)
averageTimeCycle = 0
highI = curFunc(0.22) # mA
lowI = curFunc(0.42) # mA
lfVNA = Labber.createLogFile_ForData('NBR07_autoTuneFluxSweep_Jun29_VNA',
                                     [{'name':'VNA - S21','complex':True,'vector':True,'x_name':'Frequency','x_unit':'Hz'}],
                                     step_channels=[{'name':'Phi','values':PHIS,'unit':'flux quanta'}])
lfSA = Labber.createLogFile_ForData('NBR07_autoTuneFluxSweep_Jun29_SA',
                                     [{'name':'SA - Signal','unit':'dBm','vector':True,'x_name':'Frequency','x_unit':'Hz'}],
                                     step_channels=[{'name':'PUMP','values':[False,True]},
                                                    {'name':'Phi','values':PHIS,'unit':'flux quanta'}])
tuneGuess = [1.5,-26.7e-3,8.56]
shiftGuess = 4e-4
for count,ph in enumerate(PHIS):
    now = perf_counter()
    # get the current and freq mapping from flux fits
    estimated_freq = freqFunc(ph)*1e9
    I = curFunc(ph) # mA
    logging.info(f'\n\nStarting phi = {ph:.3f}.\nEstimated frequency is {estimated_freq*1e-9:.6f} GHz.\nCurrent on coil is {I:.6f} mA.')
    # do I tune the resonator high or low?
    if ph <= 0.35:
        SMU.setValue('Source current',lowI*1e-3,rate=0.0002)
    else:
        SMU.setValue('Source current',highI*1e-3,rate=0.0002)
    
    # with resonator tuned away, take VNA data for background
    VNA.setValue('Range type','Center - Span')
    VNA.setValue('Center frequency', estimated_freq)
    VNA.setValue('Span',8e6)
    VNA.setValue('Output enabled',True)
    VNA.setValue('Output power',-20)
    VNA.setValue('# of averages',900)
    sleep(1)
    VNA.setValue('Trigger',True)
    dBG = VNA.getValue('S21')
    xBG = np.arange(dBG['t0'],dBG['t0']+dBG['shape'][0]*dBG['dt'],dBG['dt'])
    zBG = dBG['y']
    
    # bring the resonator to the correct flux
    SMU.setValue('Source current',I*1e-3,rate=0.0002)
    
    # take VNA trace of resonator
    sleep(1)
    VNA.setValue('Trigger',True)
    dData = VNA.getValue('S21')
    zData = dData['y']/zBG
    td2 = Labber.getTraceDict(zData,x0=xBG[0],x1=xBG[-1])
    lfVNA.addEntry({'VNA - S21':td2})
    
    try:
        X = xBG*1e-9
        fguessInd = np.argmin(zData.real)
        fguess = X[fguessInd]
        shiftGuess = 1.5e-9*(qp.f_n_phi(ph,0)-qp.f_n_phi(ph,1))
        pars,cov = curve_fit(sumLor,X,zData.real,p0 = [0.8,0.5,0.1,fguess,shiftGuess,0.00025],bounds=([0.1,0.002,1e-6,fguess-100e-6,0.7*shiftGuess,0.0001],[2,2,2,fguess+50e-6,1.3*shiftGuess,0.0005]))
        plt.plot(X,zData.real,label='data')
        plt.plot(X,sumLor(X,*pars),label='fit')
        plt.axvline(pars[3])
        plt.axvline(pars[3]-pars[4])
        plt.axvline(pars[3]-2*pars[4])
        plt.legend()
        plt.title(f'PHI = {ph:.3f} | shift = {pars[4]*1e6:.1f} kHz')
        plt.ylabel('S21 - real')
        plt.xlabel('Frequency [GHz]')
        plt.savefig(figpath+f'PHI_{ph*1000:3.0f}.png')
        plt.show()
        plt.close()
        f0 = pars[3]
        fd = pars[3]-pars[4]
    except:
        X = xBG*1e-9
        fguess = np.argmin(zData.real)
        try:
            logging.warning('Unable to fit shift. Trying to fit a single resonance instead.')
            pars,cov = curve_fit(Lor,X,zData.real,p0 = [0.8,X[fguess],0.00025],bounds=([0.2,4.23,0.00005],[2,4.32,0.0005]))
            plt.plot(X,zData.real,label='data')
            plt.plot(X,Lor(X,*pars),label='fit')
            plt.axvline(pars[1])
            plt.legend()
            resshift = (qp.f_n_phi(ph,0)-qp.f_n_phi(ph,1))
            f0 = pars[1]
            fd = pars[1] - 1.5e-9*resshift
            pars[2] = fd
            plt.axvline(fd)
            plt.title(f'PHI = {ph:.3f} | shift = {resshift*1e-3:.1f} kHz')
            plt.ylabel('S21 - real')
            plt.xlabel('Frequency [GHz]')
            plt.savefig(figpath+f'PHI_{ph*1000:3.0f}.png')
            plt.show()
            plt.close()
            
        except:
            logging.warning('Unable to fit resonance. Using the minimum of real part as estimated resonance.')
            plt.plot(X,zData.real,label='data')
            pars = np.array([1.5,X[fguess],0.00025])
            cov = np.ones((len(pars),len(pars)))*np.infty
            plt.plot(X,Lor(X,*pars),label='unfitted guess')
            plt.axvline(pars[1])
            plt.legend()
            resshift = (qp.f_n_phi(ph,0)-qp.f_n_phi(ph,1))
            f0 = pars[1]
            fd = pars[1] - 1.5e-9*resshift
            pars[2] = fd
            plt.axvline(fd)
            plt.title(f'*PHI = {ph:.3f} | shift = {resshift*1e-3:.1f} kHz*')
            plt.ylabel('S21 - real')
            plt.xlabel('Frequency [GHz]')
            plt.savefig(figpath+f'PHI_{ph*1000:3.0f}.png')
            plt.show()
            plt.close()
            
        
    # save a copy of the fit parameters from VNA
    with open(figpath+f'PHI_{ph*1000:3.0f}_fit.pkl','wb') as pkfile:
        pickle.dump([pars,cov],pkfile)
    logging.info(f'The 0 QP resonance is found at {f0:.6f} GHz.\nWe are driving at {fd:.6f} GHz, which is {(f0-fd)*1e6:.1f} kHz downshifted.')
    # configure SA
    SA.setValue('Center frequency',fd*1e9)
    SA.setValue('Span',0.5e6)
    SA.setValue('IF bandwidth',3e3)
    
    # turn off VNA, turn on sig gen at fd, turn on pump at 2(fd+3MHz)
    VNA.setValue('Output enabled',False)
    LO.setValue('Frequency',fd*1e9)
    LO.setValue('Output',True)
    PUMP.setValue('Frequency',2*(fd*1e9 + 10e6))
        # take a SA trace for saving
    SAsig = np.zeros(401)
    for _ in range(30):
        sleep(1)
        dSA = SA.getValue('Signal')
        SAsig += dBm2Watt(dSA['y'])
    xSA = np.arange(dSA['t0'],dSA['t0']+dSA['shape'][0]*dSA['dt'],dSA['dt'])
    SAsig /= 30
    td1 = Labber.getTraceDict(Watt2dBm(SAsig),x0=xSA[0],x1=xSA[-1])
    lfSA.addEntry({'SA - Signal':td1})
    
    
    max_ind=np.argmax(SAsig)
    max_val=np.max(SAsig)
    mask = np.logical_or(xSA < xSA[max_ind]-10e3, xSA > xSA[max_ind]+10e3)
    noise=SAsig[mask]
    avg_noise=np.mean(noise)
    snrOFF_2 = Watt2dBm(max_val)-Watt2dBm(avg_noise)
    print(f'SNR without JPA is {snrOFF_2:.3f} dB')
    logging.info(f'SNR without JPA is {snrOFF_2:.3f} dB')
    
    
    
    
    SA.setValue('IF bandwidth',1e6)
    sleep(1)
    dSA = SA.getValue('Signal')
    pSig = np.mean(dBm2Watt(dSA['y']))
    LO.setValue('Output',False)
    
    sleep(1)
    dSA = SA.getValue('Signal')
    pNoise = np.mean(dBm2Watt(dSA['y']))
    LO.setValue('Output',True)
    
    snrOFF = Watt2dBm(pSig)-Watt2dBm(pNoise)
    print(f'snrOFF = {snrOFF:.3f} dB')
    
    PUMP.setValue('Output status',True)
    # optimize JPA using minimize func
    # boundsP = Bounds(-3,7)
    # boundsF = Bounds(-27.3e-3,-26.5e-3)
    # optRes = minimize(tuneup,tuneGuess,method='L-BFGS-B',bounds=[(-3,7),(-27.3e-3,-26.5e-3)])
    # optRes = dual_annealing(tuneup, [(-3,7),(-27.3e-3,-26.5e-3)], no_local_search=True)
    Freqs = np.round(np.arange(2*(fd + 0.01),2*(fd + 0.075),0.0001),decimals=4)
    SA.setValue('IF bandwidth',1e6)
    optRes = gp_minimize(tuneup, [(0,len(Powers)-1),(0,len(Fluxes)-1),(0,len(Freqs)-1)],n_calls=200,n_initial_points=30,callback=calb)
    # print(f'Successful optimization = {optRes["success"]}\n')
    X = optRes.x
    print(f'best values are P={Powers[X[0]]},fl={Fluxes[X[1]]},fp={Freqs[X[2]]} which gives snr = {-optRes.fun:.3f}')
    logging.info(f'best values are P={Powers[X[0]]},fl={Fluxes[X[1]]},fp={Freqs[X[2]]} which gives snr = {-optRes.fun:.3f}')
    SMUj.setValue('Source current',Fluxes[X[1]],rate=0.0005)
    PUMP.setValue('Power',Powers[X[0]])
    PUMP.setValue('Frequency',Freqs[X[2]]*1e9)
    
    # take a SA trace for saving
    SA.setValue('IF bandwidth',3e3)
    sleep(2)
    SAsig = np.zeros(401)
    for _ in range(30):
        sleep(1)
        dSA = SA.getValue('Signal')
        SAsig += dBm2Watt(dSA['y'])
    xSA = np.arange(dSA['t0'],dSA['t0']+dSA['shape'][0]*dSA['dt'],dSA['dt'])
    SAsig /= 30
    td3 = Labber.getTraceDict(Watt2dBm(SAsig),x0=xSA[0],x1=xSA[-1])
    lfSA.addEntry({'SA - Signal':td3})
    max_ind=np.argmax(SAsig)
    max_val=np.max(SAsig)
    mask = np.logical_or(xSA < xSA[max_ind]-10e3, xSA > xSA[max_ind]+10e3)
    noise=SAsig[mask]
    avg_noise=np.mean(noise)
    snrON = Watt2dBm(max_val)-Watt2dBm(avg_noise)
    print(f'SNR with JPA is {snrON:.3f} dB | increases by {snrON - snrOFF_2} dB')
    logging.info(f'SNR with JPA is {snrON:.3f} dB | increases by {snrON - snrOFF_2} dB')
    
    
    # Take Alazar data power sweep
    logging.info(f'\nStarting Alazar acquisition at phi = {ph:.3f} while driving at {fd:.6f} GHz')
    for ds in np.arange(0,41,2):
        DAsetting = ds
        DA.setValue('Attenuation',ds)
        StringForFlux = r'{:3.0f}flux\DA{:2.0f}_SR{:d}MHz'.format(ph*1000,DAsetting,sampleRateMHz)
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
            f.write("PHI: "+str(ph))
        
        
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
    PUMP.setValue('Output status',False)
    DA.setValue('Attenuation',20)
    timeCycle = perf_counter()-now
    print(f'This step took {timeCycle:.6f} seconds')
    averageTimeCycle += (timeCycle-averageTimeCycle)/(count+1)
    Nremain = NPOINTS - (count+1)
    print(f'There are {Nremain} steps cycles left, estimated time remaining = {averageTimeCycle*Nremain/3600:.3f} hours')
    logging.info(f'This step took {timeCycle:.6f} seconds.\nThere are {Nremain} steps cycles left, estimated time remaining = {averageTimeCycle*Nremain/3600:.3f} hours')