# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:24:38 2022

@author: lfl
"""

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

client = Labber.connectToServer(timeout=None)
LO = client.connectToInstrument('Rohde&Schwarz RF Source',
                                dict(interface='TCPIP',address='192.168.1.128',startup='Get config'))
SA = client.connectToInstrument('HP Spectrum Analyzer',dict(interface='GPIB',address='30',startup='Get config'))
PUMP = client.connectToInstrument('SignalCore SC5511A Signal Generator',dict(name='10002F25',startup='Get config'))
DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(interface='USB',address='24679',startup='Get config'))


PUMP.setValue('Output status',False)
DA.setValue('Attenuation',20)

##############
# JPA tuneup optimization function
############
def tuneup(x):
    # record 10 traces to average
    p = Powers[x[0]]
    fp = Freqs[x[1]]
    # print(f'setting flux to {fl*1e3} mA')
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
    print(f'P = {p:.5f} | f = {fp:.5f} | snr = {snr:.4f}')
    return -snr

def calb(res):
    if -res.fun - snrOFF > 8:
        return True
    else:
        return False

######################
# bounds for JPA tuneup, Freqs changes dynamically in loop below
#####################
Powers = np.round(np.arange(3.5,6.5,0.01),decimals=2)
Freqs = np.round(np.arange(8.02,8.08,0.001),decimals=3)

fd = 4.27
# configure SA
SA.setValue('Center frequency',fd*1e9)
SA.setValue('Span',0.5e6)
SA.setValue('IF bandwidth',1e6)
    
LO.setValue('Frequency',fd*1e9)
LO.setValue('Output',True)
PUMP.setValue('Frequency',8.068e9)

sleep(2)
dSA = SA.getValue('Signal')
pSig = np.mean(dBm2Watt(dSA['y']))
LO.setValue('Output',False)

sleep(2)
dSA = SA.getValue('Signal')
pNoise = np.mean(dBm2Watt(dSA['y']))
LO.setValue('Output',True)

snrOFF = Watt2dBm(pSig)-Watt2dBm(pNoise)
print(f'snrOFF = {snrOFF:.3f} dB')

PUMP.setValue('Output status',True)


optRes = gp_minimize(tuneup, [(0,len(Powers)-1),(0,len(Freqs)-1)],n_calls=50,n_initial_points=10)
# print(f'Successful optimization = {optRes["success"]}\n')
X = optRes.x
print(f'best values are P={Powers[X[0]]},fp={Freqs[X[1]]} which gives snr = {-optRes.fun:.3f}')
PUMP.setValue('Power',Powers[X[0]])
PUMP.setValue('Frequency',Freqs[X[1]]*1e9)

SA.setValue('IF bandwidth',3e3)
PUMP.setValue('Output status',False)
SAsig = np.zeros(401)
for _ in range(30):
    sleep(1)
    dSA = SA.getValue('Signal')
    SAsig += dBm2Watt(dSA['y'])
xSA = np.arange(dSA['t0'],dSA['t0']+dSA['shape'][0]*dSA['dt'],dSA['dt'])
SAsig /= 30
max_ind=np.argmax(SAsig)
max_val=np.max(SAsig)
mask = np.logical_or(xSA < xSA[max_ind]-10e3, xSA > xSA[max_ind]+10e3)
noise=SAsig[mask]
avg_noise=np.mean(noise)
snrOFF_2 = Watt2dBm(max_val)-Watt2dBm(avg_noise)
print(f'SNR without JPA is {snrOFF_2:.3f} dB')


PUMP.setValue('Output status',True)
sleep(1)
SAsig = np.zeros(401)
for _ in range(30):
    sleep(1)
    dSA = SA.getValue('Signal')
    SAsig += dBm2Watt(dSA['y'])
xSA = np.arange(dSA['t0'],dSA['t0']+dSA['shape'][0]*dSA['dt'],dSA['dt'])
SAsig /= 30
max_ind=np.argmax(SAsig)
max_val=np.max(SAsig)
mask = np.logical_or(xSA < xSA[max_ind]-10e3, xSA > xSA[max_ind]+10e3)
noise=SAsig[mask]
avg_noise=np.mean(noise)
snrON = Watt2dBm(max_val)-Watt2dBm(avg_noise)
print(f'SNR with JPA is {snrON:.3f} dB | increases by {snrON - snrOFF_2} dB')