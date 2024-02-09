# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(path)
import example as qam
    
def SNR_BER_analysis():
    """Analyze the BER for each SNR value."""
    plt.close('all')

    n_tests_per_snr = 50
    max_snr = 20

    SNR_values = np.arange(1, max_snr + 1)
    BER_mean_acc = []

    for SNR in SNR_values:
        print("SNR:", SNR)
        BER = np.zeros(n_tests_per_snr)
        
        for j in range(len(BER)):
            BER[j] = qam.QAMsys(SNR, 0)
        
        BER_mean = np.mean(BER)
        BER_mean_acc.append(BER_mean)

    plt.figure(1)
    plt.scatter(SNR_values, BER_mean_acc)
    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel('SNR (dB)')
    plt.ylabel('BER')
    plt.grid()
    plt.show()
    
SNR_BER_analysis()