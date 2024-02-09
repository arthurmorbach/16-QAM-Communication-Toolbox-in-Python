# -*- coding: utf-8 -*-
import numpy as np

def AWGN(sig, SNR, K):
    """Adds noise to the signal based on the proposed SNR.

    Parameters
    ----------
    sig : 1D array of floats
        IF signal.
    SNR : float
        Signal to noise ratio in dB.

    Returns
    -------
    sig_n : 1D array of floats
        Signal with noise based on the SNR.
    """
    sig_abs2 = [abs(s)**2 for s in sig]
    P = (K * sum(sig_abs2)) / len(sig_abs2)
    gamma = 10**(SNR/10)
    N0 = P / gamma
    n = np.sqrt(N0 / 2) * np.random.standard_normal(len(sig))
    sig_n = sig + n
    return sig_n