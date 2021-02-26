# -*- coding: utf-8 -*-
import numpy as np
import commpy as cp
from math import pi

def data_gen(N, data_sync=0):
    """Generates an array of data. If the synchronization bits are not informed a default sequence will be used.

    Parameters
    ----------
    N : int
        Number of bits.
    data_sync : 1D array of ints. 
        Synchronization bits.
    
    Returns
    -------
    data : 1D array of ints
        Pseudo randomic data with synchronization bits.
    """    
    if data_sync == 0:
        data_sync_osc = []
        for i in range(176):
            data_sync_osc.append(1)
        data_sync_symb = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]  
        data_sync = np.concatenate((data_sync_osc, data_sync_symb), axis=None)
    data_r = np.random.rand(N - len(data_sync))
    data_r[np.where(data_r >= 0.5)] = 1
    data_r[np.where(data_r < 0.5)] = 0
    data = np.concatenate((data_sync, data_r), axis=None)
    return(data)


def slicer(data): 
    """It separates the even bits in the In-phase (I) vector, and the odd bits in the Quadrature (Q) vector.

    Parameters
    ----------
    data : 1D array of ints
        Array that will be divided.
        
    Returns
    -------
    dataI : 1D array of ints
        Array with the even position bits.
    dataQ : 1D array of ints
        Array with the odd position bits.
    """    
  
    dataI = data[slice(0, len(data), 2)]
    dataQ = data[slice(1, len(data), 2)]
    return(dataI, dataQ)


def mapper_16QAM(QAM16, data):
    """Uses the input data to index the array with the 16QAM amplitudes and phases.
    
    Parameters
    ----------
    QAM16 : 1D array of floats
        Array that will be indexed.
    data : 1D array of ints
        Data that will index the QAM16 array.
        
    Returns
    -------
    dataMapped : 1D array of floats
        Array with the mapped data.
    """    
    map0 = 2*data[slice(0, len(data), 2)] + data[slice(1, len(data), 2)]
    map0 = list(map(int, map0))
    dataMapped = []
    for i in range(len(map0)):
        dataMapped.append(QAM16[map0[i]])
    return(dataMapped)

def upsampler(Ns, K, symbols):
    """Increases the symbol samples by adding zeros between each value.

    Parameters
    ----------
    Ns : int
        Number of symbols.
    K : int
        Up-sampler factor.
    symbols : 1D array of floats
        Symbols array.

    Returns
    -------
    up : 1D array of floats
        Array with the upsampled data.
    """   
    up = np.zeros(Ns*K)
    up[slice(0, len(up), K)] = symbols
    return(up)


def shaping_filter(upsampler, Ns, alpha, Fif, Fs):
    """To give the symbols a shape, a convolution is made between the upsampled symbols and the impulse response of a square root raised cosine filter.
    It also arranges the information in a defined frequency spectrum and is projected in a way to reduce the intersymbol interference.

    Parameters
    ----------
    upsampler : int
        Upsampled symbols.
    Ns : int
        Number of symbols.
    alpha : float
        Roll-off factor. Numbers between 0-1.
    Fif : float
        Intermediary frequency.
    Fs : float
        Sampling frequency

    Returns
    -------
    shaped_signal : 1D array of floats
        Signal convoluted with a SRRC filter impulse response.
    x_axis : 1D array of floats
        Data domain.
    y_response : 1D array of floats
        Array with the amplitudes varying regarding the domain.
    """
    [x_axis, y_response] = cp.rrcosfilter(Ns, alpha, 2/Fif, Fs)
    shaped_signal = np.convolve(upsampler, y_response, 'full')
    return(shaped_signal, x_axis, y_response)


def oscillator(start, stop, step, frequency, phase=0):
    """Generates the carrier signal.

    Parameters
    ----------
    start : number
        Start of interval. The interval includes this value.
    stop : number
        End of interval. The interval does not include this value.
    step : number
        The distance between two adjacent values.  
    frequency : float
        Frequency of oscillation in Hz.
    phase : float, optional
        Phase of the sinusuidal wave. The default phase is 0.

    Returns
    -------
    Osc : 1D array of floats
        Amplitude values of the sinusoidal wave.
    time : 1D array of floats
        Data domain.
    """
    t = np.arange(start, stop, step)
    Osc = np.sin(2*pi*frequency*t + phase)
    return(Osc, t)


def mixer(signal, carrier):
    """It is a pointwise product function. In this application it mixes a signal with the carrier, shifting the frequency spectrum to a defined intermediary frequency.

    Parameters
    ----------
    signal : 1D array of floats
        Signal that will be mixed.
    carrier : 1D array of floats
        Carrier signal.

    Returns
    -------
    mix : 1D array of floats
        Mixed signal.
    """
    mix = []
    for i in range(len(signal)):
        mix.append(signal[i]*carrier[i])
    return(mix)

def combiner(signal_I, signal_Q):
    """It's a pointwase sum, combining the modulated signals in quadrature.

    Parameters
    ----------
    signalI : 1D array of floats
        In-phase symbols modulated.
    signalQ : 1D array of floats
        Quadrature symbols modulated.

    Returns
    -------
    combined_sig : 1D array of floats
        Quadrature signal.
    """
    combined_sig = []
    for i in range(len(signal_I)):
        combined_sig.append(signal_I[i] + signal_Q[i])
    return(combined_sig)
