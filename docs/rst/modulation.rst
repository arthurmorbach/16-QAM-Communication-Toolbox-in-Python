Modulation
==========
Module containing all the functions used in the modulation process of a 16 QAM communication system.

Data Generator
**************

.. automodule:: modulation
   :members: data_gen
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   import numpy as np

   def data_gen(N, data_sync=0):
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

.. image:: images/data.png
   :width: 600
   :align: center

Slicer
**************

.. automodule:: modulation
   :members: slicer
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   def slicer(data): 
      dataI = data[slice(0, len(data), 2)]
      dataQ = data[slice(1, len(data), 2)]
      return(dataI, dataQ)

.. image:: images/slicer.png
   :width: 600
   :align: center

Mapper
**************

.. automodule:: modulation
   :members: mapper_16QAM
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   def mapper_16QAM(QAM16, data):
      map0 = 2*data[slice(0, len(data), 2)] + data[slice(1, len(data), 2)]
      map0 = list(map(int, map0))
      dataMapped = []
      for i in range(len(map0)):
         dataMapped.append(QAM16[map0[i]])
      return(dataMapped)

.. image:: images/mapper.png
   :width: 600
   :align: center

.. image:: images/constellation.png
   :width: 600
   :align: center

Upsampler
**************

.. automodule:: modulation
   :members: upsampler
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   import numpy as np

   def upsampler(Ns, K, symbols):  
      up = np.zeros(Ns*K)
      up[slice(0, len(up), K)] = symbols
      return(up)

.. image:: images/upsampler.png
   :width: 600
   :align: center

Shaping Filter
**************

.. automodule:: modulation
   :members: shaping_filter
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   import numpy as np
   import commpy as cp

   def shaping_filter(upsampler, Ns, alpha, Fif, Fs):
      [x_axis, y_response] = cp.rrcosfilter(Ns, alpha, 2/Fif, Fs)
      shaped_signal = np.convolve(upsampler, y_response, 'full')
      return(shaped_signal, x_axis, y_response)
   
.. image:: images/impulse_response.png
   :width: 600
   :align: center

.. image:: images/shaping_filter.png
   :width: 600
   :align: center

Oscillator
**************

.. automodule:: modulation
   :members: oscillator
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   import numpy as np
   from math import pi
   
   def oscillator(start, stop, step, frequency, phase=0):
      t = np.arange(start, stop, step)
      Osc = np.sin(2*pi*frequency*t + phase)
      return(Osc, t)

.. image:: images/oscillator.png
   :width: 600
   :align: center

Mixer
**************

.. automodule:: modulation
   :members: mixer
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   def mixer(signal, carrier):
      mix = []
      for i in range(len(signal)):
         mix.append(signal[i]*carrier[i])
      return(mix)

.. image:: images/mixer.png
   :width: 600
   :align: center

Combiner
**************

.. automodule:: modulation
   :members: combiner
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   def combiner(signal_I, signal_Q):
      combined_sig = []
      for i in range(len(signal_I)):
         combined_sig.append(signal_I[i] + signal_Q[i])
      return(combined_sig)
   
.. image:: images/combiner.png
   :width: 600
   :align: center