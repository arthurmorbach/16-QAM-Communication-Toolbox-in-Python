BER Analysis
===================
Analyses the BER to each SNR value.

.. automodule:: BER_SNR_analysis
   :members:
   :undoc-members:
   :show-inheritance:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   import QAMsys as qam

   def SNR_BER_analysis(): 
      plt.close('all')

      SNR = 0
      n_tests_per_snr = 50
      max_snr = 20

      BER_mean_acc = []

      for i in range(max_snr):
         SNR = SNR + 1
         BER = np.zeros(n_tests_per_snr)
         BER_acc = 0
         
         for i in range(len(BER)):
               BER[i] = qam.QAMsys(SNR,0)
               BER_acc = BER_acc + BER[i]
         
         BER_mean = BER_acc/len(BER)
         BER_mean_acc.append(BER_mean)
         
         plt.figure(1)
         plt.scatter(SNR, BER_mean)
         plt.xscale('linear')
         plt.yscale('log')

      plt.xlabel('SNR (dB)')
      plt.ylabel('BER')
      plt.grid()
   
.. image:: images/SNR.png
   :width: 600
   :align: center