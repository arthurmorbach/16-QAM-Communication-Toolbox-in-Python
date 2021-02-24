Channel
=======
Module containing the functions that cover the channel characteristics.

Additive White Gaussian Noise
*****************************
.. automodule:: channel
   :members: AWGN
   :undoc-members:
   :show-inheritance:
   :noindex:

.. code-block:: python

   import numpy as np

   def AWGN(IFsig, SNR):
      dP = np.zeros(len(IFsig))
      P = 0

      for i in range(len(IFsig)):
         dP[i] = abs(IFsig[i])**2
         P = P + dP[i]

      P = P/len(IFsig)
      gamma = 10**(SNR/10)
      N0 = P/gamma
      n = ((N0/2)**(0.5))*np.random.standard_normal(len(IFsig))
      IF_n = np.zeros(len(IFsig))

      for i in range(len(IFsig)):
         IF_n[i] = IFsig[i] + n[i]

      return(IF_n)