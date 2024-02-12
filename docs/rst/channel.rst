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

   def AWGN(sig, SNR, K):
      sig_abs2 = [abs(s)**2 for s in sig]
      P = (K * sum(sig_abs2)) / len(sig_abs2)
      gamma = 10**(SNR/10)
      N0 = P / gamma
      n = np.sqrt(N0 / 2) * np.random.standard_normal(len(sig))
      sig_n = sig + n
      return sig_n