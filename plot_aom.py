#!/usr/bin/env python
# coding: utf-8
"""    Plot the Atomic Overlap Matrix (AOM) as a 2D image
"""

import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

name = "Al_fcc"
program = "fhi"
natoms = 1

for i in range(1,natoms+1):
  filename = "{0:}.{1:}{2:04d}.aom".format(name,program,i)
  print(filename)
  f = h5.File(filename)
  aom = f["dataset"][:]
  print(aom)
  print("Trace = {:5.3f}".format(np.trace(aom)))
  #plt.matshow(aom)
  plt.matshow(aom)
  plt.colorbar()
  plt.show()

