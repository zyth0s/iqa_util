#!/usr/bin/env python
# coding: utf-8
"""  Example output:
Errors in the first layer surface

Atom  Max error  Avg error  Min error
----  ---------  ---------  ---------
   1    0.01137    0.00057    0.00000
   2    0.02089    0.00045    0.00000
   3    0.02089    0.00045    0.00000
   4    0.02089    0.00045    0.00000
   5    0.02089    0.00045    0.00000
"""

import h5py as h5
import numpy as np

natoms = 4 

print("Errors in the first layer surface\n")
print("Atom  Max error  Avg error  Min error" )
print("----  ---------  ---------  ---------" )
for method in ["bs3_default","ark3","sa4","ark4","z4","dp5","ck5","f5","verner6","f8"]:
  print("{}".format(method))
  for i in range(1,natoms+1):
    #f = h5.File("benchmark_surface/ch4.adf.surf{:04d}.h5".format(i))
    f = h5.File("f8/NH3_HF_CC-PVTZ.gms.surf{:04d}.h5".format(i))
    #print(f.keys())
    bench_first_surf = f["rlimsurf"][0] # first intersection
    f = h5.File("{:s}/NH3_HF_CC-PVTZ.gms.surf{:04d}.h5".format(method,i))
    first_surf = f["rlimsurf"][0]
    error_max = np.linalg.norm(first_surf - bench_first_surf, ord=np.inf)
    error_min = np.linalg.norm(first_surf - bench_first_surf, ord=-np.inf)
    error_avg = abs(np.average(first_surf - bench_first_surf) )
    print("{0:4d}  {1:13.9f}  {2:13.9f}  {3:13.9f}".format( i, error_max, error_avg, error_min) )

