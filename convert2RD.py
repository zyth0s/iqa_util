#!/usr/bin/env python

# Convert 2RDM from Promolden to ChemInt format

from __future__ import print_function

import sys

nmopair = 4005
epsrho = 1e-12


def readlines_generator(fname):
     with open(fname) as file:
         for i in file:
             yield i

fortranfname =                  "Fe2CO9_PBE_ccpvdz.wfn.2rdm_txt"
#cxxfname     = "../pmd++_modsurf/Fe2CO9_PBE_ccpvdz.g09.2rdm"
targetfname  = "../pmd++_modsurf/Fe2CO9_PBE_ccpvdz.g09.2rdm_conv"

fortranfile = readlines_generator(fortranfname)
#cxxfile = readlines_generator(cxxfname)

# Search starting point in files
linefortran = next(fortranfile)
#linecxx = next(cxxfile)

with open(targetfname, "w") as targetf:
  targetf.write("{}\n".format(nmopair))

  while True: 

      try:
        linefortran = next(fortranfile)
      except:
        print("Finished or error")
        break

      #linecxx = next(cxxfile)
      #ei = int(linecxx.split()[0])
      #ej = int(linecxx.split()[1])
      print(linefortran)

      ei = int(linefortran.split()[0])
      ej = int(linefortran.split()[1])
      tot_f  = float(linefortran.split()[2])
      xc_f   = float(linefortran.split()[3])
      x_f    = float(linefortran.split()[4])
      coul_f = float(linefortran.split()[5])

      targetf.write("{:4d} {:4d} {:16.8f} {:16.8f} {:16.8f} {:16.8f}\n".format( 
                     ei-1, ej-1, tot_f, xc_f, x_f, coul_f))


