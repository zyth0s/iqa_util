#!/usr/bin/env python

"""
Converts basin surface files between Promolden (modular) and
ChemInt format.

Promolden                                      ChemInt
(modular)                                         ^^
   ^^                                             ||
   ||                                             ||
   vv                     surfconv.py             vv
<name>.wfn.surf-txtxxxx  ------------->  <name>.<gms/adf/fhi>.surfxxxx.h5
"""

name = "ch4"
natoms = 5
ESprogram = "gms" # gms/adf/fhi/g09

from __future__ import print_function

import h5py as h5

def readlines_generator(fname):
     with open(fname) as file:
         for i in file:
             yield i


for iatom in range(1,natoms+1):
  fname = "{}.wfn.surf-txt{0:04d}".format(name,iatom)
  print("Converting " + fname + "...")
  surfile = readlines_generator(fname)

  line = next(surfile)
  nintersecs = int(line.split()[0])

  f = h5.File("{}.{}.surf{:04d}.h5".format(name,ESprogram,iatom), "w")
  ct       = f.create_dataset("ct"       , (1,nintersecs , ) , dtype = 'f8')
  st       = f.create_dataset("st"       , (1,nintersecs , ) , dtype = 'f8')
  cp       = f.create_dataset("cp"       , (1,nintersecs , ) , dtype = 'f8')
  sp       = f.create_dataset("sp"       , (1,nintersecs , ) , dtype = 'f8')
  angw     = f.create_dataset("angw"     , (1,nintersecs , ) , dtype = 'f8')
  nlimsurf = f.create_dataset("nlimsurf" , (1,nintersecs , ) , dtype = 'i8')
  rlimsurf = f.create_dataset("rlimsurf" , (10,nintersecs , ) , dtype = 'f8')

  for line in surfile:
    if "cos(theta)" in line:
      for i in range(nintersecs):
        line = next(surfile)
        ct[0,i] = float(line.split()[0])
        st[0,i] = float(line.split()[1])
        cp[0,i] = float(line.split()[2])
        sp[0,i] = float(line.split()[3])
        angw[0,i] = float(line.split()[4])
        nlimsurf[0,i] = 1
        rlimsurf[0,i] = float(line.split()[5])
        #print(next(surfile))

