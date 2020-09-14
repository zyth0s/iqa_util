#!/usr/bin/env python

# Compare the second order reduced density  of Promolden and ChemInt

from __future__ import print_function


def readlines_generator(fname):
     with open(fname) as file:
         for i in file:
             yield i

fortranfname =                  "Fe2CO9_PBE_ccpvdz.rhout"
cxxfname     = "../pmd++_modsurf/Fe2CO9_PBE_ccpvdz.rhoout"

fortranfile = readlines_generator(fortranfname)
cxxfile     = readlines_generator(cxxfname)

# Search starting point in each file
while True:
  line = next(fortranfile)
  if "#  e1  e1  e2  e2" in line:
    break
while True:
  line = next(cxxfile)
  if "#  e1   e1   e2   e2" in line:
    next(cxxfile)
    break

linecxx = next(cxxfile)
linefortran = next(fortranfile)
while not "End computing RDM" in linecxx or not "Total elapsed time" in linefortran:
  epsilon = 1e-4
  e1_f, e1_c = int(linefortran.split()[0]), int(linecxx.split()[1])
  e2_f, e2_c = int(linefortran.split()[1]), int(linecxx.split()[2])
  e3_f, e3_c = int(linefortran.split()[2]), int(linecxx.split()[3])
  e4_f, e4_c = int(linefortran.split()[3]), int(linecxx.split()[4])
  error_text = "Error in {:4d} {:4d} {:4d} {:4d} [Fotran numbering] ->".format(e1_f, e2_f, e3_f, e4_f)
  assert e1_f == e1_c, "{} e1: [fortran] {} != {} [C++]".format(error_text,e1_f, e1_c)
  assert e2_f == e2_c, "{} e2: [fortran] {} != {} [C++]".format(error_text,e2_f, e2_c)
  assert e3_f == e3_c, "{} e3: [fortran] {} != {} [C++]".format(error_text,e3_f, e3_c)
  assert e4_f == e4_c, "{} e4: [fortran] {} != {} [C++]".format(error_text,e4_f, e4_c)
  tot_f,   tot_c = float(linefortran.split()[4]), float(linecxx.split()[5])
  xc_f,     xc_c = float(linefortran.split()[5]), float(linecxx.split()[6])
  x_f,       x_c = float(linefortran.split()[6]), float(linecxx.split()[7])
  coul_f, coul_c = float(linefortran.split()[7]), float(linecxx.split()[8])
  #assert abs(tot_f - tot_c )   < epsilon, "{} tot:  [fortran] {} != {} [C++]".format(error_text,tot_f, tot_c)
  #assert abs(xc_f - xc_c )     < epsilon, "{} xc:   [fortran] {} != {} [C++]".format(error_text,xc_f, xc_c)
  #assert abs(x_f - x_c )       < epsilon, "{} x:    [fortran] {} != {} [C++]".format(error_text,x_f, x_c)
  #assert abs(coul_f - coul_c ) < epsilon, "{} coul: [fortran] {} != {} [C++]".format(error_text,coul_f, coul_c)
  if abs(tot_f - tot_c )   > epsilon: 
    print("{} tot:  [fortran] {} != {} [C++]".format(error_text,tot_f, tot_c) )
  if abs(xc_f - xc_c )     > epsilon:
    print("{} xc:   [fortran] {} != {} [C++]".format(error_text,xc_f, xc_c) )
  if abs(x_f - x_c )       > epsilon:
    print("{} x:    [fortran] {} != {} [C++]".format(error_text,x_f, x_c) )
  if abs(coul_f - coul_c ) > epsilon:
    print("{} coul: [fortran] {} != {} [C++]".format(error_text,coul_f, coul_c) )
  #
  #print(error_text)
  #if abs(tot_f - tot_c )   > epsilon: 
  #  print("tot:  [fortran] {} != {} [C++]".format(tot_f, tot_c) )
  #if abs(xc_f - xc_c )     > epsilon:
  #  print("xc:   [fortran] {} != {} [C++]".format(xc_f, xc_c) )
  #if abs(x_f - x_c )       > epsilon:
  #  print("x:    [fortran] {} != {} [C++]".format(x_f, x_c) )
  #if abs(coul_f - coul_c ) > epsilon:
  #  print("coul: [fortran] {} != {} [C++]".format(coul_f, coul_c) )
  linecxx = next(cxxfile)
  linefortran = next(fortranfile)

