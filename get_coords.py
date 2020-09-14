# coding: utf-8

# Extract the coordinates (Bohr) of the atoms in the unit cell
# from a DGrid WFN file.

from __future__ import print_function

from sys import exit
import numpy as np

#from mendeleev import element

from select_file_dialog import gui_fname

def readlines_generator(fname):
     with open(fname) as file:
         for i in file:
             yield i

def get_coords():
  dgrid_wfn = gui_fname(filter="dgrid_wfn")
  if not dgrid_wfn:
    #print("Canceled. No file selected")
    exit()
              
  wfnf = readlines_generator(dgrid_wfn)

  xyz = np.zeros((0,3))
  for line in wfnf:
    if "coordinates" in line:
      next(wfnf)
      next(wfnf)
      next(wfnf)
      line = next(wfnf)
      while not "+------" in line:
        elSymbol = line.split()[1]
        #xyz.append([float(line.split()[3]),float(line.split()[4]),float(line.split()[5])])
        xyz_list = [float(line.split()[3]),float(line.split()[4]),float(line.split()[5])]
        xyz = np.append(xyz, [xyz_list], axis=0)
        print("{0:3s}  {1}".format(elSymbol, xyz[-1]))
        line = next(wfnf)
  return xyz, dgrid_wfn

