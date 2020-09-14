#!/usr/bin/env python
# coding: utf-8

# TODO(dmc): improve handling of symmetry generated surfaces. Now it assumes
#            that the symmetric basins have coordinate inuc+j with j = 1,2,3,...
#            i.e. atoms with same Wyckoff position are contiguous and the first
#            is taken to generate the others.

import os
import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
import h5py as h5
import numpy as np

#from select_file_dialog import gui_fname
from get_coords import get_coords


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#surf_file = gui_fname(filter="surf")
#if surf_file == "":
#  print("Canceled. No file selected")
#  exit()

xyz, dgrid_wfn = get_coords()
natoms = len(xyz)

for inuc in range(natoms):
  _surf_file = "{0}.surf{1:04d}.h5".format(dgrid_wfn,inuc+1)

  surf_file = _surf_file
  for j in range(100): # FIXME(dmc): magic number
    if j != 0:
      surf_file = _surf_file + "_sym{0:04d}".format(j)

    if not os.path.isfile(surf_file):
      if not "sym" in surf_file:
        logging.warning("There is no surface file {0:s}".format(surf_file))
      continue
    if not os.access(surf_file, os.R_OK):
      logging.warning("File {0:s} is not readable".format(surf_file))
      continue

    f = h5.File(surf_file)

    npang = len(f["ct"][0])
    nsurfpoints = sum(f["nlimsurf"][0])

    cartsurf = np.zeros((nsurfpoints,3))

    rmins = np.min(f["rlimsurf"][0])

    # Insert the points corresponding to the first intersection in each ray
    for i in range(npang):
      st = f["st"][0][i]
      cp = f["cp"][0][i]
      sp = f["sp"][0][i]
      rlim = f["rlimsurf"][0][i]
      # Spherical to cartesian
      cartsurf[i,0] = rlim * st * cp
      cartsurf[i,1] = rlim * st * sp
      cartsurf[i,2] = rlim * f["ct"][0][i] 
      # Displace the origin from (0,0,0) to the position of atom inuc
      cartsurf[i,0] += xyz[inuc+j,0]; # FIXME(dmc): +j is for symmetric generated basins
      cartsurf[i,1] += xyz[inuc+j,1]; # FIXME(dmc): +j is for symmetric generated basins
      cartsurf[i,2] += xyz[inuc+j,2]; # FIXME(dmc): +j is for symmetric generated basins

    # If there are rays with more than one intersection
    NinsertedPoints = npang
    nintersecs = 0
    while nsurfpoints > NinsertedPoints:
      nintersecs += 1
      raysWithMoreIntersections = np.where( f["nlimsurf"][0] > nintersecs )
      #print(raysWithMoreIntersections)

      for ray in range(len(raysWithMoreIntersections[0])):
        rayindex = raysWithMoreIntersections[0][ray]
        st = f["st"][0][rayindex]
        cp = f["cp"][0][rayindex]
        sp = f["sp"][0][rayindex]
        rlim = f["rlimsurf"][nintersecs][rayindex]
        cartsurf[NinsertedPoints,0] = rlim * st * cp
        cartsurf[NinsertedPoints,1] = rlim * st * sp
        cartsurf[NinsertedPoints,2] = rlim * f["ct"][0][rayindex]
        # Displace the origin from (0,0,0) to the position of atom inuc
        cartsurf[NinsertedPoints,0] += xyz[inuc,0];
        cartsurf[NinsertedPoints,1] += xyz[inuc,1];
        cartsurf[NinsertedPoints,2] += xyz[inuc,2];

        NinsertedPoints += 1


    with open(surf_file + ".csv", "w") as csvsurf:
      logging.info("SURFACE POINTS SAVED TO     {0}.csv".format(surf_file))
      for i in range(nsurfpoints):
        csvsurf.write("{0},{1},{2}\n".format(cartsurf[i,0],cartsurf[i,1],cartsurf[i,2]))

