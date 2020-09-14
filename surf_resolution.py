#!/usr/bin/env python

# Measures the angular resolution in units of surface area (Bohr^2) covered by the grid point.

from __future__ import print_function

import h5py as h5

from select_file_dialog import gui_fname
dgrid_wfn = gui_fname(filter="surf")

f = h5.File(dgrid_wfn, "r")

# Resolution closest point of the surface to the nucleus
print("Angular resolution at argmin |r_surf - R|:", end='')
min_res = f["rlimsurf"][0].min()**2 * f["angw"][0].min()
max_res = f["rlimsurf"][0].min()**2 * f["angw"][0].max()
print("   [ {0:15.7f} {1:15.7f} ] Bohr^2".format(min_res,max_res))
# Resolution farthest point of the surface to the nucleus
print("Angular resolution at argmax |r_surf - R|:", end='')
min_res = f["rlimsurf"][0].max()**2 * f["angw"][0].min()
max_res = f["rlimsurf"][0].max()**2 * f["angw"][0].max()
print("   [ {0:15.7f} {1:15.7f} ] Bohr^2".format(min_res,max_res))
