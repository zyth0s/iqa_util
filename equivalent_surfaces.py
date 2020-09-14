#!/usr/bin/env python
# coding: utf-8

# Code to generate equivalent surfaces by symmetry.
# Limitations: 
#  1. Only one intersection per ray
# REQUIRED: Find out space group, asymmetric unit, Wyckoff label of each atom
# in the asymmetric unit and the generators of the equivalent positions.
# WARNING: The atoms generated are not always the same as in *.fhi or *.fhiout
# but if they are equivalent by translation (only) then you can use them directly.
# Otherwise you need to figure out which is the correct ordering.
# Reference: Symmetry Relationships Between Crystal Structures. IUCr book by Mueller
# TODO: convert mapping into the matrix-column format (W,w).
# TODO: check that it is an isometric mapping and a crystallographic operation.
# TODO: get the geometric interpretation from (W,w)
# DONE: Symmetry operations that mix coordinates, e.g.
#         x -> -y
#         y -> z
#         z -> x
#        are allowed now (.subs() method).
# DONE: Convert cell parameters to metric tensor G and use it.
# DONE: Non-orthogonal crystallographic axes. 
#       Check with https://cci.lbl.gov/cctbx/frac_cart.html

from __future__ import print_function

import os
import sys
import math
import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
import h5py as h5
import numpy as np
from sympy import symbols
from sympy import *
# You need to fill appropriately the following data:
dgrid_wfn = "MgLiAs_SP-TiNiSi-pbe-484.fhi"
asym_ats = [1,5,9]
a_bohr = 12.73621067
b_bohr =  7.57751843
c_bohr = 14.77050051
alpha_deg = 90
beta_deg  = 90
gamma_deg = 90
# Atoms in the asymmetric unit (BOHR):
# taken from *.fhi
ats = np.array([
                [   0.48407993,    1.89437961,   -4.79125776  ], # 1
                [  -0.48473572,   -1.89437961,    4.79132792  ],
                [   5.88364676,   -1.89437961,    2.59352294  ],
                [  -5.88293977,    1.89437961,   -2.59339000  ],
                [   1.70598383,    1.89437961,    0.91416489  ], # 5
                [  -1.70617755,   -1.89437961,   -0.91417360  ],
                [   4.66214468,   -1.89437961,   -6.47110738  ],
                [  -4.66195670,    1.89437961,    6.47108389  ],
                [   3.19436161,    1.89437961,    5.56415374  ], # 9
                [  -3.19411313,   -1.89437961,   -5.56419170  ],
                [   3.17309329,   -1.89437961,   -1.82077585  ],
                [  -3.17338737,    1.89437961,    1.82064291  ]
              ], dtype=np.float64)

x,y,z = symbols('x y z')
# Operations to generate equivalent Wyckoff positions (for each atom in asymmetric unit):
# Pnma (62), 4c Wyckoff position generators
# ops[inuc][nop][component]
# use find_ucell_generators() for help.
ops = {
        1: [ # 1st atom: 4c
          [-x      , -y     , -z],
          [-x+1./2., -y     , z+1./2.],
          [x+1./2.  , -y+1./2., -z + 1./2.]
        ],
        5: [ # 2nd atom: 4c
          [-x      , -y     , -z],
          [-x+1./2., -y     , z+1./2.],
          [x+1./2.  , -y+1./2., -z + 1./2.]
        ],
        9: [ # 3rd atom: 4c
          [-x      , -y     , -z],
          [-x+1./2., -y     , z+1./2.],
          [x+1./2.  , -y+1./2., -z + 1./2.]
        ]
      }
# CAVEAT: The mappings given (in ITA A or BCS) for a certain Wyckoff position
#         are for the position of the atomic nucleus only. It may be located 
#         at a special position, whith some coordinates fixed by symmetry. 
#         Therefore, the mapping is simplified. The surface points will, however, 
#         be located at general positions - in general. So the user has to guess 
#         which mapping for the general positions, that becomes the same mapping
#         listed for the special positon when the position to map is the atomic
#         nucleus, generates the correct surface. This is checked graphically.
#         Then, use the general mapping here.
# ****************************************************************
# General position crystallographic symmetry operations
# Pnma (62), Wyckoff position generators
gops = [ # exclude identity operation
          [-x+1./2., -y      ,  z+1./2.], #1
          [-x      ,  y+1./2., -z      ], #2
          [ x+1./2., -y+1./2., -z+1./2.], #3
          [-x      , -y      , -z      ], #4
          [ x+1./2.,  y      , -z+1./2.], #5
          [ x      , -y+1./2.,  z      ], #6
          [-x+1./2.,  y+1./2.,  z+1./2.], #7
       ]

orbits = { # neq: [Wyckoff equiv atoms]
          1: [ 2, 3, 4],
          5: [ 6, 7, 8],
          9: [10,11,12]
         }
# ****************************************************************

alpha_rad = math.radians(alpha_deg)
beta_rad  = math.radians(beta_deg)
gamma_rad = math.radians(gamma_deg)
# Metric tensor (BOHR^2)
G = np.array([
              [a_bohr*a_bohr                ,
               a_bohr*b_bohr*math.cos(gamma_rad),
               a_bohr*c_bohr*math.cos(beta_rad) ],
              [a_bohr*b_bohr*math.cos(gamma_rad), 
               b_bohr*b_bohr, 
               b_bohr*c_bohr*math.cos(alpha_rad)],
              [a_bohr*c_bohr*math.cos(beta_rad) ,
               b_bohr*c_bohr*math.cos(alpha_rad), 
               c_bohr*c_bohr                ]
             ], dtype=np.float64)
logging.info("Metric tensor G:")
logging.info("  {:7.3f} {:7.3f} {:7.3f}".format(G[0][0],G[0][1],G[0][2]))
logging.info("  {:7.3f} {:7.3f} {:7.3f}".format(G[1][0],G[1][1],G[1][2]))
logging.info("  {:7.3f} {:7.3f} {:7.3f}".format(G[2][0],G[2][1],G[2][2]))
V= sqrt(np.linalg.det(G))
logging.info("Volume = {:15.4f} Bohr^3".format(V))
# Cartesian to fractional conversion, f = Mc
# https://en.wikipedia.org/wiki/Fractional_coordinates
# Derived with Gram-Schmidt?
M = np.array([
              [1./a_bohr, 
               -math.cos(gamma_rad)/(a_bohr*math.sin(gamma_rad)), 
               b_bohr*c_bohr*(math.cos(alpha_rad)*math.cos(gamma_rad)-math.cos(beta_rad))/(V*math.sin(gamma_rad))],
              [0, 
               1./(b_bohr*math.sin(gamma_rad)), 
               a_bohr*c_bohr*(math.cos(beta_rad)*math.cos(gamma_rad)-math.cos(alpha_rad))/(V*math.sin(gamma_rad))],
              [0., 
               0., 
               a_bohr*b_bohr*math.sin(gamma_rad)/V]
              ], dtype=np.float64)
# Fractional to cartesian conversion, c = M-1f
Minv = np.array([
                 [a_bohr,
                  b_bohr*math.cos(gamma_rad),
                  c_bohr*math.cos(beta_rad)],
                 [0.,
                  b_bohr*sin(gamma_rad),
                  c_bohr*(math.cos(gamma_rad)-math.cos(beta_rad)*math.cos(gamma_rad))/(sin(gamma_rad))],
                 [0.,
                  0.,
                  V/(a_bohr*b_bohr*sin(gamma_rad))]
                ], dtype=np.float64)

#########################
def create_Seitz_symbol(op):
  logging.debug("    Mapping: {}, {}, {}".format(*op))
  # Do with origin, x, y and z vectors
  o    = np.array([0.,0.,0.], dtype=np.float64)
  xvec = np.array([1.,0.,0.], dtype=np.float64)
  yvec = np.array([0.,1.,0.], dtype=np.float64)
  zvec = np.array([0.,0.,1.], dtype=np.float64)
  w = np.zeros(3)
  w[0] = op[0].subs([(x,o[0]),(y,o[1]),(z,o[2])]).evalf()
  w[1] = op[1].subs([(x,o[0]),(y,o[1]),(z,o[2])]).evalf()
  w[2] = op[2].subs([(x,o[0]),(y,o[1]),(z,o[2])]).evalf()
  W = np.zeros((3,3))
  W[0][0] = op[0].subs([(x,xvec[0]),(y,xvec[1]),(z,xvec[2])]).evalf() - w[0]
  W[1][0] = op[1].subs([(x,xvec[0]),(y,xvec[1]),(z,xvec[2])]).evalf() - w[1]
  W[2][0] = op[2].subs([(x,xvec[0]),(y,xvec[1]),(z,xvec[2])]).evalf() - w[2]
  #
  W[0][1] = op[0].subs([(x,yvec[0]),(y,yvec[1]),(z,yvec[2])]).evalf() - w[0]
  W[1][1] = op[1].subs([(x,yvec[0]),(y,yvec[1]),(z,yvec[2])]).evalf() - w[1]
  W[2][1] = op[2].subs([(x,yvec[0]),(y,yvec[1]),(z,yvec[2])]).evalf() - w[2]
  #
  W[0][2] = op[0].subs([(x,zvec[0]),(y,zvec[1]),(z,zvec[2])]).evalf() - w[0]
  W[1][2] = op[1].subs([(x,zvec[0]),(y,zvec[1]),(z,zvec[2])]).evalf() - w[1]
  W[2][2] = op[2].subs([(x,zvec[0]),(y,zvec[1]),(z,zvec[2])]).evalf() - w[2]
  logging.debug("    Seitz symbol (W|w):")
  logging.debug("      {:7.3f} {:7.3f} {:7.3f} | {:7.3f} ".format(W[0][0],W[0][1],W[0][2],w[0]))
  logging.debug("      {:7.3f} {:7.3f} {:7.3f} | {:7.3f} ".format(W[1][0],W[1][1],W[1][2],w[1]))
  logging.debug("      {:7.3f} {:7.3f} {:7.3f} | {:7.3f} ".format(W[2][0],W[2][1],W[2][2],w[2]))
  return W, w

def charazterize_cryst_symm_op(W,w):

  if not np.allclose(W.T * G * W, G):
    logging.warning("Not an isommetry")
  detW = np.linalg.det(W)
  traceW = np.trace(W)
  if np.isclose(detW,1.0): 
    logging.info("    det(W) = {:+7.3f}: Rotation".format(detW))
    if traceW != 3:
      _, u = np.linalg.eig(W)
      for i in range(3):
        if _[i] == 1:
          u = u[:,i]
    if traceW == 3:
      type = 1
    if traceW == 2:
      type = 6
    if traceW == 1:
      type = 4
    if traceW == 0:
      type = 3
    if traceW == -1:
      type = 2
    k = type
  elif np.isclose(detW,-1.0):
    logging.info("    det(W) = {:+7.3f}: Rotoinversion".format(detW))
    if traceW != -3:
      _, u = np.linalg.eig(W)
      for i in range(3):
        if _[i] == -1:
          u = u[:,i]
    if traceW == -3:
      type = -1
      k = 2
    if traceW == -2:
      type = -6
      k = 6
    if np.trace(W) == -1:
      type = -4
      k = 4
    if traceW == 0:
      type = -3
      k = 6
    if traceW == 1:
      type = -2
      k = 2
  else:
    logging.warning("    det(W) = {:+7.3f}".format(detW))

  if not np.allclose(np.linalg.matrix_power(W,k),np.eye(3)):
    logging.warning("Order is wrong!")

  _ = np.eye(3,dtype=np.float64)
  for i in range(1,k):
    _ += np.linalg.matrix_power(W,i)
  screw_glide_component = np.zeros(3,dtype=np.float64)
  screw_glide_component = np.dot(_,w) / k
  #logging.info("    Screw/glide component:  {:7.3f} {:7.3f} {:7.3f}".format(*screw_glide_component))


  logging.info("    Crystallographic symmetry operation: {:2d} (order {:2d})".format(type,k))
  if type == 1:
    if np.allclose(w,np.zeros(3)):
      logging.info("    Identity")
    else:
      logging.info("    Translation vector {:7.3f} {:7.3f} {:7.3f}".format(*w))
  if type == 6:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
    if not np.allclose(screw_glide_component,np.zeros(3)):
      logging.info("    Screw component:  {:7.3f} {:7.3f} {:7.3f}".format(*screw_glide_component))
  if type == 4:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
    if not np.allclose(screw_glide_component,np.zeros(3)):
      logging.info("    Screw component:  {:7.3f} {:7.3f} {:7.3f}".format(*screw_glide_component))
  if type == 3:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
    if not np.allclose(screw_glide_component,np.zeros(3)):
      logging.info("    Screw component:  {:7.3f} {:7.3f} {:7.3f}".format(*screw_glide_component))
  if type == 2:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
  if type == -1:
    logging.info("    Point of inversion:  {:7.3f} {:7.3f} {:7.3f}".format(0.5*w[0],0.5*w[1],0.5*w[2]))
  if type == -6:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
  if type == -4:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
  if type == -3:
    logging.info("    Fixed axis: {:7.3f} {:7.3f} {:7.3f}".format(*u))
  if type == -2:
    logging.info("    Reflection plane perpendicular to {:7.3f} {:7.3f} {:7.3f}".format(*u))
    if not np.allclose(screw_glide_component,np.zeros(3)):
      logging.info("    Glide component:  {:7.3f} {:7.3f} {:7.3f}".format(*screw_glide_component))

def gen_image_pos(originalpos,W,w, normalize=False):

  # Be careful not to modify originalpos
  startpoint = np.dot(M,originalpos)
  # Restrict to be inside unit cell
  if normalize:
    startpoint[0] -= math.floor(startpoint[0]) 
    startpoint[1] -= math.floor(startpoint[1])
    startpoint[2] -= math.floor(startpoint[2])
  # Input are fractionary coordinates
  imagepoint = np.zeros(3)
  imagepoint = np.dot(W,startpoint) + w
  imagepoint = np.dot(Minv,imagepoint)
  return imagepoint

def cart2polar(point):

  epsrad = 1e-7
  cost = 0.0
  sint = 0.0
  cosp = 0.0
  sinp = 0.0
  rxy = point[0]*point[0] + point[1]*point[1]
  if rxy < epsrad:
    if point[2] >= 0.0:
      cost = +1.0
    else:
      cost = -1
  else:
    rxy = math.sqrt(rxy)
    cost = point[2] / np.linalg.norm(point)
    sint = math.sqrt( (1.0 - cost) * (1.0 + cost) )
    cosp = point[0] / rxy
    sinp = point[1] / rxy

  return cost, sint, cosp, sinp
#########################

#from select_file_dialog import gui_fname
#from get_coords import get_coords

#surf_file = gui_fname(filter="surf")
#if surf_file == "":
#  print("Canceled. No file selected")
#  exit()

def read_surface(inuc):
  surf_file = "{0}.surf{1:04d}.h5".format(dgrid_wfn,inuc)

  if not os.path.isfile(surf_file):
    logging.warning("There is no surface file {0:s}".format(surf_file))
    sys.exit()
  if not os.access(surf_file, os.R_OK):
    logging.warning("File {0:s} is not readable".format(surf_file))
    sys.exit()

  f = h5.File(surf_file)
  return f,surf_file


def apply_mapping_to_surface(inuc, iop):

  W, w = create_Seitz_symbol(ops[inuc][iop])
  charazterize_cryst_symm_op(W,w)

  # Position of generated atom
  atnew = gen_image_pos(ats[inuc-1],W,w)
  foundLabel = False
  for i,at in enumerate(ats):
    if np.allclose(atnew, at,atol=1e-7):
      foundLabel = True
      logging.info("    Image atom is {0:4d}".format(i+1))
  if not foundLabel:
    logging.info("    Image atom is ????: {0:7.3f} {1:7.3f} {2:7.3f}".format(*atnew))

  surf_file_name = surf_file + "_sym{0:04d}".format(iop+1) 
  logging.info("    Creating surface in {} ...".format(surf_file_name))
  fnew = h5.File(surf_file_name, "w")
  ctnew       = fnew.create_dataset("ct"       , (1,npang , ) , dtype = 'f8')
  stnew       = fnew.create_dataset("st"       , (1,npang , ) , dtype = 'f8')
  cpnew       = fnew.create_dataset("cp"       , (1,npang , ) , dtype = 'f8')
  spnew       = fnew.create_dataset("sp"       , (1,npang , ) , dtype = 'f8')
  angwnew     = fnew.create_dataset("angw"     , (1,npang , ) , dtype = 'f8')
  nlimsurfnew = fnew.create_dataset("nlimsurf" , (1,npang , ) , dtype = 'i8')
  rlimsurfnew = fnew.create_dataset("rlimsurf" , (10,npang , ) , dtype = 'f8')

  # Insert the points corresponding to the first intersection in each ray
  for i in range(npang):
    st = f["st"][0][i]
    cp = f["cp"][0][i]
    sp = f["sp"][0][i]
    rlim = f["rlimsurf"][0][i]
    # Spherical to cartesian
    point = np.zeros(3)
    point[0] = rlim * st * cp
    point[1] = rlim * st * sp
    point[2] = rlim * f["ct"][0][i]
    # Change from relative origin (old nucleus) to global origin
    point += ats[inuc-1]
    # Apply sym op to point
    pointnew = gen_image_pos(point,W,w)
    #print("{0:7.3f}, {1:7.3f}, {2:7.3f} ".format(*pointnew))
    # Change from global origin to relative origin (new nucleus)
    pointnew -= atnew
    # Convert to polar
    rlimp = np.linalg.norm(pointnew)
    cost, sint, cosp, sinp = cart2polar(pointnew)
    # Save
    ctnew[0,i] = cost
    stnew[0,i] = sint
    cpnew[0,i] = cosp
    spnew[0,i] = sinp
    angwnew[0,i] = f["angw"][0][i]
    nlimsurfnew[0,i] = f["nlimsurf"][0][i]
    rlimsurfnew[0,i] = rlimp

# Helper to find out the gop that converts asym atom into other equiv atom
def find_ucell_generators(atol=1e-2):
  for i in orbits: # asym_ats
    for j in orbits[i]: # equiv_ats
      logging.info("    Atom {0:4d} generated from atom {1:4d} with op ...".format(j,i))
      for iop,op in enumerate(gops):
        W, w = create_Seitz_symbol(op)
        # Position of generated atom
        atnew = gen_image_pos(ats[i-1],W,w)
        translation = [abs(atnew[0]-ats[j-1][0]),
                       abs(atnew[1]-ats[j-1][1]),
                       abs(atnew[2]-ats[j-1][2])]
        transmod = np.zeros(3)
        if not np.allclose(translation, np.zeros(3),atol=atol):
          transmod = [max(translation[0],a_bohr) % min(translation[0],a_bohr),
                      max(translation[1],b_bohr) % min(translation[1],b_bohr),
                      max(translation[2],c_bohr) % min(translation[2],c_bohr)]
        if np.allclose(transmod, np.zeros(3),atol=1e-2):
          logging.info("                {0:12s}, {1:12s}, {2:12s}".format(*op))
          #charazterize_cryst_symm_op(W,w)

# -------------------------------------------------------------

find_ucell_generators()

for inuc in asym_ats:
  nops = len(ops[inuc])
  f, surf_file = read_surface(inuc)
  npang = len(f["ct"][0])
  for iop in range(nops):
    logging.info("* Atom {0:4d}: generating equivalent surface with op {1:4d}".format(inuc,iop+1))
    apply_mapping_to_surface(inuc, iop)


#equivats = set([v+1 for v in range(len(ats[:,1]))]).difference(set(asym_ats))
#for j in equivats:
#  #print(j)
#  logging.warn("    Atom {0:4d} generated from ...".format(j))
#  for i in asym_ats:
