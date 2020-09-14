#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function

from math import sqrt, sin, cos, radians
import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
import numpy as np

# ----------------------------------------------------------------------
# You need to fill appropriately the following data:
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
asym_ats = [1,5,9]
labels   = ["As"]*4 + ["Li"]*4 + ["Mg"]*4
# ----------------------------------------------------------------------

nats = ats.shape[0]
alpha_rad = radians(alpha_deg)
beta_rad  = radians(beta_deg)
gamma_rad = radians(gamma_deg)
# Metric tensor (BOHR^2)
G = np.array([
              [a_bohr*a_bohr                ,
               a_bohr*b_bohr*cos(gamma_rad),
               a_bohr*c_bohr*cos(beta_rad) ],
              [a_bohr*b_bohr*cos(gamma_rad), 
               b_bohr*b_bohr, 
               b_bohr*c_bohr*cos(alpha_rad)],
              [a_bohr*c_bohr*cos(beta_rad) ,
               b_bohr*c_bohr*cos(alpha_rad), 
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
               -cos(gamma_rad)/(a_bohr*sin(gamma_rad)), 
               b_bohr*c_bohr*(cos(alpha_rad)*cos(gamma_rad)-cos(beta_rad))/(V*sin(gamma_rad))],
              [0, 
               1./(b_bohr*sin(gamma_rad)), 
               a_bohr*c_bohr*(cos(beta_rad)*cos(gamma_rad)-cos(alpha_rad))/(V*sin(gamma_rad))],
              [0., 
               0., 
               a_bohr*b_bohr*sin(gamma_rad)/V]
              ], dtype=np.float64)
# Fractional to cartesian conversion, c = M-1f
Minv = np.array([
                 [a_bohr, 
                  b_bohr*cos(gamma_rad),
                  c_bohr*cos(beta_rad)],
                 [0., 
                  b_bohr*sin(gamma_rad),
                  c_bohr*(cos(gamma_rad)-cos(beta_rad)*cos(gamma_rad))/(sin(gamma_rad))],
                 [0., 
                  0., 
                  V/(a_bohr*b_bohr*sin(gamma_rad))]
                ], dtype=np.float64)

#########################

for iat in range(ats.shape[0]): 
  ats[iat] = np.dot(M, ats[iat])

nimages = 1
ats_supercell = ats

for ximage in range(-nimages,nimages+1):
  for yimage in range(-nimages,nimages+1):
    for zimage in range(-nimages,nimages+1):
      if ximage == yimage == zimage == 0:
        continue
      for i in range(ats.shape[0]):
        translation = np.array([ximage, yimage, zimage],dtype=np.float64)
        ats_supercell = np.concatenate((ats_supercell,[ats[i] + translation]),axis=0)


def distcryst(fi, G):
  return sqrt(np.dot(fi.T,np.dot(G,fi)))

nat_supercell = ats_supercell.shape[0]
pairs_dist = np.zeros( len(asym_ats)*(ats_supercell.shape[0]-1), 
                      dtype= np.dtype("i8, i8, f8, i8") )
ij = 0
for i in asym_ats:
  for j in range(nat_supercell):
    if i == j+1:
      continue
    #(i == j) ? continue : pass
    fi = ats[i-1]
    fj = ats_supercell[j]
    dij = distcryst(fi-fj,G)
    # if unique, then add
    unique_pair = True
    ids_same_dist = np.where( np.isclose(pairs_dist['f2'], dij))
    #print(ids_same_dist[0])
    if ids_same_dist[0].size != 0:
      for pair_id in ids_same_dist[0]:
        pair  = pairs_dist[pair_id]
        lbl1a = labels[pair[0] -1]
        lbl1b = labels[(pair[1]  % nats)-1]
        lbl2a = labels[i]
        lbl2b = labels[j % nats]
        #if (lbl1a == lbl2a and lbl1b == lbl2b) or (lbl1a == lbl2b and lbl1b == lbl2a):
        if (lbl1a == lbl2a and lbl1b == lbl2b):
          pairs_dist[pair_id][3] += 1
          unique_pair = False
        elif (lbl1a == lbl2b and lbl1b == lbl2a):
          unique_pair = False

    if unique_pair:
      pairs_dist[ij] = i,j+1,dij,1
      ij += 1

pairs_dist = np.delete( pairs_dist, np.where(pairs_dist["f2"] < 1e-1)) 
pairs_dist.sort(order=['f0','f2'], axis=0)

print("Labels    i      j    Distance  Mult.")
print("------  -----  -----  --------  -----")
for pair in pairs_dist:
  if pair[2] < 15:
    lbla = labels[pair[0]-1]
    lblb = labels[(pair[1] % nats)-1]
    print("{:2s}--{:2s}  ".format(lbla,lblb),end='')
    print("{:5d}--{:5d} {:8.4f} {:4d}  ".format(*pair))
  
