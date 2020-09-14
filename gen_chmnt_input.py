#!/usr/bin/env python

# Generate ChemInt input from a template

dgrid_sufix="gms"
name = "ch4"
nthreads = 6

chemint_input = """\
{0:s}.{1:s}
#-------------------------------
threads {2:4d}
loglevel debug # trace|debug|info|error|off

ode bs3 # bs3|ark3|sa4|ark4|z4|dp5|ck5|f5|ark5|verner6|f8      
#epsilon 1e-4 # better not to use # bisection epsilon = 1e-4
epsiscp 0.20 # PERFORMANCE: as big as possible
#surf 5 15 10 20 12 22 8 18 7 17
#surf 4 6 9 11 14 16 19 21
isosurface off
rmaxsurf 8.
##rmaxatom  0 10.
##rmaxatomz 0 10.
lebedev 6000  # PERFORMANCE: as small as possible
##lebedev_atom   1 6000
##lebedev_atom_z 1 6000
##rsearch 0 11 0.40 0.55 0.76 1.05 1.44 1.99 2.74 3.78 5.21 7.18 9.90

#sphere

epsdet 1e-12

#fixrmaxsurf

#scdm

# DFT
xcdft pbe 
# Choose the functional among:
# lda  | svwn  | b3lyp   | camb3lyp | bhandhlyp | pbe0 | pw86pbe | blyp   | b971
# b97  | tpssh | revtpss | tpssm    | m05       | m06  | m06-2x  | m06-hf | m06-l
# olyp | rpbe  | m11-l   | b3p86    | pw91      | bp86
# or select the X and C functionals using LibXC constants:
# ------------------------------
#    X   C   XC   Name
#  --- ---  ---   -----
#    1   7        SVWN5
#           402   B3LYP
#  106 131        BLYP
#  101 130        PBE
# ------------------------------
# Examples:
#libxc 0 101 130 # X and C separated
#libxc 0   0 402 # XC together

#nomonadic
#romberg
#multipolar 2 1.2 3 5.0
#epsbice 1e-8

# RADIAL QUADRATURE
# ------------------------------
rquad 200 1 0 # 3rd defined: outside beta
rquad 200 1 2 # 4th defined:  inside beta
mpr 2
alphaexp 0.6
#rbragg 1 0.02

# OUTSIDE BETA-SPHERE
# ------------------------------
lmax 6
lmax_atom   1 6
#lmax_atom_z 1 6
rint  0 3 # use 3rd
#rintz 0 3 # use 3rd

# INSIDE BETA-SPHERE
# ------------------------------
betasphere # / nobetasphere
betafac 1 0.6
betafac 2 0.6
betafac 3 0.6
#betafacz 1 0.6
lebedevbeta 434
#lebedevbeta_atom   1 434
#lebedevbeta_atom_z 1 434
lmaxbeta 4
rintb  0 4 # use 4th
#rintbz 0 4 # use 4th

# SELF ENERGIES
# ------------------------------
#self 1 2 3 4
#self 1 2
# INTERACTIONS
# ------------------------------
allint
#int 1 2
"""

with open(name + ".chmntin", "w") as f:
  f.write(chemint_input.format(name,dgrid_sufix,nthreads))

