#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: check it is version agnostic Python 2.x & 3.x

from __future__ import print_function, division, unicode_literals

from sys import exit
from sys import argv
import json
import copy
import textwrap
from operator import itemgetter

try:
  from matplotlib import rcParams
  rcParams["font.family"] = "serif"
  rcParams["mathtext.fontset"] = "stix"
  from matplotlib import pyplot as plt
except:
  print("Matplotlib module is not present!")
import numpy as np

import logging
logging.basicConfig()
#logging.getLogger().setLevel(logging.DEBUG)


from mendeleev import element

from select_file_dialog import gui_fname

if len(argv) == 1:
  datafile = gui_fname()
else:
  datafile = argv[1]

if not datafile:
  #raise Exception("Canceled. No file selected")
  exit()
print("Data source: " + datafile.decode())
if datafile.endswith(b"json"):
    data = json.load( open( datafile ) )
elif datafile.endswith(b"iqapack"):
    try:
        import msgpack
        data = msgpack.load( open( datafile, "rb" ), encoding="utf-8")
    except ImportError as e:
        raise ImportError(e)
elif datafile.endswith(b"json.gz"):
    try:
        import gzip
        data = json.load( gzip(open( datafile )) )
    except ImportError as e:
        raise ImportError(e)
#data = msgpack.load( open( datafile, "rb" ) )

#datafile = input("Enter the name of the JSON file containing IQA_Molecule (inside quotes) data: \n")
#try:
#  datafile = "Fe2CO8_anion2.adf.json" # FIXME: hardcoded for easing development
#  data = json.load( open( datafile ) )
#except:
#  print("Wrong! Example: 'ch4.adf.json'")
#  datafile = input("Enter the name of the JSON file containing IQA_Molecule (inside quotes) data: \n")
#  data = json.load( open( datafile ) )

#molec = type('IQA_Molecule', (object,), data["molecule_or_group"])

class IQA_Mono:
  def __init__(self,dictionary):
    for k, v in dictionary.items():
      if k == "E_SELFEL_NUCLEI":
        setattr(self, k, np.array(v) )
      else:
        setattr(self, k, v )
    if hasattr(self,"Z"):
      #setattr(self, "label", element(int(self.Z)).symbol + "_" + str(self.inuc))
      #setattr(self, "label", element(int(self.Z)).symbol + "(" + str(self.inuc) + ")")
      setattr(self, "label", "{:_<2s}{:<4d}".format(element(int(self.Z)).symbol,int(self.inuc)))
      #setattr(self, "label", "{:_<2s}{:<4d}".format(element(int(self.Z)).symbol,self.inuc))


class IQA_Pair:
  def __init__(self,dictionary):
    for k, v in dictionary.items():
        if k == "COULOMB_l" or k == "EXCHANGE_l" or k == "EXCHANGE_CORRELATION_l":
            setattr(self, k, np.array(v) )
        else:
            setattr(self, k, v )


class IQA_Molecule:
  def __init__(self,dictionary):
    setattr(self,"coarseness",[])
    try:
        rootname = ""
        if "molecule_or_group" in dictionary: 
          rootname = "molecule_or_group"
        elif "molecule" in dictionary:          
          rootname = "molecule"
        elif "fragment" in dictionary:          
          rootname = "fragment"
        else:
          raise Exception("File is empty")
        for k, v in dictionary[rootname].items():
          if k == "monos":
            setattr(self,k,[])
            for mono in v:
              self.monos.append( IQA_Mono( mono ) )
              self.coarseness.append( [mono["inuc"]] ) 
          elif k == "pairs":
            setattr(self,k,{})
            for pair in v:
              #self.pairs.append( IQA_Pair( pair ) )
              self.pairs[(pair["ia"],pair["ib"])] =  IQA_Pair( pair ) 
          else:
            setattr(self, k, v )
    except ValueError as error:
        print("Error in the data file.")
        print(repr(error))
    except Exception as error:
        print("Caught exception " + repr(error))
        print("Building molecule from scratch")
        self.monos = []
        self.pairs = {}
        self.E_KINETIC = 0.0
        self.E_POTENTIAL = 0.0
        self.E_TOTAL = 0.0
        self.E_NUC_NUC = 0.0
        self.E_NUC_EL = 0.0
        self.E_EL_EL = 0.0
        self.E_EL_EL_COULOMB = 0.0
        self.E_EL_EL_EXCHANGE_CORRELATION = 0.0
        self.E_EL_EL_EXCHANGE = 0.0
        self.E_EL_EL_CORRELATION = 0.0
        self.TWO_KINETIC_PLUS_POTENTIAL = 0.0
        self.VIRIAL_RATIO = 0.0
        self.E_NET = 0.0
        self.E_INTERACTION = 0.0
        self.E_INTERACTION_CLASSICAL = 0.0
        self.E_INTERACTION_EXCHANGE_CORRELATION = 0.0
        self.CHARGE = 0.0


  def plot_intra_convergence(self, imono):
    mono = self.monos[imono-1]
    fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1,figsize=(17,6))
    plt.suptitle( "Fragment {}".format(mono.label), fontsize=25 )
    ax1.set_yscale("log")
    #ax1.set_title( "Fragment {}".format(mono.label), fontsize=20 )
    ax1.set_ylabel(r"Intra Coulomb (Ha)", fontsize=20)
    ax1.set_xlabel(r"$l$", fontsize=20)
    ax1.tick_params(labelsize=15,width=5)
    ax1.tick_params(which='both', width=2)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=4)
    ax1.grid(True)
    ax1.plot(mono.COULOMB_l, "r-o")

    ax2.set_yscale("log")
    #ax2.set_title( "Fragment {}".format(mono.label), fontsize=20 )
    ax2.set_ylabel(r"-Intra Exchange (Ha)", fontsize=20)
    ax2.set_xlabel(r"$l$", fontsize=20)
    ax2.tick_params(labelsize=15,width=5)
    ax2.tick_params(which='both', width=2)
    ax2.tick_params(which='major', length=7)
    ax2.tick_params(which='minor', length=4)
    ax2.grid(True)
    ax2.plot([ -x for x in mono.EXCHANGE_l], "r-o")
    plt.show()

  def plot_inter_convergence(self, ia, ib):
    pair = self.pairs[(ia,ib)]
    fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1,figsize=(17,6))
    plt.suptitle( u"Interaction {0}⟷  {1}".format(self.monos[ia-1].label, self.monos[ib-1].label), fontsize=25 )
    ax1.set_yscale("log")
    ax1.set_ylabel(r"Inter Coulomb (Ha)", fontsize=20)
    ax1.set_xlabel(r"$l$", fontsize=20)
    ax1.tick_params(labelsize=15,width=5)
    ax1.tick_params(which='both', width=2)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=4)
    ax1.grid(True)
    ax1.plot(pair.COULOMB_l, "r-o")

    ax2.set_yscale("log")
    ax2.set_ylabel(r"-Inter Exchange (Ha)", fontsize=20)
    ax2.set_xlabel(r"$l$", fontsize=20)
    ax2.tick_params(labelsize=15,width=5)
    ax2.tick_params(which='both', width=2)
    ax2.tick_params(which='major', length=7)
    ax2.tick_params(which='minor', length=4)
    ax2.grid(True)
    ax2.plot([ -x for x in pair.EXCHANGE_l], "r-o")
    plt.show()

  def print_intra( self, imono, format="txt"):
    mono = self.monos[imono-1]
    print("   +++ Intra components of fragment: {0:s}".format(mono.label))
    print("   " + "-"*49)
    print("   +++ kinetic energy             = {0:+16.8f}".format(mono.E_KINETIC))
    if hasattr(mono, "E_SELF_NUC_NUC"):
        print("   +++ self nuc-nuc energy        = {0:+16.8f}".format(mono.E_SELF_NUC_NUC))
    print("   +++ self nuc-el energy         = {0:+16.8f}".format(mono.E_NUC_EL))
    for n in range(len(mono.E_SELFEL_NUCLEI)):
      if ( n < int(mono.inuc)-1 ): 
        print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
                 mono.inuc, n+1, mono.E_SELFEL_NUCLEI[n]))
      else:
        print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
                 int(mono.inuc), n+2, mono.E_SELFEL_NUCLEI[n]))
    print("   +++ coulomb energy             = {0:+16.8f}".format(mono.E_EL_EL_COULOMB))
    # TODO: if  not noxc
    print("   +++ xc      energy             = {0:+16.8f}".format(mono.E_EL_EL_EXCHANGE_CORRELATION))
    print("   +++ x       energy             = {0:+16.8f}".format(mono.E_EL_EL_EXCHANGE))
    print("   +++  c      energy             = {0:+16.8f}".format(mono.E_EL_EL_CORRELATION))
    print("   +++ Total ee rep               = {0:+16.8f}".format(mono.E_EL_EL))
    print("   +++ Localization index (xc)    = {0:+16.8f}".format(mono.LOCALIZATION_INDEX))
    print("   +++ Localization index (x )    = {0:+16.8f}".format(mono.LOCALIZATION_INDEX_EXCHANGE))
    print("   +++ Localization index ( c)    = {0:+16.8f}".format(mono.LOCALIZATION_INDEX_CORRELATION))
    print("   +++ Self-energy                = {0:+16.8f}".format(mono.E_SELF))
    #
    print("   +++ Electrons & charge         = {0:+16.8f} & {1:+16.8f}".format(
                                               mono.N_ELECTRONS,
                                               mono.CHARGE))

  def print_inter( self, ia, ib, format="txt"):
    pair = self.pairs[(ia,ib)]
    assert pair.ia == ia
    assert pair.ib == ib
    #title = "   +++ Interaction of fragments: {0:s}⟷  {1:s}".format(self.monos[ia-1].label, self.monos[ib-1].label)
    print("   +++ Interaction of ", end='')
    title = self.monos[ia-1].label
    title = textwrap.fill(title, width=80,initial_indent='', subsequent_indent=' '*22)
    print(title)
    #print("⟷ "   )
    #print(" "*22 + "with")
    title = " "*17 + "with " + self.monos[ib-1].label
    title = textwrap.fill(title, width=80,initial_indent='', subsequent_indent=' '*22)
    print(title)
    #title = self.monos[ib-1].label
    #title = textwrap.fill(title, width=80,initial_indent='', subsequent_indent=' '*33)
    #print(title)
    print(" " + "-"*89)
    print("                                                   exact       multipolar             diff")
    print(" " + "-"*89)
    if hasattr(pair, "dist"):
        print("   +++ Distance                         {0:+16.8f}".format(pair.dist))
    #print("   +++ Bicentric analysis of pair {0} {1}".format(self.coarseness[ia-1], self.coarseness[ib-1]))
    print("   +++ NN total repulsion               {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_NUC_NUC,
                                                   pair.E_NUC_NUC,
                                                   0.0))
    print("   +++ energy elec({0:3d})-nuc({1:3d})        {2:+16.8f} {3:+16.8f} {4:+16.8f}".format(
                                                         ia,
                                                         ib,
                                                         pair.E_EL_NUC,
                                                         pair.E_EL_NUC_mp,
                                                         pair.E_EL_NUC-pair.E_EL_NUC_mp))
    print("   +++ energy elec({0:3d})-nuc({1:3d})        {2:+16.8f} {3:+16.8f} {4:+16.8f}".format(
                                                         ib,
                                                         ia,
                                                         pair.E_NUC_EL,
                                                         pair.E_NUC_EL_mp,
                                                         pair.E_NUC_EL-pair.E_NUC_EL_mp))
    print("   +++ EE coulomb energy                {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_EL_EL_COULOMB,
                                                   pair.E_EL_EL_COULOMB_mp,
                                                   pair.E_EL_EL_COULOMB-pair.E_EL_EL_COULOMB_mp))
    print("   +++ Clasical interaction (xc out)    {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_INTERACTION_CLASSICAL,
                                                   pair.E_INTERACTION_CLASSICAL_mp,
                                                   pair.E_INTERACTION_CLASSICAL-pair.E_INTERACTION_CLASSICAL_mp))
    print("   +++ EE x     energy                  {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_EL_EL_EXCHANGE,
                                                   pair.E_EL_EL_EXCHANGE_mp,
                                                   pair.E_EL_EL_EXCHANGE-pair.E_EL_EL_EXCHANGE_mp))
    print("   +++ EE  c    energy                  {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_EL_EL_CORRELATION,
                                                   pair.E_EL_EL_CORRELATION_mp,
                                                   pair.E_EL_EL_CORRELATION-pair.E_EL_EL_CORRELATION_mp))
    print("   +++ EE xc    energy                  {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_EL_EL_EXCHANGE_CORRELATION,
                                                   pair.E_EL_EL_EXCHANGE_CORRELATION_mp,
                                                   pair.E_EL_EL_EXCHANGE_CORRELATION-pair.E_EL_EL_EXCHANGE_CORRELATION_mp))
    print("   +++ EE total repulsion               {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_EL_EL,
                                                   pair.E_EL_EL_mp,
                                                   pair.E_EL_EL-pair.E_EL_EL_mp))
    print("   +++ Total interaction                {0:+16.8f} {1:+16.8f} {2:+16.8f}".format(
                                                   pair.E_INTERACTION,
                                                   pair.E_INTERACTION_mp,
                                                   pair.E_INTERACTION-pair.E_INTERACTION_mp))
    print("   +++ Delocalization index (xc)        {0:+16.8f}".format(
                                                   pair.DELOCALIZATION_INDEX))
    print("   +++ Delocalization index (x )        {0:+16.8f}".format(
                                                   pair.DELOCALIZATION_INDEX_EXCHANGE))
    print("   +++ Delocalization index ( c)        {0:+16.8f}".format(
                                                   pair.DELOCALIZATION_INDEX_CORRELATION))

  def print_intra_table(self):
    with open("results.dat","a") as f:
      print("Basin  Label      E_kinetic         E_nuc-el       E_Coulomb_ee     E_Coulomb_total",file=f)
      print("----   ------   --------------   --------------   --------------    --------------",file=f)
      for mono in self.monos:
          print("{0:4d}   {1:6s}   {2:14.7f}   {3:14.7f}   {4:14.7f}   {5:14.7f}".format(
            int(mono.inuc),str(mono.label),mono.E_KINETIC,mono.E_NUC_EL,mono.E_EL_EL_COULOMB,mono.E_SELF),file=f)

  def print_inter_table(self):
    with open("results.dat","a") as f:
      print("Interaction     Label        Dist.         E_nn            E_ne^AB          E_ne^BA        E_Coulomb_ee    E_Coulomb_total",file=f)
      print("---------   -------------   -------   --------------   --------------   --------------    --------------    --------------",file=f)
      pairs_idx = [x for _,x in sorted(zip([pair.dist for pair in self.pairs.values()],self.pairs.keys()))]
      for pair_idx in pairs_idx:
        pair = self.pairs[pair_idx]
        print("{0:4d}-{1:4d}   {2:6s}-{3:6s}   {4:7.4f}   {5:14.7f}   {6:14.7f}   {7:14.7f}   {8:14.7f}   {9:14.7f}".format(
              pair.ia,
              pair.ib,
              self.monos[pair.ia-1].label,
              self.monos[pair.ib-1].label,
              pair.dist,
              pair.E_NUC_NUC,
              pair.E_NUC_EL,
              pair.E_EL_NUC,
              pair.E_EL_EL_COULOMB,
              pair.E_INTERACTION),file=f)

  def print_inter_table_by_coordination(self):
    with open("results.dat","a") as f:
      print("Interaction     Label        Dist.         E_nn            E_ne^AB          E_ne^BA        E_Coulomb_ee    E_Coulomb_total",file=f)
      print("---------   -------------   -------   --------------   --------------   --------------    --------------    --------------",file=f)
      for i in range(1,len(self.monos)+1):
        for j in range(1,len(self.monos)+1):
          if i == j:
            continue
          if i > j:
            #i,j = j,i
            pair = self.pairs[(j,i)]
            #print("{0:4d}-{1:4d}".format(i, j), file=f)
            #continue
          else:
            pair = self.pairs[(i,j)]
          print("{0:4d}-{1:4d}   {2:6s}-{3:6s}   {4:7.4f}   {5:14.7f}   {6:14.7f}   {7:14.7f}   {8:14.7f}   {9:14.7f}".format(
                i,
                j,
                self.monos[i-1].label,
                self.monos[j-1].label,
                pair.dist,
                pair.E_NUC_NUC,
                pair.E_NUC_EL,
                pair.E_EL_NUC,
                pair.E_EL_EL_COULOMB,
                pair.E_INTERACTION),file=f)
        print(" ",file=f)

  def plot_DI_table(self):
    natoms = len(self.monos)
    npairs = len(self.pairs)
    DI_table = np.zeros((natoms,natoms))
    INT_table = np.zeros((npairs,npairs))

    for i in range(natoms):
      DI_table[i,i] = self.monos[i].LOCALIZATION_INDEX

    for pair in self.pairs:
      i = pair.ia - 1
      j = pair.ib - 1
      DI_table[i,j] = pair.DELOCALIZATION_INDEX
      INT_table[i,j] = pair.E_INTERACTION

    #plt.imshow(DI_table, interpolation='nearest', origin='lower')
    #plt.colorbar()
    plt.matshow(DI_table)
    plt.show()

    plt.matshow(INT_table)
    plt.show()

  def plot_eself_pie(self):

    sizes = []
    labels = []
    for mono in self.monos:
      sizes.append(mono.E_SELF/self.E_NET)
      labels.append(mono.inuc)

    plt.pie(sizes, labels=labels, autopct='%1.0f%%', 
            pctdistance=1.1, labeldistance=1.2)
    plt.show()

  def check_consistency(self):
      for i in range(len(self.monos)):
          assert self.monos[i].inuc == i+1

  def demo(self):

    self.plot_intra_convergence(1)
    self.plot_inter_convergence(1,2)
    self.plot_DI_table()
    self.plot_eself_pie()
    self.print_mono(1)
    self.print_pair(1,2)

    lmax = 0
    for qlm in m.monos[0].ELEC_POP_MULTIPOLES_Qlm:
      l = qlm["tuple_element0"]
      lmax = max(l, lmax)

    ind = np.arange(lmax+1)
    width = 0.25
    fig, ax = plt.subplots()

    grouped_by_eme = [0.0 for x in range(lmax+1)]
    for qlm in m.monos[1].ELEC_POP_MULTIPOLES_Qlm:
      ele = qlm["tuple_element0"]
      eme = qlm["tuple_element1"]
      if eme == -1:
        grouped_by_eme[ele] = qlm["tuple_element2"]

    pm1 = ax.bar(ind, grouped_by_eme, width, color="r", bottom=0.0)

    for qlm in m.monos[1].ELEC_POP_MULTIPOLES_Qlm:
      ele = qlm["tuple_element0"]
      eme = qlm["tuple_element1"]
      if eme == 0:
        grouped_by_eme[ele] = qlm["tuple_element2"]

    p0 = ax.bar(ind+width, grouped_by_eme, width, color="tab:gray", bottom=0.0)

    for qlm in m.monos[1].ELEC_POP_MULTIPOLES_Qlm:
      ele = qlm["tuple_element0"]
      eme = qlm["tuple_element1"]
      if eme == 1:
        grouped_by_eme[ele] = qlm["tuple_element2"]

    p1 = ax.bar(ind+2*width, grouped_by_eme, width, color="b", bottom=0.0)


    plt.axhline(linewidth=2, color="k")

    ax.set_title(r"$Q_{lm}$")
    ax.set_xlabel(r"$l$")
    ax.set_xticks(ind + width )
    ax.set_xticklabels( ind )
    ax.autoscale_view()
    plt.show()

#  def _coarsening(self, groupsidx, group_other=False):
#    """
#    Performs a coarse graining of basins
#    E.g. (with groupsidx = [[1,2]], group_other = False )
#    1 and 2 --> 1
#    3       --> 2
#    4       --> 3
#    E.g. (with groupsidx = [[1,2]], group_other = True )
#    1 and 2 --> 1
#    3 and 4 --> 2
#    Besides, reordering of fragments can be done
#    E.g. (with groupsidx = [[2],[1]], group_other = False )
#    1       --> 2
#    2       --> 1
#    3       --> 3
#    4       --> 4
#    """
#    # Notice that there is no empty constructor for the class IQA_Molecule
#    # It would be simpler if there were an empty constructor but that would imply
#    # redundancy with the C++ code, which means syncronizaton, ...
#    # Naive solution:
#    # Merge group data in first element mono/pair of group
#    # Move interactions intragroup to mono terms
#    coarsed_mol = copy.deepcopy(self)
#
#    # Sort elements in each group (repeated elements are kept)
#    for groupidx in groupsidx:
#        groupidx.sort()
#    # Fix groups to include all basins (exhaustive)
#    # and change indexing from [1,] to [0,]
#    basins = [item for sublist in groupsidx for item in sublist]
#    # check that it is sorted in strict increasing order (no repeated elements)
#    basins_set = set(basins)
#    for basin in basins_set:
#        if basins.count(basin) > 1:
#            raise Exception("{:4d} is duplicated".format(basin))
#
#    #assert all(basins[i] < basins[i+1] for i in range(len(basins)-1)) 
#    all_basins = set(range(1,len(self.monos)+1))
#    other_basins = list(all_basins.difference(basins_set))
#    assert all(other_basins[i] < other_basins[i+1] for i in range(len(other_basins)-1)) 
#    if group_other:
#        groupsidx.append( other_basins )
#    else:
#        for basin in other_basins:
#            groupsidx.append( [basin] )
#    coarsed_mol.coarseness = copy.deepcopy(groupsidx)
#
#    for n in range(len(groupsidx)):
#        groupsidx[n] = [ x - 1 for x in groupsidx[n]]
#
#    other_groupsidx = copy.deepcopy(groupsidx)
#    for igroup, groupidx in enumerate(groupsidx):
#        #ndel += len(groupidx[1:])
#        del other_groupsidx[0]
#        #print(groupidx, other_groupsidx)
#        # merge in first element of group
#        group = coarsed_mol.monos[groupidx[0]]
#        if not hasattr(group, "E_SELF_NUC_NUC"):
#            setattr(group, "E_SELF_NUC_NUC", 0.0)
#        for imono in groupidx:
#            # Sum intra terms to the group
#            if imono != groupidx[0]:
#                groupelm = coarsed_mol.monos[imono]
#                group.label                        += ", " + groupelm.label
#                group.E_KINETIC                    += groupelm.E_KINETIC
#                group.Z                            += groupelm.Z
#                group.CHARGE                       += groupelm.CHARGE
#                group.N_ELECTRONS                  += groupelm.N_ELECTRONS
#
#                group.E_NUC_EL                     += groupelm.E_NUC_EL
#                group.E_EL_EL_COULOMB              += groupelm.E_EL_EL_COULOMB
#                group.E_EL_EL_EXCHANGE_CORRELATION += groupelm.E_EL_EL_EXCHANGE_CORRELATION
#                group.E_EL_EL_EXCHANGE             += groupelm.E_EL_EL_EXCHANGE
#                group.E_EL_EL_CORRELATION          += groupelm.E_EL_EL_CORRELATION
#                group.E_EL_EL                      += groupelm.E_EL_EL
#                group.E_SELF                       += groupelm.E_SELF
#
#            # Sum interactions within the group
#            for imono2 in groupidx:
#                if imono >= imono2: # Do not count interactions twice
#                    continue        # elements within the group must be sorted
#                pair = self.pairs[(imono+1,imono2+1)]
#                assert pair.ia == imono+1
#                assert pair.ib == imono2+1
#                group.E_EL_EL_COULOMB                += pair.E_EL_EL_COULOMB
#                group.E_EL_EL_EXCHANGE_CORRELATION   += pair.E_EL_EL_EXCHANGE_CORRELATION
#                group.E_EL_EL_EXCHANGE               += pair.E_EL_EL_EXCHANGE
#                group.E_EL_EL_CORRELATION            += pair.E_EL_EL_CORRELATION
#                group.E_EL_EL                        += pair.E_EL_EL
#                group.E_SELF_NUC_NUC                 += pair.E_NUC_NUC
#                group.E_NUC_EL                       += pair.E_NUC_EL + pair.E_EL_NUC
#                group.E_SELF                         += pair.E_INTERACTION
#                coarsed_mol.pairs.pop((imono+1,imono2+1))
#                #group.E_SELF                         += pair.E_NUC_NUC + pair.E_EL_EL + \
#                #                                        pair.E_NUC_EL  + pair.E_EL_NUC
#    #            for n in range(len(group.E_SELFEL_NUCLEI)):
#    #              if ( n < group.inuc-1 ): 
#    #                print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
#    #                         mono.inuc, n+1, mono.E_SELFEL_NUCLEI[n]))
#    #              else:
#    #                print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
#    #                         mono.inuc, n+2, mono.E_SELFEL_NUCLEI[n]))
#            for igroup2, other_groupidx in enumerate(other_groupsidx):
#                if groupidx[0] > other_groupidx[0]:
#                    #groupidx[0], other_groupidx[0] = other_groupidx[0], groupidx[0]
#                    groups_pair = coarsed_mol.pairs[(other_groupidx[0]+1,groupidx[0]+1)]
#                    assert groups_pair.ia == other_groupidx[0]+1 and groups_pair.ib == groupidx[0]+1
#                else:
#                    groups_pair = coarsed_mol.pairs[(groupidx[0]+1,other_groupidx[0]+1)]
#                    assert groups_pair.ia == groupidx[0]+1 and groups_pair.ib == other_groupidx[0]+1
#                #print(groups_pair.ia-1, groups_pair.ib-1)
#                for other_imono in other_groupidx:
#                    if imono != groupidx[0] or other_imono != other_groupidx[0]:
#                        if imono > other_imono:
#                            #imono, other_imono = other_imono, imono
#                            pair = self.pairs[(other_imono+1,imono+1)]
#                        else:
#                            pair = self.pairs[(imono+1,other_imono+1)]
#                        groups_pair.E_NUC_NUC                        += pair.E_NUC_NUC
#                        groups_pair.E_EL_NUC                         += pair.E_EL_NUC
#                        groups_pair.E_NUC_EL                         += pair.E_NUC_EL
#                        groups_pair.E_EL_EL_COULOMB                  += pair.E_EL_EL_COULOMB
#                        groups_pair.E_INTERACTION_CLASSICAL          += pair.E_INTERACTION_CLASSICAL
#                        groups_pair.E_EL_EL_EXCHANGE_CORRELATION     += pair.E_EL_EL_EXCHANGE_CORRELATION
#                        groups_pair.E_EL_EL_EXCHANGE                 += pair.E_EL_EL_EXCHANGE
#                        groups_pair.E_EL_EL_CORRELATION              += pair.E_EL_EL_CORRELATION
#                        groups_pair.E_EL_EL                          += pair.E_EL_EL
#                        groups_pair.E_INTERACTION                    += pair.E_INTERACTION
#                        groups_pair.DELOCALIZATION_INDEX             += pair.DELOCALIZATION_INDEX
#                        groups_pair.DELOCALIZATION_INDEX_EXCHANGE    += pair.DELOCALIZATION_INDEX_EXCHANGE
#                        groups_pair.DELOCALIZATION_INDEX_CORRELATION += pair.DELOCALIZATION_INDEX_CORRELATION
#                        groups_pair.N_ELECTRON_PAIRS                 += pair.N_ELECTRON_PAIRS
#                        groups_pair.E_EL_NUC_mp                      += pair.E_EL_NUC_mp
#                        groups_pair.E_NUC_EL_mp                      += pair.E_NUC_EL_mp
#                        groups_pair.E_EL_EL_COULOMB_mp               += pair.E_EL_EL_COULOMB_mp
#                        groups_pair.E_INTERACTION_CLASSICAL_mp       += pair.E_INTERACTION_CLASSICAL_mp
#                        groups_pair.E_EL_EL_EXCHANGE_CORRELATION_mp  += pair.E_EL_EL_EXCHANGE_CORRELATION_mp
#                        groups_pair.E_EL_EL_EXCHANGE_mp              += pair.E_EL_EL_EXCHANGE_mp
#                        groups_pair.E_EL_EL_CORRELATION_mp           += pair.E_EL_EL_CORRELATION_mp
#                        groups_pair.E_EL_EL_mp                       += pair.E_EL_EL_mp
#                        groups_pair.E_INTERACTION_mp                 += pair.E_INTERACTION_mp
#                        if imono > other_imono:
#                            coarsed_mol.pairs.pop((other_imono+1,imono+1))
#                        else:
#                            coarsed_mol.pairs.pop((imono+1,other_imono+1))
#
#
#                groups_pair.COULOMB_l              = []
#                groups_pair.EXCHANGE_CORRELATION_l = []
#                groups_pair.EXCHANGE_l             = []
#
#
#        group.rmins = np.nan # No meaning at all for a group
#        group.rmaxs = np.nan # No meaning at all for a group
#        group.LOCALIZATION_INDEX             = np.nan # groupelm.LOCALIZATION_INDEX
#        group.LOCALIZATION_INDEX_EXCHANGE    = np.nan # groupelm.LOCALIZATION_INDEX_EXCHANGE
#        group.LOCALIZATION_INDEX_CORRELATION = np.nan # groupelm.LOCALIZATION_INDEX_CORRELATION
#        group.ELEC_POP_MULTIPOLES_Qlm            += [{
#                                                          "tuple_element0": 0,
#                                                          "tuple_element1": 0,
#                                                          "tuple_element2": group.CHARGE
#                                                      }]
#        group.COULOMB_l              = []
#        group.EXCHANGE_CORRELATION_l = []
#        group.EXCHANGE_l             = []
#        
#
#    # Clean dispensable data
#    # 1) clean basins
#    # dynamic means that as an item is removed the indices shift -1
#    # this shift stacks up to the number of deleted items
#    dynamic_dispensable_basins = []
#    ndel = 0
#    for sublist in groupsidx:
#        for item in sublist[1:]:
#            dynamic_dispensable_basins.append(item-ndel)
#            ndel += 1
#    for dispensable_basin in dynamic_dispensable_basins:
#        del coarsed_mol.monos[dispensable_basin]
#    for igroup in range(len(groupsidx)):
#        coarsed_mol.monos[igroup].inuc = igroup+1
#
#    # Sort only items containing the final groups to have the user specified order
#    desired_group_order = [x[0] for x in groupsidx]
#    current_order = sorted(desired_group_order)
#    final_order = []
#    for idx in desired_group_order:
#      final_order.append(current_order.index(idx))
#    coarsed_mol.monos = list(itemgetter(*final_order)(coarsed_mol.monos))
#
#
#
#    # 2) clean interactions
#    other_groupsidx = copy.deepcopy(groupsidx)
#    for igroup, groupidx in enumerate(groupsidx):
#        del other_groupsidx[0]
#        for igroup2, other_groupidx in enumerate(other_groupsidx):
#            #print(igroup, igroup2)
#            #print("intergroup",igroup+1,igroup2+1+igroup+1)
#            if groupidx[0] > other_groupidx[0]:
#                #groupidx[0], other_groupidx[0] = other_groupidx[0], groupidx[0]
#                coarsed_mol.pairs[(igroup+1,igroup2+1+igroup+1)] = coarsed_mol.pairs.pop((other_groupidx[0]+1,groupidx[0]+1))
#            else:
#                coarsed_mol.pairs[(igroup+1,igroup2+1+igroup+1)] = coarsed_mol.pairs.pop((groupidx[0]+1,other_groupidx[0]+1))
#            coarsed_mol.pairs[(igroup+1,igroup2+1+igroup+1)].ia = igroup+1
#            coarsed_mol.pairs[(igroup+1,igroup2+1+igroup+1)].ib = igroup2+1+igroup+1
#
#    # Set total values and check consistency
#    coarsed_mol.CHARGE    = 0.0
#    coarsed_mol.E_KINETIC = 0.0
#    coarsed_mol.E_NUC_NUC = 0.0
#    for mono in coarsed_mol.monos:
#        coarsed_mol.CHARGE    += mono.CHARGE
#        coarsed_mol.E_KINETIC += mono.E_KINETIC
#        coarsed_mol.E_NUC_NUC += mono.E_SELF_NUC_NUC
#    for pair in coarsed_mol.pairs.values():
#      coarsed_mol.E_NUC_NUC += pair.E_NUC_NUC
#    assert coarsed_mol.CHARGE    - self.CHARGE    < 1e-4
#    assert coarsed_mol.E_KINETIC - self.E_KINETIC < 1e-4
#    assert coarsed_mol.E_NUC_NUC - self.E_NUC_NUC < 1e-4
#
#    print("The following groups have been created:")
#    coarsed_mol.content()
#
#    return coarsed_mol
    
  def coarsening(self, groupsidx, group_other=False):
    """
    Performs a coarse graining of basins
    E.g. (with groupsidx = [[1,2]], group_other = False )
    1 and 2 --> 1
    3       --> 2
    4       --> 3
    E.g. (with groupsidx = [[1,2]], group_other = True )
    1 and 2 --> 1
    3 and 4 --> 2
    Besides, reordering of fragments can be done
    E.g. (with groupsidx = [[2],[1]], group_other = False )
    1       --> 2
    2       --> 1
    3       --> 3
    4       --> 4
    """
    # Brighter solution:
    coarsed_mol = IQA_Molecule({})

    # Sort elements in each group (repeated elements are kept)
    for groupidx in groupsidx:
        groupidx.sort()
    # Fix groups to include all basins (exhaustive)
    # and change indexing from [1,] to [0,]
    basins = [item for sublist in groupsidx for item in sublist]
    # check that it is sorted in strict increasing order (no repeated elements)
    basins_set = set(basins)
    for basin in basins_set:
        if basins.count(basin) > 1:
            raise Exception("{:4d} is duplicated".format(basin))

    #assert all(basins[i] < basins[i+1] for i in range(len(basins)-1)) 
    all_basins = set(range(1,len(self.monos)+1))
    other_basins = list(all_basins.difference(basins_set))
    #assert all(other_basins[i] < other_basins[i+1] for i in range(len(other_basins)-1)) 
    if group_other:
        groupsidx.append( sorted(other_basins) )
    else:
        for basin in sorted(other_basins):
            groupsidx.append( [basin] )
    coarsed_mol.coarseness = copy.deepcopy(groupsidx)

    logging.debug("Interpreted coarsening: {}".format(groupsidx))
    for n in range(len(groupsidx)):
        groupsidx[n] = [ x - 1 for x in groupsidx[n]]

    other_groupsidx = copy.deepcopy(groupsidx)
    for igroup, groupidx in enumerate(groupsidx):
        del other_groupsidx[0]
        coarsed_mol.monos.append(copy.deepcopy(self.monos[groupidx[0]]))
        group = coarsed_mol.monos[igroup]
        if not hasattr(group, "E_SELF_NUC_NUC"):
            setattr(group, "E_SELF_NUC_NUC", 0.0)
        logging.debug("Adding internal terms of group {:8d}...".format(igroup+1))
        for imono in groupidx:
            # Sum intra terms to the group
            logging.debug("                                          {:4d}       (mono )".format(imono+1))
            if imono != groupidx[0]:
                groupelm = self.monos[imono]
                group.label                        += " " + groupelm.label
                group.E_KINETIC                    += groupelm.E_KINETIC
                group.Z                            += groupelm.Z
                group.CHARGE                       += groupelm.CHARGE
                group.N_ELECTRONS                  += groupelm.N_ELECTRONS

                group.E_NUC_EL                     += groupelm.E_NUC_EL
                group.E_EL_EL_COULOMB              += groupelm.E_EL_EL_COULOMB
                group.E_EL_EL_EXCHANGE_CORRELATION += groupelm.E_EL_EL_EXCHANGE_CORRELATION
                group.E_EL_EL_EXCHANGE             += groupelm.E_EL_EL_EXCHANGE
                group.E_EL_EL_CORRELATION          += groupelm.E_EL_EL_CORRELATION
                group.E_EL_EL                      += groupelm.E_EL_EL
                group.E_SELF                       += groupelm.E_SELF

                group.LOCALIZATION_INDEX             += groupelm.LOCALIZATION_INDEX
                group.LOCALIZATION_INDEX_EXCHANGE    += groupelm.LOCALIZATION_INDEX_EXCHANGE
                group.LOCALIZATION_INDEX_CORRELATION += groupelm.LOCALIZATION_INDEX_CORRELATION

            # Sum interactions within the group
            for imono2 in groupidx:
                if imono >= imono2: # Do not count interactions twice
                    continue        # elements within the group must be sorted
                logging.debug("                                          {:4d} {:4d}  (bicen)".format(imono+1,imono2+1))
                pair = self.pairs[(imono+1,imono2+1)]
                assert pair.ia == imono+1
                assert pair.ib == imono2+1
                group.E_EL_EL_COULOMB                += pair.E_EL_EL_COULOMB
                group.E_EL_EL_EXCHANGE_CORRELATION   += pair.E_EL_EL_EXCHANGE_CORRELATION
                group.E_EL_EL_EXCHANGE               += pair.E_EL_EL_EXCHANGE
                group.E_EL_EL_CORRELATION            += pair.E_EL_EL_CORRELATION
                group.E_EL_EL                        += pair.E_EL_EL
                group.E_SELF_NUC_NUC                 += pair.E_NUC_NUC
                group.E_NUC_EL                       += pair.E_NUC_EL + pair.E_EL_NUC
                group.E_SELF                         += pair.E_INTERACTION
                #group.E_SELF                         += pair.E_NUC_NUC + pair.E_EL_EL + \
                #                                        pair.E_NUC_EL  + pair.E_EL_NUC
                group.LOCALIZATION_INDEX             += 0.5 * pair.DELOCALIZATION_INDEX
                group.LOCALIZATION_INDEX_EXCHANGE    += 0.5 * pair.DELOCALIZATION_INDEX_EXCHANGE
                group.LOCALIZATION_INDEX_CORRELATION += 0.5 * pair.DELOCALIZATION_INDEX_CORRELATION
                # TODO: add E_SELFEL_NUCLEI
    #            for n in range(len(group.E_SELFEL_NUCLEI)):
    #              if ( n < group.inuc-1 ): 
    #                print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
    #                         mono.inuc, n+1, mono.E_SELFEL_NUCLEI[n]))
    #              else:
    #                print("   +++ energy elec({0:3d})-nuc({1:3d})  = {2:+16.8f}".format(
    #                         mono.inuc, n+2, mono.E_SELFEL_NUCLEI[n]))
        for igroup2, other_groupidx in enumerate(other_groupsidx):
            logging.debug("Adding interactions of groups {:4d} {:4d}...".format(igroup+1,igroup+1+igroup2+1))
            if groupidx[0] > other_groupidx[0]:
                #groupidx[0], other_groupidx[0] = other_groupidx[0], groupidx[0]
                #coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)] = copy.deepcopy(self.pairs[(other_groupidx[0]+1,groupidx[0]+1)])
                #coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)] = self.pairs.pop((other_groupidx[0]+1,groupidx[0]+1))
                coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)] = copy.deepcopy(self.pairs[(other_groupidx[0]+1,groupidx[0]+1)])
                groups_pair = coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)]
                assert groups_pair.ia == other_groupidx[0]+1 and groups_pair.ib == groupidx[0]+1
                logging.debug("                                          {:4d} {:4d}".format(other_groupidx[0]+1,groupidx[0]+1))
            else:
                #coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)] = self.pairs.pop((groupidx[0]+1,other_groupidx[0]+1))
                coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)] = copy.deepcopy(self.pairs[(groupidx[0]+1,other_groupidx[0]+1)])
                groups_pair = coarsed_mol.pairs[(igroup+1,igroup+1+igroup2+1)]
                assert groups_pair.ia == groupidx[0]+1 and groups_pair.ib == other_groupidx[0]+1
                logging.debug("                                          {:4d} {:4d}".format(groupidx[0]+1,other_groupidx[0]+1))
            groups_pair.ia = igroup+1
            groups_pair.ib = igroup+1+igroup2+1
            for imono in groupidx:
                for other_imono in other_groupidx:
                    if imono != groupidx[0] or other_imono != other_groupidx[0]:
                        if imono > other_imono:
                            #imono, other_imono = other_imono, imono
                            pair = self.pairs[(other_imono+1,imono+1)]
                            logging.debug("                                          {:4d} {:4d}".format(other_imono+1,imono+1))
                        else:
                            pair = self.pairs[(imono+1,other_imono+1)]
                            logging.debug("                                          {:4d} {:4d}".format(imono+1,other_imono+1))
                        groups_pair.E_NUC_NUC                        += pair.E_NUC_NUC
                        groups_pair.E_EL_NUC                         += pair.E_EL_NUC
                        groups_pair.E_NUC_EL                         += pair.E_NUC_EL
                        groups_pair.E_EL_EL_COULOMB                  += pair.E_EL_EL_COULOMB
                        groups_pair.E_INTERACTION_CLASSICAL          += pair.E_INTERACTION_CLASSICAL
                        groups_pair.E_EL_EL_EXCHANGE_CORRELATION     += pair.E_EL_EL_EXCHANGE_CORRELATION
                        groups_pair.E_EL_EL_EXCHANGE                 += pair.E_EL_EL_EXCHANGE
                        groups_pair.E_EL_EL_CORRELATION              += pair.E_EL_EL_CORRELATION
                        groups_pair.E_EL_EL                          += pair.E_EL_EL
                        groups_pair.E_INTERACTION                    += pair.E_INTERACTION
                        groups_pair.DELOCALIZATION_INDEX             += pair.DELOCALIZATION_INDEX
                        groups_pair.DELOCALIZATION_INDEX_EXCHANGE    += pair.DELOCALIZATION_INDEX_EXCHANGE
                        groups_pair.DELOCALIZATION_INDEX_CORRELATION += pair.DELOCALIZATION_INDEX_CORRELATION
                        groups_pair.N_ELECTRON_PAIRS                 += pair.N_ELECTRON_PAIRS
                        groups_pair.E_EL_NUC_mp                      += pair.E_EL_NUC_mp
                        groups_pair.E_NUC_EL_mp                      += pair.E_NUC_EL_mp
                        groups_pair.E_EL_EL_COULOMB_mp               += pair.E_EL_EL_COULOMB_mp
                        groups_pair.E_INTERACTION_CLASSICAL_mp       += pair.E_INTERACTION_CLASSICAL_mp
                        groups_pair.E_EL_EL_EXCHANGE_CORRELATION_mp  += pair.E_EL_EL_EXCHANGE_CORRELATION_mp
                        groups_pair.E_EL_EL_EXCHANGE_mp              += pair.E_EL_EL_EXCHANGE_mp
                        groups_pair.E_EL_EL_CORRELATION_mp           += pair.E_EL_EL_CORRELATION_mp
                        groups_pair.E_EL_EL_mp                       += pair.E_EL_EL_mp
                        groups_pair.E_INTERACTION_mp                 += pair.E_INTERACTION_mp


                groups_pair.COULOMB_l              = []
                groups_pair.EXCHANGE_CORRELATION_l = []
                groups_pair.EXCHANGE_l             = []


        group.rmins = np.nan # No meaning at all for a group
        group.rmaxs = np.nan # No meaning at all for a group
        group.LOCALIZATION_INDEX             = np.nan # groupelm.LOCALIZATION_INDEX
        group.LOCALIZATION_INDEX_EXCHANGE    = np.nan # groupelm.LOCALIZATION_INDEX_EXCHANGE
        group.LOCALIZATION_INDEX_CORRELATION = np.nan # groupelm.LOCALIZATION_INDEX_CORRELATION
        # TODO: Sum rule for population
        group.ELEC_POP_MULTIPOLES_Qlm            += [{
                                                          "tuple_element0": 0,
                                                          "tuple_element1": 0,
                                                          "tuple_element2": group.CHARGE
                                                      }]
        group.COULOMB_l              = []
        group.EXCHANGE_CORRELATION_l = []
        group.EXCHANGE_l             = []
        


    # Set total values and check consistency
    #coarsed_mol.CHARGE    = 0.0
    #coarsed_mol.E_KINETIC = 0.0
    #coarsed_mol.E_NUC_NUC = 0.0
    for mono in coarsed_mol.monos:
        coarsed_mol.CHARGE    += mono.CHARGE
        coarsed_mol.E_KINETIC += mono.E_KINETIC
        coarsed_mol.E_NUC_NUC += mono.E_SELF_NUC_NUC
        coarsed_mol.E_NUC_EL += mono.E_NUC_EL
        coarsed_mol.E_EL_EL   += mono.E_EL_EL
        coarsed_mol.E_EL_EL_COULOMB   += mono.E_EL_EL_COULOMB
        coarsed_mol.E_EL_EL_EXCHANGE_CORRELATION += mono.E_EL_EL_EXCHANGE_CORRELATION
        coarsed_mol.E_EL_EL_EXCHANGE += mono.E_EL_EL_EXCHANGE
        coarsed_mol.E_EL_EL_CORRELATION += mono.E_EL_EL_CORRELATION
        coarsed_mol.E_NET += mono.E_SELF
    for pair in coarsed_mol.pairs.values():
      coarsed_mol.E_NUC_NUC += pair.E_NUC_NUC
      coarsed_mol.E_NUC_EL += pair.E_NUC_EL + pair.E_EL_NUC
      coarsed_mol.E_EL_EL   += pair.E_EL_EL
      coarsed_mol.E_EL_EL_COULOMB   += pair.E_EL_EL_COULOMB
      coarsed_mol.E_EL_EL_EXCHANGE_CORRELATION += pair.E_EL_EL_EXCHANGE_CORRELATION
      coarsed_mol.E_EL_EL_EXCHANGE += pair.E_EL_EL_EXCHANGE
      coarsed_mol.E_EL_EL_CORRELATION += pair.E_EL_EL_CORRELATION
      coarsed_mol.E_INTERACTION     += pair.E_INTERACTION
      coarsed_mol.E_INTERACTION_CLASSICAL     += pair.E_INTERACTION_CLASSICAL

    coarsed_mol.E_POTENTIAL= coarsed_mol.E_NUC_NUC + coarsed_mol.E_NUC_EL + coarsed_mol.E_EL_EL
    coarsed_mol.E_TOTAL= coarsed_mol.E_KINETIC + coarsed_mol.E_POTENTIAL
    coarsed_mol.TWO_KINETIC_PLUS_POTENTIAL = 2*coarsed_mol.E_KINETIC + coarsed_mol.E_POTENTIAL
    coarsed_mol.VIRIAL_RATIO =  coarsed_mol.E_POTENTIAL /coarsed_mol.E_KINETIC 
    coarsed_mol.E_INTERACTION_EXCHANGE_CORRELATION += coarsed_mol.E_INTERACTION - coarsed_mol.E_INTERACTION_CLASSICAL

    assert abs(coarsed_mol.CHARGE    - self.CHARGE)    < 1e-4, "The total charge is not preserved"
    assert abs(coarsed_mol.E_KINETIC - self.E_KINETIC) < 1e-4, "The total kinetic energy is not preserved"
    assert abs(coarsed_mol.E_NUC_NUC - self.E_NUC_NUC) < 1e-4, "The total nuclear-nuclear repulsion is not preserved"
    assert abs(coarsed_mol.E_NUC_EL - self.E_NUC_EL) < 1e-4, "The total nuclear-electron energy is not preserved"
    assert abs(coarsed_mol.E_EL_EL_COULOMB - self.E_EL_EL_COULOMB) < 1e-4, "The total Coulomb energy is not preserved"
    assert abs(coarsed_mol.E_EL_EL_EXCHANGE - self.E_EL_EL_EXCHANGE) < 1e-4, "The total exchange energy is not preserved"
    assert abs(coarsed_mol.E_EL_EL_EXCHANGE_CORRELATION - self.E_EL_EL_EXCHANGE_CORRELATION) < 1e-4, "The total exchange-correlation energy is not preserved"
    assert abs(coarsed_mol.E_EL_EL_CORRELATION - self.E_EL_EL_CORRELATION) < 1e-4, "The total correlation energy is not preserved"
    assert abs(coarsed_mol.E_POTENTIAL - self.E_POTENTIAL) < 1e-4, "The total potential energy is not preserved"
    assert abs(coarsed_mol.E_TOTAL - self.E_TOTAL) < 1e-4, "The total energy is not preserved"

    print("The following groups have been created:")
    coarsed_mol.content()

    return coarsed_mol
    
  def content(self):

    print("{:-^5}   {:-^30}".format("",""))
    print("{:5s}   {:^30s}".format("Group","Basins in each group"))
    print("{:-^5}   {:-^30}".format("",""))
    for igroup, group in enumerate(self.monos):
        label = textwrap.fill(group.label, width=80,initial_indent='', subsequent_indent=' '*8)
        print("{:5d}   {}".format(igroup+1,label))
    print("{:-^5}   {:-^30}".format("",""))
    print("NOTE: Basin indices are allways the same as in the DGRID_WFN file.")

  def get_fragment(self,i):

    return self.monos[i-1]
  
  def get_interaction(self,ia,ib):

    return self.pairs[(ia,ib)]

    


m = IQA_Molecule( data );

m.monos[0].__dict__.keys()
m.pairs[(1,2)].__dict__.keys()


#m.demo(plt)


