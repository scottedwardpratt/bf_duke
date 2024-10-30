#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy
pardict={}

run_number=sys.argv[1]
qualifier=sys.argv[2]

#run_number=input("Enter run number: ")
#qualifier=input("Enter qualifier: (e.g. alicePbPb_cent10_15) ")

if qualifier == 'alicePbPb_cent0_5':
  b=2.467
elif qualifier == 'alicePbPb_cent5_10':
  b=4.272
elif qualifier == 'alicePbPb_cent10_15':
  b=5.515
elif qualifier == 'alicePbPb_cent15_20':
  b=6.526
else:
  print('qualifier='+qualifier+', not recognized\n')
  exit()

parfilename="hydrodata/modelruns/run"+str(run_number)+"/mod_parameters.txt"


if os.path.isfile(parfilename):
  print("model parameter file="+parfilename+" exists\n")
else:
  print("model parameter file="+parfilename+" does not exist!\n")
  exit()

parfile = open(parfilename, "r")

for oneline in parfile:
  words=oneline.split()
  key=words[0]
  value=words[1]
  pardict[words[0]]=words[1]

print(pardict)

visc_shear_min=float(pardict['VISCOSITY_SHEARMIN'])
#SLOPE1=(0.16-visc_shear_min)/(0.158-0.17)
SHEAR_SLOPE1=0.0

#command='./trento-music.sh'+' '+run_number+' '+qualifier+' '+str(b)+' '+pardict['VISCOSITY_BULKMAX']+' '+pardict['VISCOSITY_SHEARMIN']+' '+str(SHEAR_SLOPE1)+' ' +pardict['VISCOSITY_SHEARSLOPE']+' '+pardict['TRENTO_NORM']+' '+pardict['TRENTO_P']

command='./trento-music.sh'+' '+run_number+' '+qualifier+' '+str(b)+' '+pardict['VISCOSITY_SHEARMIN']+' '+pardict['TRENTO_NORM']+' '+pardict['TRENTO_P']

print('command='+command)

subprocess.run(command, shell=True)
