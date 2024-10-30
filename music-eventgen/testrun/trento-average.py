#!/usr/bin/env python3 

import os
import numpy as np
import sys
import subprocess
directory=str(sys.argv[1])
files=os.listdir(directory)
df=np.loadtxt(directory+files[0],skiprows=8)
files_num=len(files)
if int(files_num) == int(sys.argv[2]):
	for file in files[1:]:
		df+=np.loadtxt(directory+file,skiprows=8)
	df=df/files_num
	np.savetxt("./trento_data/average.dat",df)
else:
	print("Extra trento input files")


#transforms averaged files to MUSIC format

in_name="./trento_data/average.dat"
out_name="./trento_data/music-input-from-trento.dat"

grid_n=int(sys.argv[3])
grid_step=float(sys.argv[4])

header_bash_command = "echo '# dummy 0.0 etamax= 1 xmax= "+str(grid_n)+" ymax= "+str(grid_n)+" deta= 10. dx= "+str(grid_step)+" dy= "+str(grid_step)+"' > "+out_name
subprocess.run(header_bash_command, shell=True)

conversion_bash_command="Ns=\""+str(grid_n)+"\"; tail -n+9 "+in_name +" | perl -pe 's/\s+/\\n/g' | awk -v Ns=${Ns} 'BEGIN {for(x=0;x<Ns;x++){for(y=0;y<Ns;y++){T00[x+Ns*y]=0.0;}} Index=0;} {T00[Index]=$1; Index++;} END {for(x=0;x<Ns;x++){for(y=0;y<Ns;y++){print 0, (x-"+str(grid_n/2.)+")*"+str(grid_step)+","+"(y-"+str(grid_n/2.)+")*"+str(grid_step)+",T00[x+Ns*y],1,0,0,0, 0,0,0,0,0,0,0,0,0,0;} }}' >> "+out_name
subprocess.run(conversion_bash_command, shell=True)
