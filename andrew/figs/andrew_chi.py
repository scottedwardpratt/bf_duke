import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

mydata = np.loadtxt('eos_vs_epsilon.txt',skiprows=1,unpack=True)
epsilon=mydata[0]
T=1000*mydata[1]
s=mydata[2]
chiuu=mydata[3]
chiud=mydata[4]
chius=mydata[5]
chiss=mydata[6]

mydata = np.loadtxt('eos_hadron_vs_epsilon.txt',skiprows=1,unpack=True)
epsilon_hadron=mydata[0]
T_hadron=1000*mydata[1]
s_hadron=mydata[2]
chiuu_hadron=mydata[3]
chiud_hadron=mydata[4]
chius_hadron=mydata[5]
chiss_hadron=mydata[6]

mydata = np.loadtxt('eos_lattice_vs_epsilon.txt',skiprows=1,unpack=True)
epsilon_lattice=mydata[0]
T_lattice=1000*mydata[1]
s_lattice=mydata[2]
chiuu_lattice=mydata[3]
chiud_lattice=mydata[4]
chius_lattice=mydata[5]
chiss_lattice=mydata[6]

plt.plot(T,chiuu,linestyle='-',linewidth=4,color='r',markersize=3, marker=None)
plt.plot(T_hadron ,chiuu_hadron,linestyle='-',linewidth=2,color='g')
plt.plot(T_lattice ,chiuu_lattice,linestyle='-',linewidth=2,color='b')

#plt.plot(T,chiud/s,linestyle='-',linewidth=2,color='r',markersize=3, marker='o', markerfacecolor='r', markeredgecolor='r')
#plt.plot(T,chius/s,linestyle='-',linewidth=2,color='b',markersize=3, marker='o', markerfacecolor='b', markeredgecolor='b')
#plt.plot(T,chiss/s,linestyle='-',linewidth=2,color='k',markersize=3, marker='o', markerfacecolor='k', markeredgecolor='k')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_yticks(np.arange(-40.0,80.0,2), minor=False)
ax.set_yticklabels(np.arange(-40.0,80.0,2), minor=False, family='serif')
ax.set_yticks(np.arange(-40.0,80.0,1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.ylim(-1,3)

ax.set_xticks(np.arange(0,500,50), minor=False)
ax.set_xticklabels(np.arange(0,500,50), minor=False, family='serif')
ax.set_xticks(np.arange(0,500,10), minor=True)
plt.xlim(75,410.0)
#ax.set_xticks(0.1:1.0:10.0:100.0, minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=18, weight='normal')
plt.ylabel('$\chi$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('chi_andrew.pdf',format='pdf')
os.system('open -a Preview chi_andrew.pdf')
#plt.show()
quit()
