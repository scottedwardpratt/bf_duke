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

rcParams['lines.markersize'] ** 1

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
sdens=mydata[2]
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

mydata = np.loadtxt('eos_lattice_vs_epsilon.txt',skiprows=0,unpack=True)
epsilon_lattice=mydata[0]
T_lattice=1000*mydata[1]
s_lattice=mydata[2]
chiuu_lattice=mydata[3]
chiud_lattice=mydata[4]
chius_lattice=mydata[5]
chiss_lattice=mydata[6]

xline=[158,158]
yline=[-0.1,0.15]
plt.plot(xline,yline,linestyle='--',color='k',lw='2')
msize=2.0

plt.plot(T,chiuu/sdens,linestyle='-',linewidth=4,color='r',markersize=3, marker=None)
plt.scatter(T_hadron,chiuu_hadron/s_hadron,msize,color='r',marker='o')
plt.scatter(T_lattice ,chiuu_lattice/s_lattice,msize,color='r')

plt.plot(T,chiud/sdens,linestyle='-',linewidth=4,color='g',markersize=3, marker=None)
plt.scatter(T_hadron,chiud_hadron/s_hadron,msize,color='g',marker='o')
plt.scatter(T_lattice ,chiud_lattice/s_lattice,msize,color='g')

plt.plot(T,chius/sdens,linestyle='-',linewidth=4,color='b',markersize=3, marker=None)
plt.scatter(T_hadron,chius_hadron/s_hadron,msize,color='b',marker='o')
plt.scatter(T_lattice ,chius_lattice/s_lattice,msize,color='b')

plt.plot(T,chiss/sdens,linestyle='-',linewidth=4,color='k',markersize=3, marker=None)
plt.scatter(T_hadron,chiss_hadron/s_hadron,msize,color='k',marker='o')
plt.scatter(T_lattice ,chiss_lattice/s_lattice,msize,color='k')


#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_yticks(np.arange(-1.0,1.0,0.05), minor=False)
ax.set_yticklabels(np.arange(-1.0,1.0,0.05), minor=False, family='serif')
ax.set_yticks(np.arange(-1.0,1.0,0.01), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.ylim(-0.1,0.15)

ax.set_xticks(np.arange(0,500,50), minor=False)
ax.set_xticklabels(np.arange(0,500,50), minor=False, family='serif')
ax.set_xticks(np.arange(0,500,10), minor=True)
plt.xlim(120,410.0)
#ax.set_xticks(0.1:1.0:10.0:100.0, minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=18, weight='normal')
plt.ylabel('$\chi/s$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('chiratio_andrew.pdf',format='pdf')
os.system('open -a Preview chiratio_andrew.pdf')
#plt.show()
quit()
