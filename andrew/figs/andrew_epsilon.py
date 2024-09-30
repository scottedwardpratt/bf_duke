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

plt.plot(T,epsilon,linestyle='-',linewidth=2,color='r',markersize=3, marker='o', markerfacecolor='g', markeredgecolor='g')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_yticks(np.arange(0,61,10), minor=False)
ax.set_yticklabels(np.arange(0,61,10), minor=False, family='serif')
ax.set_yticks(np.arange(0,61,5), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.ylim(0.0,43)

ax.set_xticks(np.arange(0,500,50), minor=False)
ax.set_xticklabels(np.arange(0,500,50), minor=False, family='serif')
ax.set_xticks(np.arange(0,500,10), minor=True)
plt.xlim(75,410.0)
#ax.set_xticks(0.1:1.0:10.0:100.0, minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=18, weight='normal')
plt.ylabel('$\epsilon$ (GeV/fm$^3$)',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('epsilon_andrew.pdf',format='pdf')
os.system('open -a Preview epsilon_andrew.pdf')
#plt.show()
quit()
