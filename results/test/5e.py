import re
import os
import operator
from subprocess import Popen, PIPE
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

from numpy import *

#energy, probability = loadtxt("probability.txt", unpack=True, skiprows=1)

output = Popen(["ls"], stdout=PIPE).communicate()[0]
txtfiles = re.findall("5_T_.*\.txt",output,re.IGNORECASE)
print txtfiles
i = -1

kb = 1.3806488e-23
E0eV = 1.0318e-2     
E0 = 1.60217657e-19*E0eV
T0 = E0/kb

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

for txtfile in txtfiles:
    with open(txtfile,"r") as infile:
        i+=1
        data = {}
        #infile.readline()
        identifiers = infile.readline().split()
        for identifier in identifiers:
            data[identifier] = []

        lines = infile.readlines()
        for line in lines: 
            values = line.split()
            for identifier,value in zip(identifiers,values):
                data[identifier].append(float(value))

        for key in data.keys():
            data[key] = array(data[key])

        #print data.keys()

        start = txtfile.find("T_")
        stopp = txtfile.find(".txt")
        labell = txtfile[start+2:stopp]
        print labell

        figure(1, figsize=(9,7))
       
        plot(data["t"],data["Temperature"]/data["Temperature"][0], label = "T$_{ini}$ = " + labell + " K")

ax = subplot(111)       

title("Development of temperature")
ylabel(r"$T/T_{ini}$")
xlabel("t")
xlim([0,1.5e-12])
# Shrink current axis's height by 10% on the bottom
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                 box.width, box.height * 0.9])
#
## Put a legend below current axis
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.085),
#          fancybox=True, shadow=True, ncol=5, fontsize=12)
legend(fancybox=True, framealpha=0.5)
savefig("../../figures/temp_development.pdf")
show()
