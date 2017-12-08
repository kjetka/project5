import re
import os
import operator
from subprocess import Popen, PIPE
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

from numpy import *

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

output = Popen(["ls"], stdout=PIPE).communicate()[0]
txtfiles = re.findall("nearTc_T_.*\.txt",output,re.IGNORECASE)
print txtfiles
i = -1

markers = ["o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "*", "h", "H", "+", "x","X","D"]
colors = ["b", "g", "r", "c", "m", "y","b", "g", "r", "c", "m", "y","b", "g", "r", "c", "m", "y"]

diffusion = zeros(len(txtfiles))
temp = zeros(len(txtfiles))

kb = 1.3806488e-23
E0eV = 1.0318e-2     
E0 = 1.60217657e-19*E0eV
T0 = E0/kb

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
        
        a,b = polyfit(data["t"][-len(data["t"])/2:], data["MSD"][-len(data["t"])/2:], 1)
        D_= a/6*1e-16 #from AA to cm        
        #D_= a/6
        D = "%0.3g"%D_        
        
                
        
        diffusion[i] = D_
        
        temp[i] = sum(data["Temperature"][500:600])/100.0
        
        figure(1,figsize=(8,6)) 
        plot(data["t"], data["MSD"], colors[i]+markers[i], Markersize=4,markeredgecolor=colors[i],label = r"T = %.0f K"%temp[i])
        plot(data["t"], a*data["t"]+b, colors[i], label = "D = "+ D + " cm$^2$/s")


ax2 = subplot(111)

figure(1)

title("Displacement - diffusivity")
ylabel(r"$<r^2(t)> \,[\AA^2]$",fontsize=16)
xlabel("time [s]")
legend(loc=2,fontsize=12, fancybox=True, framealpha=0.5)
savefig("../../figures/nearTc.pdf")

figure(2, figsize=(8,6))
title("Diffusion constant")
gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
ylabel("D [cm$^2$/s]")
xlabel("temperature [K]")
plot(temp, diffusion, 'o')
savefig("../../figures/diffusion_temp_nearTc.pdf")

show()






