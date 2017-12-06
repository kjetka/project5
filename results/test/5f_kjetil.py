import re
import os
import operator
from subprocess import Popen, PIPE
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

from numpy import *

#energy, probability = loadtxt("probability.txt", unpack=True, skiprows=1)

output = Popen(["ls"], stdout=PIPE).communicate()[0]
txtfiles = re.findall("5d_T_.*\.txt",output,re.IGNORECASE)
print txtfiles
i = -1

markers = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+", "x","X","D"]

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

        figure(1)
        #plot(data["Temperatures"][3:-5], data["E_avg"][3:-5], label = "L = " + labell)
        plot(data["t"], data["MSD"], markers[i], label = "T = " + labell)
        
 


legend()
ylim([0,2])
show()


