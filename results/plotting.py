import re
import os
from subprocess import Popen, PIPE
from matplotlib.pyplot import *
from numpy import *

#energy, probability = loadtxt("probability.txt", unpack=True, skiprows=1)


output = Popen(["ls"], stdout=PIPE).communicate()[0]
txtfiles = re.findall(".*\.txt",output,re.IGNORECASE)
print txtfiles

for txtfile in txtfiles:
    with open(txtfile,"r") as infile:
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

        print data.keys()
       
     