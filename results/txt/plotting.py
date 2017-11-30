import re
import os
from subprocess import Popen, PIPE
from matplotlib.pyplot import *
from numpy import *


def sorting(txtfiles):
    i=0
    new_txtfiles = txtfiles[:]
    temperatures = []
    for file in txtfiles:

        #print data.keys()
        start = file.find("T")
        stopp = file.find(".txt")
        temperatureK = float(file[start+2:stopp])
        temperatures.append(temperatureK)
        new_txtfiles[i] = temperatureK
        i+=1

    sorts = argsort(new_txtfiles)
    ny = [txtfiles[sorts[i]] for i in range(len(sorts))]
    temps = [int(temperatures[sorts[i]]) for i in range(len(sorts))]
    return ny, temps



output = Popen(["ls"], stdout=PIPE).communicate()[0]

txtfiles = re.findall(".*\.txt",output,re.IGNORECASE)

files, temperatures = sorting(txtfiles)
print files

i=-1
for file in files:
    with open(file,"r") as infile:
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

        print data.keys()

        figure(1)

        plot(data["timesteps"],data["TotalE"], label = '%0.f'%temperatures[i])
        figure(2)
        plot(data["timesteps"],data["Temperature"], label = '%0.f'%temperatures[i])

legend()
show()
