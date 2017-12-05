from numpy import *
from matplotlib.pyplot import *

timesteps,t,Temperature,Kin,Pot,TotalE, MSD = loadtxt("txt/5d_T_100.txt", unpack=True, skiprows=1)


plot(t, MSD)
title("Energy probability at T = 2.4 K" )
ylabel(r"$<r(t)^2>$")
xlabel("t")
#xticks(sorted(energy[:-1]))
savefig("diffusion_T_1000.pdf")
show()


        
        
        