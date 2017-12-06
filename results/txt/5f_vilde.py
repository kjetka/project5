from numpy import *
from matplotlib.pyplot import *

Ts = [100, 200, 300, 400, 500]

for T in Ts:
    timesteps, t, Temperature, Kin, Pot, TotalE, MSD = loadtxt("d_T_%.0g.txt"%T, unpack=True, skiprows=1)

title("Diffusion constant" )
ylabel(r"$<r^2(t)>$")
xlabel("t")
#xticks(sorted(energy[:-1]))
savefig("diffusion_constants.pdf")
show()