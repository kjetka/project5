from matplotlib.pyplot import *
from scipy.optimize import curve_fit

from numpy import *


def LJ(x):
	epsilon = 1.00054
	sigma = 3.405
	return 4*epsilon*(	(sigma/x)**12 - (sigma/x)**6	)



def F(x):
	epsilon = 1.00054
	sigma = 3.405
	return 24*epsilon*(	2*(sigma/x)**12 - (sigma/x)**6	)*1/x

r = linspace(1,10,10000)
ev = 1.60217662e-19
MD = 1.651e-21
print ev/MD

plot(r,LJ(r), label = 'Potential')
plot(r, F(r), label = 'Force')
ylim([-2,2])
legend()
grid('on')
xlabel('atomistic distance $r_{ij}$, $\AA$')
ylabel('Energy, MD energy')
title('Lennard Jones potential for Ar')
savefig('LJ.pdf')
