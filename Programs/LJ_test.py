from __future__ import division 
from numpy import *
from scipy.integrate import quad 
from matplotlib.pyplot import *


def LJ(rij):
	rij = float(rij)
	m_epsilon = 1.10054
	sigma = 3.405
	F = 24*m_epsilon*( 2* (sigma/rij)**12 -(sigma/rij)**6   )/rij
	U = 4*m_epsilon*( (sigma/rij)**12 -(sigma/rij)**6   )

	return U



sigma = 3.405
figure()
r = linspace(1,10,1000)
U = [LJ(r[i]) for i in range(len(r))]
plot(r,U, label = 'Lennard Jones potentital')
xlabel('Distance between particles, $\AA$')
ylabel('Energy, $1.651e-21 J$')
legend()
title('Lennard Jones potentital')
#axvline(2*sigma)
ylim([-1.5,10])
grid('on')
savefig('../figures/LJpot.pdf')

show()