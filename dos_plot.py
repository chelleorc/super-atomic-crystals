import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np


# load data
energy, dos, idos = np.loadtxt('mos2_dos.dat', unpack=True)

# make plot
plt.figure(figsize = (12, 6))
plt.plot(energy, dos, linewidth=0.75, color='red')
plt.yticks([])
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
plt.axvline(x=6.642, linewidth=0.5, color='k', linestyle=(0, (8, 10)))
plt.xlim(-16, 6)
plt.ylim(0, )

# get Fermi energy from non-scf output file
plt.fill_between(energy, 0, dos, where=(energy < 0.9206), facecolor='red', alpha=0.25)
plt.text(6, 1.7, 'Fermi energy', fontsize= 'medium', rotation=90)
plt.title("DOS vs. Energy(eV)")
plt.savefig('mos2_dos.png')
plt.show()