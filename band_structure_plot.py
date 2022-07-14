import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np


plt.rcParams["figure.dpi"]=150
plt.rcParams["figure.facecolor"]="white"
plt.rcParams["figure.figsize"]=(8, 6)

# load data
data = np.loadtxt('MoS2_bands.dat.gnu')
print(data)
# first column of dat.gnu band data
k = np.unique(data[:, 0])


# second column of dat.gnu band data
bands = np.reshape(data[:, 1], (-1, len(k)))
for band in range(len(bands)):
    plt.plot(k, bands[band, :], linewidth=1, alpha=0.5, color='k')

plt.xlim(min(k), max(k))
print(min(k))
print(max(k))

# Set Fermi energy level to 0eV (check bands.x output file)
plt.axhline(0.0000, linestyle=(0, (5, 5)), linewidth=0.75, color='k', alpha=0.5)
# High symmetry k-points (check bands.x output file)
plt.axvline(0.0000, linewidth=0.75, color='k', alpha=0.5)
plt.axvline(0.5774, linewidth=0.75, color='k', alpha=0.5)
plt.axvline(1.0865, linewidth=0.75, color='k', alpha=0.5)
plt.axvline(2.1049, linewidth=0.75, color='k', alpha=0.5)


# text labels and ticks at each high symmetry point on x-axis
plt.xticks(ticks= [0.0000,0.5774,1.0865,2.1049], \
           labels=['$\Gamma$','M', 'K','$\Gamma$'])
plt.ylabel("Energy (eV)")
plt.text(1.3, -0.3, 'Fermi energy', fontsize='small')

# plt.ylim([-2,5])
# save file as png
plt.savefig('MoS2_bands.png')
plt.show()