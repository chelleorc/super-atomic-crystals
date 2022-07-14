import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import seaborn as sb


# list of headers for each data column
# col = ['energy cutoff (Ry)', 'total energy (Ry)', 'wall time (s)', 'kpoints']

# open file as csv which is tab separated, and ignores blank lines and comments
data = pd.read_csv('_data/mos2_kpoints.csv', 
                    sep=',',skip_blank_lines=True,header=0,names=["Energy Cutoff", "Total Energy","Wall Time","K Points"])

# conversion factor for eV
rydberg = 0.0734986176

# print dataframe of csv file
print(data)

total_energy = data["Total Energy"] * rydberg
sb.set_context("notebook", font_scale=1.0, rc={"lineswidth": 2.5})
sb.set_style("whitegrid")

# plot total energy data
fig = sb.scatterplot(data=data, x="K Points", y=total_energy).figure
#data.plot.scatter(x = 3, y = 1, figsize=(20,10),s=300)
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.6f'))

plt.title("Total Energy vs K Points")
plt.tight_layout()

plt.show()

fig.savefig("etot_vs_kpoints.png")

