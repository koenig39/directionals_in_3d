import lasio
import pandas as pd
import matplotlib.pyplot as plt

filename = "29757-CBL.las"
las = lasio.read(filename)
plt.figure(figsize=(16, 8))
for item in las.curves:
    if item.mnemonic == "DEPT":
        depthUnit = item.unit
        continue
    plt.plot(las[item.mnemonic],las["DEPT"],label=f"{item.descr} / {item.unit}")

# plt.xscale('log')        

plt.xlabel('Measured Value')
plt.ylabel(f'Depth (MD) {depthUnit}')
plt.title(f'LAS log: {filename}')

# # Adding grid
plt.grid(True)
plt.gca().invert_yaxis()
plt.legend()
plt.show()
    
