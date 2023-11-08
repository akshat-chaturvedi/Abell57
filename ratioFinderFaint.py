import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd

dat1 = Table.read("abell57CSFaint.emission", format="ascii.tab")
dat1.rename_column("O  2 3728.81A", "O2 3728")
dat1.rename_column("Ne 3 3868.76A", "Ne3 3868")
dat1.rename_column("O  3 4363.21A", "O3 4363")
dat1.rename_column("O  3 4960.91A", "O3 4961")
dat1.rename_column("O  3 5006.84A", "O3 5007")
dat1.rename_column("Blnd 5875.66A", "He1 5876")
dat1.rename_column("H  1 6562.80A", "HAlpha")
dat1.rename_column("H  1 4861.32A", "HBeta")
dat1.rename_column("h  1 4101.73A", "HDelta")
dat1.rename_column("H  1 4340.47A", "HGamma")
dat1.rename_column("H  1 3889.02A", "HZeta")
dat1.rename_column("H  1 3970.08A", "HEpsilon")
dat1.rename_column("H  1 3798.94A", "H10")
dat1.rename_column("He 1 4471.68A", "He1 4471")
dat1.rename_column("He 2 4685.80A", "He2 4685")
dat1.rename_column("He 2 5411.52A", "He2 5411")
dat1.rename_column("N  2 6583.45A", "N2 6583")
dat1.rename_column("Ar 4 4711.35A", "Ar4 4711")
dat1.rename_column("O  2 3726.00A", "O2 3726")
dat1.rename_column("O  2 3729.00A", "O2 3729")
dat1.rename_column("S  2 6716.44A", "S2 6716")
dat1.rename_column("Ne 3 3967.47A", "Ne3 3968")
dat1.rename_column("Ar 4 4740.12A", "Ar4 4740")
dat1.rename_column("He 1 6678.15A", "He1 6678")
dat1.rename_column("H  1 3835.39A", "HEta")

lineRatios = []
colNames = []
for colName in dat1.colnames:
    lineRatios.append(np.mean(dat1[colName][:100]/dat1["HBeta"][:100]))
    colNames.append(colName)

# print(colNames[1:])
# print(lineRatios[1:])
colNames = np.asarray(colNames[1:])
lineRatios = np.asarray(lineRatios[1:])
# fileRatio = np.asarray(fileRatio)

# breakpoint()

alphaBetaRatio = np.mean(dat1["HAlpha"][:100]/dat1["HBeta"][:100])
O3Ratio = np.mean((dat1["O3 5007"][:100]+dat1["O3 4961"][:100])/dat1["O3 4363"][:100])
O3HBetaRatio = np.mean(dat1["O3 5007"][:100]/dat1["HBeta"][:100])
He5876HBetaRatio = np.mean(dat1["He1 5876"][:100]/dat1["HBeta"][:100])
Ne3689HBetaRatio = np.mean(dat1["Ne3 3868"][:100]/dat1["HBeta"][:100])
ArIV4711HBetaRatio = np.mean(dat1["Ar4 4711"][:100]/dat1["HBeta"][:100])
O2DoubletRatio = np.mean((dat1["O2 3726"][:100]+dat1["O2 3729"][:100])/dat1["HBeta"][:100])
N2HBetaRatio = np.mean(dat1["N2 6583"][:100]/dat1["HBeta"][:100])
S2HBetaRatio = np.mean(dat1["S2 6716"][:100]/dat1["HBeta"][:100])
# breakpoint()

with open("abell57CSFaint.in", "r") as f:
    lines = f.readlines()
    temp = lines[2].split(" ")[-2]
    density = lines[4].split("=")[-1]
    density = float(density.split(" ")[0])
    metal = float(lines[7].split("     ")[1])
    runParams = lines[:9]

with open("abell57CSFaintParams.dat", "w") as f:
    for param in runParams:
        f.write(param)
    f.write("\n")
    f.write("="*50+"\n")

df = pd.DataFrame({'Species': colNames, 'Relative Flux': lineRatios*100})
df.to_csv('abell57CSFaintParams.dat', sep=',', mode='a', float_format='%.5f', header=True, index=False)

'''
hebAlphaBetaRatio = 2.906
hebO3HBetaRatio = 10.09
hebHe5876HBetaRatio = 0.112
hebNe3689HBetaRatio = 0.881
hebArIV4711HBetaRatio = 0.053
hebO2DoubletRatio = 2
hebN2HBetaRatio = 0.056
hebS2HBetaRatio = 0.01

ratios = ["HAlpha", "O III 5007", "He I 5876", "Ne III 3689", "Ar IV 4711", "O2 Doublet", "N II 6584", "S II 6717"]
values = [alphaBetaRatio, O3Ratio, O3HBetaRatio, He5876HBetaRatio, Ne3689HBetaRatio, ArIV4711HBetaRatio, O2DoubletRatio,
          N2HBetaRatio, S2HBetaRatio]

# values1 = [alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio]
# values2 = [hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio]

ratioDict = {
    'CLOUDY': (alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio, Ne3689HBetaRatio, ArIV4711HBetaRatio, O2DoubletRatio,
               N2HBetaRatio, S2HBetaRatio),
    'HEB': (hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio, hebNe3689HBetaRatio, hebArIV4711HBetaRatio,
            hebO2DoubletRatio, hebN2HBetaRatio, hebS2HBetaRatio)
}

x = np.arange(len(ratios))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0


fig, ax = plt.subplots(layout='constrained')
for dataset, ratio in ratioDict.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, ratio, width, label=dataset)
    ax.bar_label(rects, padding=5, rotation=90)
    multiplier += 1

ax.set_ylabel(r'Line Strength Relative to H$\beta$')
ax.set_title('Emission Line Ratios')
ax.text(4.5, 9, f'T={temp} LOG\n den = {density:.5f} LOG\n Metallicity = {metal:.5f} LOG', fontsize=12, ha='center', va='center')
ax.set_xticks(x + width/2, ratios, fontsize="small")
ax.legend(loc='upper right')
ax.set_ylim(0, 12)
fig.savefig("ratiosPlotFaint1.jpg", bbox_inches="tight", dpi=300)

'''