import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd

dat1 = Table.read("abell57_090.emission", format="ascii.tab")
# print(dat1)
dat1.rename_column("O  2 3728.81A", "O2")
dat1.rename_column("Ne 3 3868.76A", "Ne3 3869")
dat1.rename_column("O  3 4363.21A", "O3 4363")
dat1.rename_column("O  3 4960.91A", "O3 4961")
dat1.rename_column("O  3 5006.84A", "O3 5007")
dat1.rename_column("Blnd 5875.66A", "He1 5876")
dat1.rename_column("H  1 6562.80A", "HAlpha")
dat1.rename_column("H  1 4861.32A", "HBeta")
dat1.rename_column("h  1 4101.73A", "HDelta")
dat1.rename_column("h  1 4340.46A", "HGamma")
dat1.rename_column("He 1 6678.15A", "He1 6678")
dat1.rename_column("N  2 6583.45A", "N2 6583")
dat1.rename_column("Ne 3 3967.47A", "Ne3 3968")
dat1.rename_column("H  1 3889.02A", "HZeta")
dat1.rename_column("H  1 3970.08A", "HEpsilon")
dat1.rename_column("O  1 6300.30A", "O1 6300")

lineRatios = []
colNames = []
for colName in dat1.colnames:
    lineRatios.append(np.mean(dat1[colName][:100]/dat1["HBeta"][:100]))
    colNames.append(colName)

# print(colNames[1:])
# print(lineRatios[1:])
colNames = np.asarray(colNames[1:])
lineRatios = np.asarray(lineRatios[1:])

with open("abell57_ebv090.in", "r") as f:
    lines = f.readlines()
    temp = lines[2].split(" ")[-2]
    density = lines[4].split("=")[-1]
    density = float(density.split(" ")[0])
    metal = float(lines[8].split("     ")[1])
    runParams = lines[:9]

with open("abell57_ebv090_Params.dat", "w") as f:
    for param in runParams:
        f.write(param)
    f.write("\n")
    f.write("="*50+"\n")

df = pd.DataFrame({'Species': colNames, 'Relative Flux': lineRatios*100})
df.to_csv('abell57_ebv090_Params.dat', sep=',', mode='a', float_format='%.5f', header=True, index=False)

alphaBetaRatio = np.mean(dat1["HAlpha"]/dat1["HBeta"])
O3Ratio = np.mean((dat1["O3 5007"]+dat1["O3 4961"])/dat1["O3 4363"])
O3HBetaRatio = np.mean(dat1["O3 5007"]/dat1["HBeta"])
He5876HBetaRatio = np.mean(dat1["He1 5876"]/dat1["HBeta"])
Ne33869HBetaRatio = np.mean(dat1["Ne3 3869"]/dat1["HBeta"])
O34363HBetaRatio = np.mean(dat1["O3 4363"]/dat1["HBeta"])
O34960HBetaRatio = np.mean(dat1["O3 4961"]/dat1["HBeta"])
Ne33968HBetaRatio = np.mean(dat1["Ne3 3968"]/dat1["HBeta"])
gammaBetaRatio = np.mean(dat1["HGamma"]/dat1["HBeta"])
deltaBetaRatio = np.mean(dat1["HDelta"]/dat1["HBeta"])
He6678HBetaRatio = np.mean(dat1["He1 6678"]/dat1["HBeta"])
O16300HBetaRatio = np.mean(dat1["O1 6300"]/dat1["HBeta"])
N26584HBetaRatio = np.mean(dat1["N2 6583"]/dat1["HBeta"])


with open("abell57_ebv090.in", "r") as f:
    lines = f.readlines()
    temp = lines[2].split(" ")[-2]
    density = lines[4].split("=")[-1]

hebAlphaBetaRatio = 2.908
hebO3HBetaRatio = 1.095
hebHe5876HBetaRatio = 0.148
hebNe33869HBetaRatio = 2.651
hebO34363HBetaRatio = 0.724
hebO34960HBetaRatio = 0.398
hebNe33968HBetaRatio = 0.685
hebGammaBetaRatio = 0.451
hebDeltaBetaRatio = 0.237
hebHe6678HBetaRatio = 0.014
hebO16300HBetaRatio = 0.019
hebN26584HBetaRatio = 0.042

'''
H-a              6564.25     290.8
H-b              4862.32     100.0
H-g              4341.10      45.1
H-d              4102.32      23.7
[O III] 4363     4363.88      72.4
[O III] 4959     4959.85      39.8
[O III] 5007     5007.69     110.2
[Ne III] 3869    3869.49     265.1
He I 5876        5876.62      14.8
[N II] 6584      6584.93       4.2
He I 6678        6679.38       1.4
'''

ratios = ["HAlpha", "O III 5007", "He I 5876", "Ne III 3869", "O III 4363", "O III 4960", "Ne III 3968", "HGamma",
          "HDelta", "He 6678", "O I 6300", "N II 6584"]
values = [alphaBetaRatio, O3Ratio, O3HBetaRatio, He5876HBetaRatio, Ne33869HBetaRatio, O34363HBetaRatio, O34960HBetaRatio,
          Ne33968HBetaRatio, gammaBetaRatio, deltaBetaRatio, He6678HBetaRatio, O16300HBetaRatio, N26584HBetaRatio]

# values1 = [alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio]
# values2 = [hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio]

ratioDict = {
    'CLOUDY': (alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio, Ne33869HBetaRatio, O34363HBetaRatio, O34960HBetaRatio,
               Ne33968HBetaRatio, gammaBetaRatio, deltaBetaRatio, He6678HBetaRatio, O16300HBetaRatio, N26584HBetaRatio),
    'HEB': (hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio, hebNe33869HBetaRatio, hebO34363HBetaRatio,
            hebO34960HBetaRatio, hebNe33968HBetaRatio, hebGammaBetaRatio, hebDeltaBetaRatio, hebHe6678HBetaRatio,
            hebO16300HBetaRatio, hebN26584HBetaRatio)
}

x = np.arange(len(ratios))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for dataset, ratio in ratioDict.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, ratio, width, label=dataset)
    ax.bar_label(rects, padding=3, rotation=90, fontsize="small")
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r'Line Strength Relative to H$\beta$')
ax.set_title('CEK Emission Line Ratios $E(B-V)$ = 0.90')
ax.text(10, 2, f'T={temp} LOG\n den = {density}', fontsize=12, ha='center', va='center')
ax.set_xticks(x + width/2, ratios, rotation=90)
ax.legend(loc='upper right')
ax.set_ylim(0, 4)
fig.savefig("ratiosPlot_ebv090.jpg", bbox_inches="tight", dpi=300)