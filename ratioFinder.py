import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd

dat1 = Table.read("abell57_2.emission", format="ascii.tab")
# print(dat1)
dat1.rename_column("O  2 3728.81A", "O2")
dat1.rename_column("Ne 3 3868.76A", "Ne3")
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

with open("abell57.in", "r") as f:
    lines = f.readlines()
    temp = lines[2].split(" ")[-2]
    density = lines[4].split("=")[-1]
    density = float(density.split(" ")[0])
    metal = float(lines[7].split("     ")[1])
    runParams = lines[:9]

with open("abell57Params.dat", "w") as f:
    for param in runParams:
        f.write(param)
    f.write("\n")
    f.write("="*50+"\n")

df = pd.DataFrame({'Species': colNames, 'Relative Flux': lineRatios*100})
df.to_csv('abell57Params.dat', sep=',', mode='a', float_format='%.5f', header=True, index=False)

'''

alphaBetaRatio = np.mean(dat1["HAlpha"]/dat1["HBeta"])
O3Ratio = np.mean((dat1["O3 5007"]+dat1["O3 4961"])/dat1["O3 4363"])
O3HBetaRatio = np.mean(dat1["O3 5007"]/dat1["HBeta"])
He5876HBetaRatio = np.mean(dat1["He1 5876"]/dat1["HBeta"])

with open("abell57.in", "r") as f:
    lines = f.readlines()
    temp = lines[2].split(" ")[-2]
    density=lines[4].split("=")[-1]

hebAlphaBetaRatio = 2.81
hebO3HBetaRatio = 1.311
hebHe5876HBetaRatio = 0.168

ratios = ["HAlpha", "O III 5007", "He I 5876"]
values = [alphaBetaRatio, O3Ratio, O3HBetaRatio, He5876HBetaRatio]

# values1 = [alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio]
# values2 = [hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio]

ratioDict = {
    'CLOUDY': (alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio),
    'HEB': (hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio)
}

x = np.arange(len(ratios))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

# fig, ax = plt.subplots(layout='constrained')
#
# for dataset, ratio in ratioDict.items():
#     offset = width * multiplier
#     rects = ax.bar(x + offset, ratio, width, label=dataset)
#     ax.bar_label(rects, padding=3)
#     multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel(r'Line Strength Relative to H$\beta$')
# ax.set_title('Emission Line Ratios')
# ax.text(2, 7, f'T={temp} LOG\n den = {density}', fontsize=12, ha='center', va='center')
# ax.set_xticks(x + width/2, ratios)
# ax.legend(loc='upper right')
# ax.set_ylim(0, 9)
# fig.savefig("ratiosPlot2.jpg", bbox_inches="tight", dpi=300)

dat2 = Table.read("abell57.spectrum", format="ascii.tab")
# a = dat2["line"][0]
print(dat2["Cont  nu"][5381])
wavelength = (1240/dat2["Cont  nu"])*10
lineNames = dat2['lineID']
lineNames = np.asarray(lineNames)
relFlux = dat2["total"]/(dat2["total"][5381])
ind = np.where(np.logical_and(relFlux > 0.05, wavelength < 9500))[0]
fig1, ax1 = plt.subplots()
ax1.scatter(wavelength[5381], relFlux[5381], color='red', s=5)
wavelength = wavelength[ind]
relFlux = relFlux[ind]
ax1.scatter(wavelength, relFlux, s=5)
for i, label in enumerate(lineNames[ind]):
    if relFlux[i] > 0.5 and relFlux[i] < 3.5:
        ax1.text(wavelength[i], relFlux[i]+0.1, label, fontsize=8, ha='center', va='bottom', rotation=90)
    if relFlux[i] > 3.5:
        ax1.text(wavelength[i], relFlux[i] - 0.5, label, fontsize=8, ha='center', va='bottom', rotation=90)
ax1.set_xlim(-0.01,10000)
ax1.set_xlabel(r"Wavelength [$\mathrm{\AA}$]")
ax1.set_ylabel(r"Relative Flux (F/F$_{\mathrm{H}_{\beta}}$)")
ax1.set_title("Spectrum of Abell 57")
# ax1.set_ylim(-1,1)
fig1.savefig("spectrum.jpg", bbox_inches="tight", dpi=300)

relFlux = np.array(relFlux, dtype=np.float64)

dat3 = np.rec.fromarrays([lineNames[ind], wavelength, relFlux])
np.savetxt("abell57Spectrum.dat", dat3, fmt="%s\t%.5f\t%.5f", delimiter="\t",
           header="LineID\tWavelength[A]\tRelative Flux")
'''
"""
Line              Wavelength    Flux          Obs Ratio    De-reddened Ratio
Halpha          6564.25        13.9             772.              281.0
Hbeta           4862.32         1.80            100.              100.0
Hgamma      4341.10         0.658:           36.6:             62.3
Hdelta          4102.32         0.262            14.6              31.9
Hepsilon        3970.51        0.210:           11.7:             28.9
Hzeta           3889.69         0.277:           15.4:             40.8
[O III] 4363    4363.88         0.966            53.7              89.2
[O III] 4959    4959.85         0.901            50.1              45.9
[O III] 5007    5007.69         2.69            149.              131.1
[O II]  3727    3727.               …                …                  ...
[Ne III] 3869   3869.49         2.36            131.              352.4
[Ne III] 3968   3967.93         0.525:           29.2:             72.3
He I 4471       4472.46:        0.0758::          4.2::             6.2
He I 5876       5876.62         0.602            33.4              16.8
[O I] 6300      6299.64:        0.0785::          4.4::             1.8
[N II] 6584     6584.93         0.0887::          4.9::             1.8
He I 6678       6679.38:        0.142::           7.9::             2.7
"""
"""
[O II]          3727.28:        0.170::           5.2::
[O II]          3730.04:        0.141::           4.3::
H10             3799.21:        0.137             4.2
H9              3836.50:        0.132:            4.0:
[Ne III]        3869.62         2.04             62.4
Hzeta           3889.80         0.393            12.0
[Ne III]        3968.36:        0.467:           14.3:
Hepsilon        3970.97:        0.319             9.8
Hdelta          4102.49         0.624            19.1
Hgamma          4341.26         1.16             35.5
[O III]         4363.96         0.316             9.7
He I            4472.15         0.110             3.4
He II           4686.79         0.879            26.9
[Ar IV]         4712.33         0.163:            5.0:
[Ar IV]         4742.15:        0.0930::          2.8::
Hbeta           4862.49         3.27            100
[O III]         4960.17        11.2             342
[O III]         5008.13        34.5            1055
He II           5412.61:        0.106:            3.2:
He I            5877.11         0.464            14.2
Halpha          6564.07        13.5             413
[N II]          6585.13         0.261             8.0
He I            6679.40:        0.128:            3.9:
[S II]          6717.75:        0.0479::          1.5::
"""
