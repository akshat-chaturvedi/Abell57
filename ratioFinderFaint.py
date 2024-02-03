import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd
import plotly.graph_objects as go


# dat1 = Table.read("abell57CSFaint.emission", format="ascii.tab")
dat1 = Table.read("abell57FaintPN.emission", format="ascii.tab")
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
dat1.rename_column("He 1 3888.60A", "He1 3889")
dat1.rename_column("S  2 6730.81A", "S2 6731")

lineRatios = []
colNames = []
for colName in dat1.colnames:
    lineRatios.append(np.mean(dat1[colName][10:30]/dat1["HBeta"][10:30]))
    colNames.append(colName)

# print(colNames[1:])
# print(lineRatios[1:])
colNames = np.asarray(colNames[1:])
lineRatios = np.asarray(lineRatios[1:])
# fileRatio = np.asarray(fileRatio)

# breakpoint()

alphaBetaRatio = np.mean(dat1["HAlpha"][10:30]/dat1["HBeta"][10:30])
O3Ratio = np.mean((dat1["O3 5007"][10:30]+dat1["O3 4961"][10:30])/dat1["O3 4363"][10:30])
O3HBetaRatio = np.mean(dat1["O3 5007"][10:30]/dat1["HBeta"][10:30])
He5876HBetaRatio = np.mean(dat1["He1 5876"][10:30]/dat1["HBeta"][10:30])
Ne3869HBetaRatio = np.mean(dat1["Ne3 3868"][10:30] / dat1["HBeta"][10:30])
ArIV4711HBetaRatio = np.mean(dat1["Ar4 4711"][10:30]/dat1["HBeta"][10:30])
O2DoubletRatio = np.mean((dat1["O2 3726"][10:30]+dat1["O2 3729"][10:30])/dat1["HBeta"][10:30])
N2HBetaRatio = np.mean(dat1["N2 6583"][10:30]/dat1["HBeta"][10:30])
S2HBetaRatio = np.mean(dat1["S2 6716"][10:30]/dat1["HBeta"][10:30])
He4686HBetaRatio = np.mean(dat1["He2 4685"][10:30]/dat1["HBeta"][10:30])
S2DoubletHBetaRatio = np.mean((dat1["S2 6731"][10:30]+dat1["S2 6716"][10:30])/dat1["HBeta"][10:30])
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
# df.to_csv('abell57CSFaintParams.dat', sep=',', mode='a', float_format='%.5f', header=True, index=False)
df.to_csv('abell57FaintPNParams1.txt', sep='\t', mode='w', float_format='%.5f', header=True, index=False)

hebAlphaBetaRatio = 2.87*100
hebO3HBetaRatio = 9.84*10
hebHe5876HBetaRatio = 0.118*100
hebNe3869HBetaRatio = 0.872*100
hebArIV4711HBetaRatio = 0.055*100
hebO2DoubletRatio = 0.196*100
hebN2HBetaRatio = 0.047*100
hebS2HBetaRatio = 0.26*100
hebHe4686HBetaRatio = 0.052*100

# ratios = ["HAlpha", "O III 5007", "He I 5876", "Ne III 3869", "Ar IV 4711", "O2 Doublet", "N II 6584", "S II 6716"]
# values = [alphaBetaRatio, O3Ratio, O3HBetaRatio, He5876HBetaRatio, Ne3869HBetaRatio, ArIV4711HBetaRatio, O2DoubletRatio,
#           N2HBetaRatio, S2HBetaRatio]

# ratios1 = ["HAlpha", "[O III] 5007", "He I 5876", "[Ne III] 3869", "[Ar IV] 4711", "[S II] 6717", "He II 4686", "[O II] Doublet"]
ratios1 = ["[O II]\n Doublet", "[Ne III]\n 3869", "He II\n 4686", "[Ar IV]\n 4711", "[O III]\n 5007", "He I\n 5876", r"H$\alpha$"+"\n 6562",
           "[S II]\n Doublet"]
#values1 = [alphaBetaRatio, O3Ratio, O3HBetaRatio, He5876HBetaRatio, Ne3869HBetaRatio, ArIV4711HBetaRatio, S2HBetaRatio,
           #He4686HBetaRatio, O2DoubletRatio]

# values1 = [alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio]
# values2 = [hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio]

ratioDict = {
    # 'CLOUDY': (alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio, Ne3869HBetaRatio, ArIV4711HBetaRatio, O2DoubletRatio,
    #            N2HBetaRatio, S2HBetaRatio),
    # 'OBSERVED': (hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio, hebNe3869HBetaRatio, hebArIV4711HBetaRatio,
    #         hebO2DoubletRatio, hebN2HBetaRatio, hebS2HBetaRatio)
    # 'CLOUDY': (alphaBetaRatio, O3HBetaRatio, He5876HBetaRatio, Ne3869HBetaRatio, ArIV4711HBetaRatio,
    #            S2HBetaRatio, He4686HBetaRatio, O2DoubletRatio),
    # 'OBSERVED': (hebAlphaBetaRatio, hebO3HBetaRatio, hebHe5876HBetaRatio, hebNe3869HBetaRatio, hebArIV4711HBetaRatio,
    #              hebHe4686HBetaRatio, hebS2HBetaRatio, hebO2DoubletRatio)
    'CLOUDY': (float(f"{O2DoubletRatio*100:.1f}"), float(f"{Ne3869HBetaRatio*100:.1f}"), float(f"{He4686HBetaRatio*100:.2f}"),
               float(f"{ArIV4711HBetaRatio*100:.1f}"), float(f"{O3HBetaRatio*10:.1f}"), float(f"{He5876HBetaRatio*100:.1f}"),
               float(f"{alphaBetaRatio*100:.0f}"), float(f"{S2DoubletHBetaRatio*100:.1f}")),
    'OBSERVED': (hebO2DoubletRatio, hebNe3869HBetaRatio, hebS2HBetaRatio, hebArIV4711HBetaRatio, hebO3HBetaRatio,
                 hebHe5876HBetaRatio, hebAlphaBetaRatio, hebHe4686HBetaRatio)
    }

x = np.arange(len(ratios1))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

label_format = '{:g}'
# Add labels on top of the bars using ax.bar_label()
# ax.bar_label(bars, labels=[label_format(value) for value in data.values()])

#
# fig, ax = plt.subplots(layout='constrained')
# for dataset, ratio in ratioDict.items():
#     offset = width * multiplier
#     rects = ax.bar(x + offset, ratio, width, label=dataset)
#     ax.bar_label(rects, padding=5, rotation=90, fontsize="medium")
#     multiplier += 1
#
# ax.set_ylabel(r'Line Strength Relative to H$\beta$', fontsize="x-large")
# ax.set_title(r'Emission Line Ratios (H$\beta$=100)', fontsize="x-large")
# # ax.text(4.5, 9, f'T={temp} LOG\n den = {density:.5f} LOG\n Metallicity = {metal:.5f} LOG', fontsize=12, ha='center', va='center')
# ax.set_xticks(x + width/2, ratios1, fontsize="large")
# plt.yticks(fontsize=14)
# ax.legend(loc='upper left')
# ax.set_ylim(0, 320)
# fig.savefig("ratiosPlotFaintPN.jpg", bbox_inches="tight", dpi=300)



# print(np.mean(dat1["HAlpha"][10:30]))

# arr1 = []
# arr2 = []
# arr3 = []
# arr4 = []
# for i in range(1, 20):
#     arr1.append(abs(hebO3HBetaRatio-np.mean(dat1["O3 5007"][:i]/dat1["HBeta"][:i])))
#     arr2.append(abs(hebHe4686HBetaRatio-np.mean(dat1["He2 4685"][:i]/dat1["HBeta"][:i])))
#     arr3.append(abs(hebS2HBetaRatio-np.mean(dat1["S2 6716"][:i]/dat1["HBeta"][:i])))
#
# a = np.where(arr1==min(arr1))[0]
# b = np.where(arr2==min(arr2))[0]
# c = np.where(arr3==min(arr3))[0]

# print(a)
# print(b)
# print(c)

#print(np.mean((dat1["O2 3726"][10:30]/dat1["HBeta"][10:30])))
#print(np.mean((dat1["O2 3729"][10:30]/dat1["HBeta"][10:30])))

cloudy_values = ratioDict['CLOUDY']
cloudyVals = [19.3, 89.0, 5.1, 3.7, 1113, 12.8, 276, 4.0]
observed_values = ratioDict['OBSERVED']
obsVals = [19.6, 87.2, 26.0, 5.5, 984, 11.8, 287, 5.2]
# Bar chart creation
fig = go.Figure()

# Adding bars for CLOUDY
fig.add_trace(go.Bar(
    x=['[O II] Doublet', '[Ne III] 3689', 'He II 4686', '[Ar IV] 4711', '[O III] 5007', 'He I 5876', 'Hα', '[S II] Doublet'],
    y=cloudy_values,
    name='CLOUDY',
    text=[val for val in cloudyVals],  # Display values on top of bars
    textposition='outside',
    textfont=dict(size=30),  # Place text outside the bars
    textangle=270,  # Rotate text vertically
))

# Adding bars for OBSERVED
fig.add_trace(go.Bar(
    x=['[O II] Doublet', '[Ne III] 3689', 'He II 4686', '[Ar IV] 4711', '[O III] 5007', 'He I 5876', 'Hα', '[S II] Doublet'],
    y=observed_values,
    name='OBSERVED',
    text=[val for val in obsVals],  # Display values on top of bars
    textposition='outside',
    textfont=dict(size=30),  # Place text outside the bars
    textangle=270,  # Rotate text vertically
))

# Layout configuration
fig.update_layout(
    barmode='group',  # Display bars in groups
    title='Emission-Line Ratios (Hβ = 100)',
    title_font=dict(size=20),
    xaxis=dict(title="Species", tickfont=dict(size=15), title_font=dict(size=20)),
    yaxis=dict(title=r'Relative Line Strength', tickfont=dict(size=20), range=[0, 320], title_font=dict(size=20)),
    legend=dict(x=0.7, y=0.95, xanchor='right', yanchor='top', font=dict(size=15)), # Adjust legend position
    height= 500,  # Set the overall height of the plot
    width=600,  # Set the overall width of the plot
    # margin=dict(l=50, r=50, b=50, t=50)  # Adjust margin to leave space for labels
    paper_bgcolor='rgba(0,0,0,0)',  # Set background color to transparent
    plot_bgcolor='rgba(0,0,0,0)'  # Set plot area background color to transparent
    #config=dict(toImageButtonOptions=dict(format='jpg', scale=3))
)

# Show the plot
# fig.show()
#
cloudy_values = ratioDict['CLOUDY']
observed_values = ratioDict['OBSERVED']

# Elements
elements = ratios1

# Bar chart creation
fig1, ax1 = plt.subplots()

# Bar width
bar_width = 0.35

# Bar positions
cloudy_positions = np.arange(len(elements))
observed_positions = cloudy_positions + bar_width

# Plotting bars for CLOUDY
ax1.bar(cloudy_positions, cloudy_values, bar_width, label='CLOUDY')

# Plotting bars for OBSERVED
ax1.bar(observed_positions, observed_values, bar_width, label='OBSERVED')

# Adding text on top of bars
for i, value in enumerate(cloudyVals):
    if value < 500:
        ax1.text(cloudy_positions[i], value + 5, f'{value}', ha='center', va='bottom', fontsize=12, rotation=90)
    else:
        ax1.text(cloudy_positions[i], value/10 + 5, f'{value}', ha='center', va='bottom', fontsize=12, rotation=90)
for i, value in enumerate(obsVals):
    if value < 500:
        ax1.text(observed_positions[i], value + 5, f'{value}', ha='center', va='bottom', fontsize=12, rotation=90)
    else:
        ax1.text(observed_positions[i], value/10 + 5, f'{value}', ha='center', va='bottom', fontsize=12, rotation=90)
# Adding labels and title
ax1.set_xlabel('Species', fontsize=16)
ax1.set_ylabel('Line Strength Relative to Hβ = 100', fontsize=16)
ax1.set_title('Planetary Nebula Emission-Line Ratios', fontsize=18)
plt.tick_params(axis='y', which='major', labelsize=14)
plt.tick_params(axis='x', which='major', labelsize=11.5)

# Adding x-axis ticks and labels
ax1.set_xticks(cloudy_positions + bar_width / 2)
ax1.set_xticklabels(elements)

ax1.set_ylim(0, 320)
# Adding legend
ax1.legend(loc="upper center")

# Show the plot
plt.savefig("PNPLOTSTRIAL.pdf", bbox_inches="tight", dpi=300)