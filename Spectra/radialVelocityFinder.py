import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd
from scipy.interpolate import interp1d
from scipy import signal

dat1 = Table.read("Abell57_combined_weighted_spectrum.dat", format="ascii")  # Observed HET Spectrum
# Spacing for dat1 = 0.7A
dat2 = Table.read("dahot_hhe_d9060.flux", format="ascii")  # Theoretical Spectrum
dat3 = Table.read("a57_reddened_model_flux.txt", format="ascii")  # Theoretical Spectrum -> Reddened
# Spacing for dat3 is variable
dat4 = Table.read("dahot_hhe_d9060_opt.WFP", format="ascii")  # Theoretical Spectrum -> Widened
dat5 = Table.read("a57_reddened_synth_spec.txt", format="ascii")  # Theoretical Spectrum -> Widened and Reddened
dat6 = Table.read("a57_flattened.dat", format="ascii")  # Observed Spectrum -> Flattened

for i in range(len(dat1["f_lam"])):
    if np.isnan(dat1["f_lam"][i]):
        dat1["f_lam"][i] = 0

ind = np.where(dat1["f_lam"] != 0)[0]
ind1 = np.where(dat6["col2"] < 1.025)[0]  # Filtering out the emission lines
ind2 = np.where(np.logical_and(dat6["col1"] > 4000, dat6["col1"] < 4200))[0]

shiftedLambda = dat5["lambda"]+(100000/3e8)*dat5["lambda"]
dat5["lambdaShift"] = shiftedLambda

# Adding a starting radial velocity shift
initRadVel = 103000
radVelLambda = dat4["lambda"]+(initRadVel/3e8)*dat4["lambda"]
dat4["lambdaShift"] = radVelLambda

# interpWavelengths = np.arange(3650.0, 6949.8, 0.7)
# interpFunc = interp1d(dat4["lambdaShift"], dat4["relFLambda"], kind='linear', fill_value='extrapolate')
# newFluxes = interpFunc(interpWavelengths)


fig, ax = plt.subplots()

# ax.plot(dat6["col1"][ind1], dat6["col2"][ind1], label="Observed")
ax.plot(dat6["col1"], dat6["col2"], label="Observed")
ax.plot(dat4["lambda"], dat4["relFLambda"], label="Theoretical")
ax.plot(dat4["lambdaShift"], dat4["relFLambda"], label="Theoretical Shifted", alpha=0.8)
ax.set_xlim(4250, 4450)
ax.set_ylim(0.7, 1.2)
ax.legend(loc="upper left")
ax.set_ylabel("Flux")
ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
fig.savefig("a57_rad_vel.jpg", dpi=300)

for i in range(len(dat6["col2"])):
    if dat6["col2"][i] > 1.025:
        dat6["col2"][i] = 0
    else:
        pass

corArray = []
RadVelArray = []
for RadVel in range(-5000000, 5000000, 1000):
    radVelLambda = dat4["lambda"] + (RadVel / 3e8) * dat4["lambda"]
    dat4["lambdaShift"] = radVelLambda
    interpWavelengths = np.arange(4000, 4200, 0.7)
    interpFunc = interp1d(dat4["lambdaShift"], dat4["relFLambda"], kind='linear', fill_value='extrapolate')
    newFluxes = interpFunc(interpWavelengths)
    corVel = []
    for i in range(len(dat6["col2"][ind2])):
        corVel.append(dat6["col2"][ind2][i]*newFluxes[i])
        # print(corVel)
    corVelSum = sum(corVel)
    corArray.append(corVelSum)
    RadVelArray.append(RadVel)

maxCorVelInd = np.where(corArray == max(corArray))[0]
print(RadVelArray[maxCorVelInd[0]])
maxRadVel = RadVelArray[maxCorVelInd[0]]/1000
maxCorVelVal = max(corArray)

RadVelArray = np.asarray(RadVelArray)

fig, ax = plt.subplots()

ax.plot(RadVelArray/1000, corArray, label="Correlation")
ax.set_ylabel("Correlation")
ax.set_xlabel("Radial Velocity [km/s]")
ax.vlines(13.8, 220, 240, color="k", linestyle="--", label="Best Fit Radial Velocity")
# ax.text(maxRadVel-20, maxCorVelVal+0.01, f"{maxRadVel:.2f}")
ax.set_ylim(234.5, 235.5)
ax.set_xlim(-500, 500)
ax.legend(loc="upper right", fontsize="small")
fig.savefig("CorrelationFunction1.jpg", bbox_inches='tight', dpi=300)

# # H Gamma Flux Difference
# lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4337, dat5["lambdaShift"] < 4345))[0]
# lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4337, dat1["wavelength"] < 4345))[0]
# interpWavelengths = np.arange(4337, 4345, 0.7)
# interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
# newFluxes = interpFunc(interpWavelengths)
# obsFlux = dat1["f_lam"][ind][lineIndex2]
# modFlux = newFluxes*34e-33
# # print(len(obsFlux))
# # print(len(modFlux))
# HGammaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# # print(HGammaFluxDiff)
# # print(sum(dat1["f_lam"][ind][lineIndex2]))
# # print(sum(dat5["fLambda"][lineIndex]*31e-33))
#
# # H Delta Flux Difference
# lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4100, dat5["lambdaShift"] < 4105))[0]
# lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4100, dat1["wavelength"] < 4105))[0]
# interpWavelengths = np.arange(4100, 4105, 0.7)
# interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
# newFluxes = interpFunc(interpWavelengths)
# obsFlux = dat1["f_lam"][ind][lineIndex2]
# modFlux = newFluxes*35e-33
# HDeltaFluxDiff = sum((obsFlux-modFlux)*0.7)
