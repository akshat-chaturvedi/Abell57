import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd
from scipy.interpolate import interp1d

dat1 = Table.read("Abell57_combined_weighted_spectrum.dat", format="ascii")  # Observed HET Spectrum
# Spacing for dat1 = 0.7A
dat2 = Table.read("dahot_hhe_d9060.flux", format="ascii")  # Theoretical Spectrum
dat3 = Table.read("a57_reddened_model_flux.txt", format="ascii")  # Theoretical Spectrum -> Reddened
# Spacing for dat3 is variable
dat4 = Table.read("dahot_hhe_d9060_opt.WFP", format="ascii")  # Theoretical Spectrum -> Widened
dat5 = Table.read("a57_reddened_synth_spec.txt", format="ascii")  # Theoretical Spectrum -> Widened and Reddened

for i in range(len(dat1["f_lam"])):
    if np.isnan(dat1["f_lam"][i]):
        dat1["f_lam"][i] = 0

ind = np.where(dat1["f_lam"] != 0)[0]

# Accounting for Doppler Effect - radial velocity of 100 km/s as per Klaus
shiftedLambda = dat5["lambda"]+(100000/3e8)*dat5["lambda"]
dat5["lambdaShift"] = shiftedLambda


# fig, ax = plt.subplots()
# ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind])
# # ax.plot(dat3["col1"], dat5["fLambda"]*35e-33)
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(3550, 7000)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# fig.savefig("a57_observed.jpg", dpi=300)
#
# fig, ax = plt.subplots()
# ax.plot(dat2["lambda/A"][ind], dat2["F(lambda)"][ind])
# # ax.plot(dat3["col1"], dat5["fLambda"])
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(850, 3550)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# # ax.set_ylim(0,10e16)
# fig.savefig("a57_theoretical.jpg", dpi=300)
#
# fig, ax = plt.subplots()
# ax.plot(dat3["col1"], dat5["fLambda"])
# # ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind]*(25e30))
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(850, 7000)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# # ax.set_xlim(4860,4870)
# fig.savefig("a57_reddened.jpg", dpi=300)

# fig, ax = plt.subplots()
# ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind], label="Observed")
# ax.plot(dat3["col1"], dat5["fLambda"]*40.75e-33, label="Reddened Theoretical")
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(6480, 6660)
# ax.set_ylim(1e-16, 0.3e-15)
# ax.legend(loc="upper left")
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# fig.savefig("a57_scaled.jpg", dpi=300)

fig, ax = plt.subplots()
ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind], label="Observed")
# ax.plot(dat5["lambdaShift"], dat5["fLambda"]*37.8e-33, label="Reddened Theoretical")
ax.plot(dat5["lambdaShift"], dat5["fLambda"]*33e-33, label="Theoretical")
# ax.scatter(dat5["lambdaShift"], dat5["fLambda"]*37.8e-33, label="Reddened Theoretical (points)", color="black")
# plt.plot(dat2["lambda/A"], dat2["F(nu)"])
ax.set_xlim(3700, 3750)
ax.set_ylim(0.4e-15, 0.6e-15)
ax.vlines(3727, 0.2e-15, 0.8e-15, color="black", linestyle="--")
ax.vlines(3730.5, 0.2e-15, 0.6e-15, color="red", linestyle="--")
ax.legend(loc="upper left")
ax.set_ylabel("Flux")
ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
fig.savefig("a57_scaled.jpg", dpi=300)

# O II 3726 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3723.5, dat5["lambdaShift"] < 3730.5))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3723.5, dat1["wavelength"] < 3730.5))[0]
interpWavelengths = np.arange(3723.5, 3730.5, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*33e-33
# print(len(obsFlux))
# print(len(modFlux))
OII3726FluxDiff = sum(abs(obsFlux-modFlux)*0.7)
# print(OII3726FluxDiff)

# O II 3728 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3727, dat5["lambdaShift"] < 3730.5))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3727, dat1["wavelength"] < 3730.5))[0]
interpWavelengths = np.arange(3727, 3730.5, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*33e-33
# print(len(obsFlux))
# print(len(modFlux))
OII3728FluxDiff = sum(abs(obsFlux-modFlux)*0.7)
# print(OII3726FluxDiff)

# H Alpha Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 6550, dat5["lambdaShift"] < 6580))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6550, dat1["wavelength"] < 6580))[0]
interpWavelengths = np.arange(6550, 6580, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*40e-33
# print(len(obsFlux))
# print(len(modFlux))
HAlphaFluxDiff = sum((obsFlux-modFlux)*0.7)
# print(HAlphaFluxDiff)

# H Beta Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4855, dat5["lambdaShift"] < 4870))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4855, dat1["wavelength"] < 4870))[0]
interpWavelengths = np.arange(4855, 4870, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35.5e-33
# print(len(obsFlux))
# print(len(modFlux))
HBetaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HBetaFluxDiff)

# H Gamma Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4337, dat5["lambdaShift"] < 4345))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4337, dat1["wavelength"] < 4345))[0]
interpWavelengths = np.arange(4337, 4345, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*34e-33
# print(len(obsFlux))
# print(len(modFlux))
HGammaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HGammaFluxDiff)
# print(sum(dat1["f_lam"][ind][lineIndex2]))
# print(sum(dat5["fLambda"][lineIndex]*31e-33))

# H Delta Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4100, dat5["lambdaShift"] < 4105))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4100, dat1["wavelength"] < 4105))[0]
interpWavelengths = np.arange(4100, 4105, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35e-33
HDeltaFluxDiff = sum((obsFlux-modFlux)*0.7)
# print(HDeltaFluxDiff)

# H Epsilon Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3970, dat5["lambdaShift"] < 3974))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3970, dat1["wavelength"] < 3974))[0]
interpWavelengths = np.arange(3970, 3974, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*34e-33
HEpsilonFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HEpsilonFluxDiff)

# H Zeta Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3886, dat5["lambdaShift"] < 3892))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3886, dat1["wavelength"] < 3892))[0]
interpWavelengths = np.arange(3886, 3892, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*33.5e-33
HZetaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HZetaFluxDiff)

'''
# He I 4471 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4469, dat5["lambdaShift"] < 4474))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4469, dat1["wavelength"] < 4474))[0]
interpWavelengths = np.arange(4469, 4474, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = sum(dat1["f_lam"][ind][lineIndex2]*0.7)
modFlux = sum(dat5["fLambda"][lineIndex]*34e-33*0.7)
HeI4471FluxDiff = abs(obsFlux - modFlux)
print(HeI4471FluxDiff)
'''

# He I 5876 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 5867, dat5["lambdaShift"] < 5882))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5867, dat1["wavelength"] < 5882))[0]
interpWavelengths = np.arange(5867, 5882, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
# print(len(obsFlux))
# print(len(modFlux))
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*38.5e-33
HeI5876FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HeI5876FluxDiff)

# He I 6678 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"]-6 > 6675, dat5["lambdaShift"]-6 < 6685))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6675, dat1["wavelength"] < 6685))[0]
interpWavelengths = np.arange(6675, 6685, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex]-6, dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*40.25e-33
HeI6678FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HeI6678FluxDiff)

# O III 4363 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4361, dat5["lambdaShift"] < 4370))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4361, dat1["wavelength"] < 4370))[0]
interpWavelengths = np.arange(4361, 4370, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*34e-33
OIII4363FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(OIII4363FluxDiff)

# Ne III 3868 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3865, dat5["lambdaShift"] < 3872))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3865, dat1["wavelength"] < 3872))[0]
interpWavelengths = np.arange(3865, 3872, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*33.25e-33
NeIII3868FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(NeIII3868FluxDiff)

# O III 4960 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4955, dat5["lambdaShift"] < 4966))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4955, dat1["wavelength"] < 4966))[0]
interpWavelengths = np.arange(4955, 4966, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35.5e-33
OIII4960FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(OIII4363FluxDiff)

# O III 5007 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 5003, dat5["lambdaShift"] < 5012))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5003, dat1["wavelength"] < 5012))[0]
interpWavelengths = np.arange(5003, 5012, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35.5e-33
OIII5007FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(OIII4363FluxDiff)

# Ne III 3968 Flux Difference - Blended with H Epsilon
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 3966, dat5["lambdaShift"] < 3970))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3966, dat1["wavelength"] < 3970))[0]
interpWavelengths = np.arange(3966, 3970, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*34e-33
NeIII3968FluxDiff = sum((obsFlux-modFlux)*0.7)

# N II 6583 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 6580, dat5["lambdaShift"] < 6589))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6580, dat1["wavelength"] < 6589))[0]
# breakpoint()
interpWavelengths = np.arange(6580, 6589, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*40e-33
NII6583FluxDiff = sum((obsFlux-modFlux)*0.7)

# He II 4686 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 4683, dat5["lambdaShift"] < 4690))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4683, dat1["wavelength"] < 4690))[0]
# breakpoint()
interpWavelengths = np.arange(4683, 4690, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35e-33
HeII4686FluxDiff = sum((obsFlux-modFlux)*0.7)

# He II 5411 Flux Difference
lineIndex = np.where(np.logical_and(dat5["lambdaShift"] > 5409, dat5["lambdaShift"] < 5420))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5409, dat1["wavelength"] < 5420))[0]
interpWavelengths = np.arange(5409, 5420, 0.7)
interpFunc = interp1d(dat5["lambdaShift"][lineIndex], dat5["fLambda"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*36.75e-33
HeII5411FluxDiff = sum((obsFlux-modFlux)*0.7)

speciesList = (["H Alpha", "H Beta", "H Gamma", "H Delta", "H Epsilon", "H Zeta", "He I 5876", "He I 6678",
                "O III 4363", "Ne III 3868", "O III 4960", "O III 5007", "N II 6583", "Ne III 3968", "He II 4686",
                "He II 5411", "O II 3726", "O II 3728"])

fluxArray = np.asarray([HAlphaFluxDiff, HBetaFluxDiff, HGammaFluxDiff, HDeltaFluxDiff, HEpsilonFluxDiff, HZetaFluxDiff,
                        HeI5876FluxDiff, HeI6678FluxDiff, OIII4363FluxDiff, NeIII3868FluxDiff, OIII4960FluxDiff,
                        OIII5007FluxDiff, NII6583FluxDiff, NeIII3968FluxDiff, HeII4686FluxDiff, HeII5411FluxDiff,
                        OII3726FluxDiff, OII3728FluxDiff])

relFluxArray = fluxArray/HBetaFluxDiff*100

relFluxArray = ["{:.2f}".format(x) for x in relFluxArray]
fluxArray = ["{:.3e}".format(x) for x in fluxArray]

commentArray = ["", "", "", "", "HEpsilon is blended with [Ne III] 3967", "", "Adjacent to Na D Line",
                "Unsure, had to blueshift model further", "", "",
                "", "", "Almost blended with H Alpha, includes some flux from that", "Blended with H Epsilon", "", "",
                "", ""]

df = pd.DataFrame({'Species': speciesList, 'Relative Flux': relFluxArray, 'Absolute Flux': fluxArray,
                   "Comments": commentArray})

df.to_csv('measuredFluxes.txt', sep='\t', mode='w', header=True, index=False)
