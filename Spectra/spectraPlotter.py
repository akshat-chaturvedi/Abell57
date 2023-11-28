import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd
from scipy.interpolate import interp1d

dat1 = Table.read("Abell57_combined_weighted_spectrum.dat", format="ascii")  # Observed HET Spectrum
# Spacing for dat1 = 0.7A
dat2 = Table.read("dahot_hhe_d9060.flux", format="ascii")  # Theoretical Spectrum
dat3 = Table.read("a57_reddened_model_flux.txt", format="ascii")  # Theoretical Spectrum -> De-reddened
# Spacing for dat3 is variable
dat4 = Table.read("dahot_hhe_d9060_opt.WFP", format="ascii")  # Theoretical Spectrum -> Widened


for i in range(len(dat1["f_lam"])):
    if np.isnan(dat1["f_lam"][i]):
        dat1["f_lam"][i] = 0

ind = np.where(dat1["f_lam"] != 0)[0]

# fig, ax = plt.subplots()
# ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind])
# # ax.plot(dat3["col1"], dat3["col2"]*35e-33)
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(3550, 7000)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# fig.savefig("a57_observed.jpg", dpi=300)
#
# fig, ax = plt.subplots()
# ax.plot(dat2["lambda/A"][ind], dat2["F(lambda)"][ind])
# # ax.plot(dat3["col1"], dat3["col2"])
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(850, 3550)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# # ax.set_ylim(0,10e16)
# fig.savefig("a57_theoretical.jpg", dpi=300)
#
# fig, ax = plt.subplots()
# ax.plot(dat3["col1"], dat3["col2"])
# # ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind]*(25e30))
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(850, 7000)
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# # ax.set_xlim(4860,4870)
# fig.savefig("a57_reddened.jpg", dpi=300)

# fig, ax = plt.subplots()
# ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind], label="Observed")
# ax.plot(dat3["col1"], dat3["col2"]*40.75e-33, label="Reddened Theoretical")
# # plt.plot(dat2["lambda/A"], dat2["F(nu)"])
# ax.set_xlim(6480, 6660)
# ax.set_ylim(1e-16, 0.3e-15)
# ax.legend(loc="upper left")
# ax.set_ylabel("Flux")
# ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
# fig.savefig("a57_scaled.jpg", dpi=300)

# H Alpha Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 6510, dat3["col1"] < 6610))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6510, dat1["wavelength"] < 6610))[0]
interpWavelengths = np.arange(6510, 6610, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)

obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*40.75e-33
# print(len(obsFlux))
# print(len(modFlux))
HAlphaFluxDiff = sum((obsFlux-modFlux)*0.7)
# print(HAlphaFluxDiff)

# H Beta Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4840, dat3["col1"] < 4890))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4840, dat1["wavelength"] < 4890))[0]
interpWavelengths = np.arange(4840, 4890, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*36.5e-33
# print(len(obsFlux))
# print(len(modFlux))
HBetaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HBetaFluxDiff)

# H Gamma Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4338, dat3["col1"] < 4347))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4338, dat1["wavelength"] < 4347))[0]
interpWavelengths = np.arange(4338, 4347, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35e-33
# print(len(obsFlux))
# print(len(modFlux))
HGammaFluxDiff = sum((obsFlux-modFlux)*0.7)
# print(HGammaFluxDiff)
# print(sum(dat1["f_lam"][ind][lineIndex2]))
# print(sum(dat3["col2"][lineIndex]*31e-33))

# H Delta Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4100, dat3["col1"] < 4105))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4100, dat1["wavelength"] < 4105))[0]
interpWavelengths = np.arange(4100, 4105, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*36.5e-33
HDeltaFluxDiff = sum((obsFlux-modFlux)*0.7)
# print(HDeltaFluxDiff)

# H Epsilon Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 3969, dat3["col1"] < 3974))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3969, dat1["wavelength"] < 3974))[0]
interpWavelengths = np.arange(3969, 3974, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35e-33
HEpsilonFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HEpsilonFluxDiff)

# H Zeta Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 3886, dat3["col1"] < 3892))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3886, dat1["wavelength"] < 3892))[0]
interpWavelengths = np.arange(3886, 3892, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35e-33
HZetaFluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HZetaFluxDiff)

'''
# He I 4471 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4469, dat3["col1"] < 4474))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4469, dat1["wavelength"] < 4474))[0]
interpWavelengths = np.arange(4469, 4474, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = sum(dat1["f_lam"][ind][lineIndex2]*0.7)
modFlux = sum(dat3["col2"][lineIndex]*34e-33*0.7)
HeI4471FluxDiff = abs(obsFlux - modFlux)
print(HeI4471FluxDiff)
'''

# He I 5876 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 5867, dat3["col1"] < 5882))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5867, dat1["wavelength"] < 5882))[0]
interpWavelengths = np.arange(5867, 5882, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
# print(len(obsFlux))
# print(len(modFlux))
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*39.5e-33
HeI5876FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HeI5876FluxDiff)

# He I 6678 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"]-5 > 6675, dat3["col1"]-5 < 6685))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6675, dat1["wavelength"] < 6685))[0]
interpWavelengths = np.arange(6675, 6685, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex]-5, dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*41.3e-33
HeI6678FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(HeI6678FluxDiff)

# O III 4363 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4360, dat3["col1"] < 4370))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4360, dat1["wavelength"] < 4370))[0]
interpWavelengths = np.arange(4360, 4370, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35.5e-33
OIII4363FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(OIII4363FluxDiff)

# Ne III 3868 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 3865, dat3["col1"] < 3875))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3865, dat1["wavelength"] < 3875))[0]
interpWavelengths = np.arange(3865, 3875, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*34.5e-33
NeIII3868FluxDiff = sum((obsFlux-modFlux[:-1])*0.7)
# print(NeIII3868FluxDiff)

# O III 4960 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4944, dat3["col1"] < 4969))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4944, dat1["wavelength"] < 4969))[0]
# breakpoint()
interpWavelengths = np.arange(4944, 4969, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*36.5e-33
OIII4960FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(OIII4363FluxDiff)

# O III 5007 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 5000, dat3["col1"] < 5015))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5000, dat1["wavelength"] < 5015))[0]
interpWavelengths = np.arange(5000, 5015, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*36.5e-33
OIII5007FluxDiff = sum((obsFlux-modFlux)*0.7)
# print(OIII4363FluxDiff)

# Ne III 3968 Flux Difference - Blended with H Epsilon
lineIndex = np.where(np.logical_and(dat3["col1"] > 3967, dat3["col1"] < 3975))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 3967, dat1["wavelength"] < 3975))[0]
interpWavelengths = np.arange(3967, 3975, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*35.5e-33
NeIII3968FluxDiff = sum((obsFlux-modFlux)*0.7)

# N II 6583 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 6574, dat3["col1"] < 6590))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 6574, dat1["wavelength"] < 6590))[0]
# breakpoint()
interpWavelengths = np.arange(6574, 6590, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*41e-33
NII6583FluxDiff = sum((obsFlux-modFlux)*0.7)

# He II 4686 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 4680, dat3["col1"] < 4693))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 4680, dat1["wavelength"] < 4693))[0]
# breakpoint()
interpWavelengths = np.arange(4680, 4693, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*37e-33
HeII4686FluxDiff = sum((obsFlux-modFlux)*0.7)

fig, ax = plt.subplots()
ax.plot(dat1["wavelength"][ind], dat1["f_lam"][ind], label="Observed")
ax.plot(dat3["col1"], dat3["col2"]*37.8e-33, label="Reddened Theoretical")
ax.plot(dat4["lambda"], dat4["fLambda"]*37.8e-33, label="Reddened Theoretical")
#ax.scatter(dat3["col1"], dat3["col2"]*37.8e-33, label="Reddened Theoretical (points)", color="black")
# plt.plot(dat2["lambda/A"], dat2["F(nu)"])
ax.set_xlim(5400, 5430)
ax.set_ylim(0.2e-15, 4e-15)
ax.legend(loc="upper left")
ax.set_ylabel("Flux")
ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
fig.savefig("a57_scaled.jpg", dpi=300)

# He II 5411 Flux Difference
lineIndex = np.where(np.logical_and(dat3["col1"] > 5412, dat3["col1"] < 5414))[0]
lineIndex2 = np.where(np.logical_and(dat1["wavelength"] > 5412, dat1["wavelength"] < 5414))[0]
interpWavelengths = np.arange(5412, 5414, 0.7)
interpFunc = interp1d(dat3["col1"][lineIndex], dat3["col2"][lineIndex], kind='linear', fill_value='extrapolate')
newFluxes = interpFunc(interpWavelengths)
obsFlux = dat1["f_lam"][ind][lineIndex2]
modFlux = newFluxes*37e-33
HeII5411FluxDiff = sum((obsFlux-modFlux)*0.7)

speciesList = (["H Alpha", "H Beta", "H Gamma", "H Delta", "H Epsilon", "H Zeta", "He I 5876", "He I 6678",
                "O III 4363", "Ne III 3868", "O III 4960", "O III 5007", "N II 6583", "Ne III 3968", "He II 4686",
                "He II 5411"])

fluxArray = np.asarray([HAlphaFluxDiff, HBetaFluxDiff, HGammaFluxDiff, HDeltaFluxDiff, HEpsilonFluxDiff, HZetaFluxDiff,
                        HeI5876FluxDiff, HeI6678FluxDiff, OIII4363FluxDiff, NeIII3868FluxDiff, OIII4960FluxDiff,
                        OIII5007FluxDiff, NII6583FluxDiff, NeIII3968FluxDiff, HeII4686FluxDiff, HeII5411FluxDiff])

relFluxArray = fluxArray/HBetaFluxDiff*100

relFluxArray = ["{:.2f}".format(x) for x in relFluxArray]
fluxArray = ["{:.3e}".format(x) for x in fluxArray]

commentArray = ["", "", "", "", "HEpsilon is blended with [Ne III] 3967", "", "Unknown absorption feature in model",
                "Unsure, had to shift model to the blueward", "", "",
                "Had to include feature at 4944 A due to interpolation requirements", "",
                "Almost blended with H Alpha, includes some flux from that", "Blended with H Epsilon", "",
                "Very narrow emission line"]

df = pd.DataFrame({'Species': speciesList, 'Relative Flux': relFluxArray, 'Absolute Flux': fluxArray,
                   "Comments": commentArray})

df.to_csv('measuredFluxes.txt', sep='\t', mode='w', header=True, index=False)
