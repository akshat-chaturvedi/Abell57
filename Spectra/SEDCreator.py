import numpy as np
import pandas
from astropy.table import Table

dat1 = Table.read("dahot_hhe_d9060.flux", format="ascii")  # Full SED
dat2 = dat1["lambda/A", "F(lambda)"]

dat2 = dat2[np.unique(dat2["lambda/A"], axis=0, return_index=True)[1]]

dat2.write("abell57Unique.sed", format="ascii.tab", delimiter="\t", overwrite=True)