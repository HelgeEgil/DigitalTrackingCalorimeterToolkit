# This pandas program is supposed to replace the remaining energy (row 6) in findManyRangesDegrader.csv with the corresponding value in EnergyAfterDegraderPSTAR.csv (row 2).
# The connecting value is the degrader thickness - row 1 in EAD and row 1 in fMRD.

import numpy as np
import pandas as pd

mrd = pd.read_csv("findManyRangesDegraderCarbon.csv", sep=" ", names=["waterdegrader", "absorber", "nomrange", "nomstrag", "inelastic", "remainingenergy", "energystrag"])
ead = pd.read_csv("EnergyAfterDegraderPSTAR.csv", sep=" ", names=["waterdegrader", "remainingenergy", "energystrag"])

for index, row in mrd.iterrows():
   absorber = row["absorber"]
   waterdegrader = row["waterdegrader"]
   
   good_row = ead["waterdegrader"] == waterdegrader
   energy_new = ead[good_row].values[0][1]

   selection1 = mrd["absorber"] == absorber
   selection2 = mrd["waterdegrader"] == waterdegrader

   mrd.loc[selection1 & selection2, "remainingenergy"] = energy_new

mrd.to_csv("findManyRangesDegraderCarbonModified.csv", sep=" ", header=False, index=False)
