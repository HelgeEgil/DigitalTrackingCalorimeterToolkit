import numpy as np
import pandas as pd

headers = ["datapoints", "Q1BK", "Q2BK", "Q3BK", "Q1BKinv", "Q2BKinv", "Q3BKinv", 
        "Q1Ulmer", "Q2Ulmer", "Q3Ulmer", "Q1Ulmerinv", "Q2Ulmerinv", "Q3Ulmerinv", 
        "Q1Spline", "Q2Spline", "Q3Spline", "Q1Splineinv", "Q2Splineinv", "Q3Splineinv", 
        "Q1Linear", "Q2Linear", "Q3Linear", "Q1Linearinv", "Q2Linearinv", "Q3Linearinv"]

csvInput = pd.read_csv("../../DTCToolkit/OutputFiles/MedianValuesForParameterization.csv", sep=" ", names=headers)
uniqueDatapoints = set(csvInput["datapoints"].values)

pd_list = []
for datapoints in uniqueDatapoints:
    rowfilter = csvInput["datapoints"] == datapoints
    pd_list.append(csvInput[rowfilter].mean().to_frame().transpose())

outputdf = pd.concat(pd_list)
outputdf.to_csv("../../DTCToolkit/OutputFiles/MedianValuesForParameterization_merged.csv", sep=" ", index=False, header=False)
