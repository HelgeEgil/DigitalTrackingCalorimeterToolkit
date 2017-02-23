# Program to calculate the range straggling for different materials, for different ranges
# And check if the WEPL converted range straggling is similar to the range straggling for water


# 1) parameters

from math import *

def getRange(energy, material):
    return material.alpha * energy ** material.p

def getRangeWEPL(Rr, material, water):
    aw = water.alpha
    ar = material.alpha
    pw = water.p
    pr = material.p

    return aw * (Rr / ar) ** (pw / pr)

def getSigma(range, material):
    b = 3 - 2./material.p
    apr = material.alpha_prime
    pr = material.p
    ar = material.alpha

    a = apr * (pr ** 2 * ar ** (2./pr)) / b
    
    return sqrt(a * range ** b)

def getSigmaWEPL(Rr, material, water):
    sigma = getSigma(R, material)
    aw = water.alpha
    ar = material.alpha
    pw = water.p
    pr = material.p
    Rw = getRangeWEPL(Rr, material, water)

    sigmaWEPL = 2 * aw * ((Rr + sigma/2.) / ar) ** (pw / pr) - Rw
    
    return sigmaWEPL

class Param:
    def __init__(self, p, alpha, alpha_prime, name):
        self.p = p
        self.alpha = alpha
        self.alpha_prime = alpha_prime
        self.name = name

water = Param(1.725517, 0.027681, 0.008700, "Water")
aluminum = Param(1.707283, 0.014467, 0.020382, "Aluminum")
tungsten = Param(1.630000, 0.004320, 0.121442, "Tungsten")

energies = [0.1, 2, 10, 50, 100, 150, 200, 250]
materials = [water, aluminum, tungsten]

for energy in energies:
    print("\nFor {} MeV: ".format(energy))
    Rw = getRange(energy, water)
    sigmaWater = getSigma(Rw, water)

    for material in materials:
        R = getRange(energy, material)
        sigma = getSigma(R, material)
        sigmaWEPL = getSigmaWEPL(R, material, water)
        name = material.name
        per = sigma / R * 100
        ratio = getRangeWEPL(R, material, water) / R

        print("For {}:\t Range = {:6.2f} mm, sigma = {:5.4f} mm ({:.2f} %), sigmaWEPL = {:6.2f} mm, ratio = {:4.2f}".format(name, R, sigma, per, sigmaWEPL, ratio))



