import numpy as np
from matplotlib import pyplot as plt

with open('output/testSingles.dat', 'r') as ifile:
    line = ifile.readline().split()
    print line
    print len(line)
    timeid, eventid, sourceid, sourcex, sourcey, sourcez, \
        volumeid1, volumeid2, volumeid3, volumeid4, volumeid5, volumeid6, \
        timestamp, edep, singlex, singley, singlez, \
        ncomptonphantom, ncomptondetector, \
        nrayleighphantom, nrayleighdetector, \
        namecomptonphantom, namerayleighphantom \
         = line

    print "XYZ of event {}: ({}, {}, {})".format(eventid, singlex, singley, singlez)

