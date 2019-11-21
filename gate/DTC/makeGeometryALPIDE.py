import matplotlib
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Tkinter import *
import os, ttk
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler

nm = 0.000001
um = 0.001
mm = 1.
cm = 10.
m  = 1000.
kYellow = "yellow"
kGrey = "grey"
kBlack = "black"
kRed = "red"
kGreen = "green"
kBlue = "blue"
kNoColor = False
kAluminium = "Aluminium"
kTungsten = "Tungsten"
kWater = "Water"
kScintillator = "myScintillator"
kFR4 = "FR4"
kKapton = "Kapton"
kEpoxy = "Epoxy"
kKapton = "Kapton"
kGlue = "Glue"
kChip = "Silicon"
kAir = "Air"

kNoMaterial = False
kNoMother = False
kVisible = True
kNotVisible = False
kIsFirstAbsorber = True

class Pos:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
        self.type = "Pos"

    def isPos(self):
        return x+y+z and True or False

    def __str__(self):
        return "({},{},{})".format(self.x, self.y, self.z)

    def __add__(self,other):
        
        if other.type == "Pos":
            x = self.x + other.x
            y = self.y + other.y
            z = self.z + other.z

        elif other.type == "Size":
            x = self.x + other.dx
            y = self.y + other.dy
            z = self.z + other.dz

        elif other.type == "Box":
            x = self.x + other.pos.x
            y = self.y + other.pos.y
            z = self.z + other.pos.z

        return Pos(x,y,z)

    def __sub__(self, other):
        if other.type == "Pos":
            x = self.x - other.x
            y = self.y - other.y
            z = self.z - other.z

        elif other.type == "Size":
            x = self.x - other.dx
            y = self.y - other.dy
            z = self.z - other.dz
            
        return Pos(x,y,z)

class Size:
    def __init__(self, dx = 0, dy = 0, dz = 0):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.type = "Size"

    def __str__(self):
        return "({},{},{})".format(self.dx, self.dy, self.dz)

    def __add__(self, other):
        dx = self.dx + other.dx
        dy = self.dy + other.dy
        dz = self.dz + other.dz
        return Size(dx, dy, dz)

    def __sub__(self, other):
        dx = self.dx - other.dx
        dy = self.dy - other.dy
        dz = self.dz - other.dz
        return Size(dx, dy, dz)

    def __div__(self, f):
        dx = self.dx / f
        dy = self.dy / f
        dz = self.dz / f
        return Size(dx, dy, dz)

    def __mul__(self, f):
        dx = self.dx * f
        dy = self.dy * f
        dz = self.dz * f
        return Size(dx, dy, dz)

class Rotation:
    def __init__(self, name, filename, rotation, translationA, translationB):
        self.name = name
        self.filename = filename
        self.rotation = rotation
        self.type = "Rotation"
        self.rotationmatrix = "1 0 0"
        self.translationA = translationA
        self.translationB = translationB

    def makeString(self):
        string = "Time s\n"
        string += "Rotation deg\n"
        string += "Translation mm\n"
        string += "{s} {deg} {rotationmatrix} {translation}\n".format(s="0", deg=0, rotationmatrix=self.rotationmatrix, translation=self.translationA)
        string += "{s} {deg} {rotationmatrix} {translation}\n".format(s="0", deg=self.rotation, rotationmatrix=self.rotationmatrix, translation=self.translationB)
        string += "\n"

        return string

    def printBox(self):
        print(self.makeString())

    def writeBox(self):
        f = open(self.filename, "w")
        f.write(self.makeString())
        f.close()

class Copy:
    def __init__(self, name, nlayers, translation):
        self.type = "Copy"
        self.name = name
        self.nlayers = nlayers
        self.translation = translation

    def makeString(self):
        string =  "/gate/{}/repeaters/insert linear\n".format(self.name)
        string +=  "/gate/{}/linear/autoCenter false\n".format(self.name) # this addition of line cost me 5 hours of debugging. DON'T REMOVE
        string += "/gate/{}/linear/setRepeatNumber {}\n".format(self.name, self.nlayers)
        string += "/gate/{}/linear/setRepeatVector {} mm\n".format(self.name, self.translation)

        string += "\n"

        return string

    def printBox(self):
        print(self.makeString())

    def writeBox(self):
        outputFile = open("Module.mac", "a")
        outputFile.write(self.makeString())
        outputFile.close()

class Box:
    def __init__(self, name = False, mother = False, pos = False, size = False, material = kNoMaterial, visible = kNotVisible, color = kNoColor, genericRepeater = False):
        self.name = name
        self.mother = mother
        self.pos = pos
        self.size = size
        self.material = material
        self.visible = visible
        self.type = "Box"
        self.color = color
        self.genericRepeater = genericRepeater

    def makeString(self):
        string = ""
        if self.mother:
            string += "/gate/{}/daughters/name {}\n".format(self.mother, self.name)

            if "scanner" in self.name or "phantom" in self.name:
                string += "/gate/{}/daughters/systemType scanner\n".format(self.mother)

            string += "/gate/{}/daughters/insert box\n".format(self.mother)
        
        string += "/gate/{}/geometry/setXLength {} mm\n".format(self.name, self.size.dx)
        string += "/gate/{}/geometry/setYLength {} mm\n".format(self.name, self.size.dy)
        string += "/gate/{}/geometry/setZLength {} mm\n".format(self.name, self.size.dz)

        if self.pos.isPos:
            string += "/gate/{}/placement/setTranslation {} {} {} mm\n".format(self.name, self.pos.x, self.pos.y, self.pos.z)
        if self.material:
            string += "/gate/{}/setMaterial {}\n".format(self.name, self.material)
        if self.color:
            string += "/gate/{}/vis/setColor {}\n".format(self.name, self.color)
        if not self.visible:
            string += "/gate/{}/vis/setVisible false\n".format(self.name)

        if self.genericRepeater:
            string += "\n/gate/{}/repeaters/insert genericRepeater\n".format(self.name)
            string += "/gate/{}/genericRepeater/setPlacementsFilename Module.placements\n".format(self.name)

        string += "\n"

        return string

    def printBox(self):
        print(self.makeString())

    def writeBox(self):
        outputFile = open("Module.mac", "a")
        outputFile.write(self.makeString())
        outputFile.close()

    def move(self, pos):
        self.pos += pos

    def addRotation(self):
        self.genericRepeater = True

    def rotatePosition(self, mat, deg, gapY, minY, maxY, minZ, maxZ, i=0, cop=0, fromZ=0):
        if mat == "1 0 0" and deg == 180:
            meanZ = maxZ
            meanY = (maxY + minY)/2.
            delta = maxZ - minZ

            if cop:
                z += (i+1) * float(cop.translation.split()[2])

            newPos = Pos(self.pos.x, 2*meanY - self.pos.y - gapY, 2*meanZ - self.pos.z)

            return newPos

        else:
            return self.pos

    def Contains(self, other):
        x1 = self.pos.x - self.size.dx/2.
        x2 = self.pos.x + self.size.dy/2.
        y1 = self.pos.y - self.size.dy/2.
        y2 = self.pos.y + self.size.dy/2.

        x1p = other.pos.x - other.size.dx/2.
        x2p = other.pos.x + other.size.dy/2.
        y1p = other.pos.y - other.size.dy/2.
        y2p = other.pos.y + other.size.dy/2.

        return x1 < x2p and x2 > x1p and y1 < y2p and y2 > y1p

    def __str__(self):
        x1 = self.pos.x - self.size.dx/2.
        x2 = self.pos.x + self.size.dy/2.
        y1 = self.pos.y - self.size.dy/2.
        y2 = self.pos.y + self.size.dy/2.
        z1 = self.pos.z - self.size.dz/2.
        z2 = self.pos.z + self.size.dz/2.

        return "BOX {} ({} -> {}, {} -> {}, {} -> {})".format(self.name, x1, x2, y1, y2, z1, z2)

class Module:
#    def __init__(self, name, mother, pos, absorberMaterial, nlayers = 40, firstAbsorberMaterial = kNoMaterial, absorberThickness = 1.5):
    def __init__(self, name, mother, pos, absorberMaterial, nlayers = 40, absorberThickness = 1.5):
        self.objects = []
        self.nlayers = 40
        self.pos = pos
        self.name = name
        self.chipSize = Size(90 * mm, 30 * mm)
        self.chipGap = Size(0, 0)
        self.chipPos = Pos(0, 0)

        self.absorberThickness = absorberThickness
        self.absorberMaterial = absorberMaterial
#        self.firstAbsorberThickness = absorberThickness
        self.absorberSize = Size(self.chipSize.dx, self.chipSize.dy, self.absorberThickness)
        self.glueThickness = 5 * um
        self.kaptonThickness = 75 * um
        self.spacerThickness = 260 * um
        self.AlCableThickness = 150 * um
        self.activeChipThickness = 50 * um
        self.absorberAirGapThickness = 1480 * um
        self.topFDI_AlThickness = 30 * um
        self.topFDI_PiThickness = 20 * um
        self.bottomFDI_AlThickness = 100 * um
        self.bottomFDI_PiThickness = 20 * um
        self.fromZ = 0
#         self.thickness = 2*self.glueThickness + self.kaptonThickness + self.AlCableThickness + self.activeChipThickness + self.absorberThickness + self.absorberAirGapThickness
#         self.moduleThickness = 2 * self.glueThickness + self.kaptonThickness + self.AlCableThickness + self.activeChipThickness
        
        self.thickness = 5 * self.glueThickness + self.absorberThickness + self.spacerThickness + self.activeChipThickness + self.absorberAirGapThickness + 
                        self.topFDI_AlThickness + self.topFDI_PiThickness + self.kaptonThickness + self.bottomFDI_AlThickness + self.bottomFDI_PiThickness

        self.thickness = 5 * self.glueThickness + self.absorberThickness + self.spacerThickness + self.activeChipThickness + self.absorberAirGapThickness + 
                        self.topFDI_AlThickness + self.topFDI_PiThickness + self.kaptonThickness + self.bottomFDI_AlThickness + self.bottomFDI_PiThickness

        self.size = Size(10*cm, 4*cm, self.thickness)
        self.mother = mother
#        self.firstAbsorberMaterial = firstAbsorberMaterial

    def getLeftSide(self):
        return self.pos - self.size / 2.

    def getRightSide(self):
        return self.pos + self.size / 2.

    def getLeftSideNoAbsorber(self):
        for obj in self.objects:
            if not obj.name.endswith("Layer") and obj.name != self.name and not "Absorber" in obj.name:
                return obj.pos.z - obj.size.dz / 2

    def getFirstObjectLeftSide(self):
        return self.object[0].pos.z - self.object[0].size.dz/2.

    def getLastObjectCenter(self):
        sumPoint = self.getLeftSide()
        for obj in self.objects:
            sumPoint += obj

        return sumPoint

    def getLastObjectRightSide(self):
        rightSide = -1000
        for obj in self.objects:
            if obj.name != self.name and not obj.name.endswith("Layer") and obj.type == "Box" and not "Absorber" in obj.name:
                rightSide = max(rightSide, obj.pos.z + obj.size.dz/2)

        return rightSide

    def getLayerThickness(self):
        rightSide = -1000
        leftSide = 1000
        for obj in self.objects:
            if obj.name != self.name and not "Layer" in obj.name and obj.type == "Box":
                rightSide = max(rightSide, obj.pos.z + obj.size.dz/2.)
                leftSide = min(leftSide, obj.pos.z - obj.size.dz/2.)

        return rightSide - leftSide

    def getThisObjectCenter(self, box):
        """Add objects from left to the right IFF the object overlaps in x,y with this object"""

        sumPoint = self.getLeftSide()
        sumPoint.z += box.size.dz / 2.

        print("Searching for overlaps with {}.".format(box))
        for obj in self.objects:
            if obj.type != "Box": continue
            if obj.Contains(box) and obj.name != self.name and not "Layer" in obj.name:
                sumPoint.z = obj.pos.z + obj.size.dz/2. + box.size.dz/2.
                print("{} overlaps with {}! New z is {}. ".format(box, obj, sumPoint.z))

        z = sumPoint.z
        print("From original position {}, new z is {}.".format(box, z))
        return z
    
    def getPosDepth(self, x, y, size, name):
        pos = Pos(x, y)
        xyBox = Box(name, kNoMother, pos, size)
        pos.z = self.getThisObjectCenter(xyBox)
        return pos

    def addObject(self, name, pos, size, material, visible, color):
        self.objects.append(Box(name, self.name, pos, size, material, visible, color))
        print("\033[1mAdded {} with mean pos {}\033[0m".format(self.objects[-1], self.objects[-1].pos))
   
    def addAbsorber(self, name, pos, size, material, visible, color):
        self.objects.append(Box(name, "Layer", pos, size, material, visible, color))
        print("\033[1mAdded {} with mean pos {}\033[0m".format(self.objects[-1], self.objects[-1].pos))

    def makeMotherModule(self):
        size = copy.copy(self.size)
        size.dz = self.moduleThickness
        pos = self.getPosDepth(0, 0, size, self.name)
        self.objects.append(Box(self.name, "Layer", pos, size))
        print("Added {}".format(self.objects[-1]))

    def makeMotherLayer(self):
        self.objects.append(Box("Layer", self.mother, self.pos, self.size))
        print("Added {}".format(self.objects[-1]))

    def makeAbsorber(self, n, isFirstAbsorber=False):
        name = "Absorber" + str(n)
        size = self.chipSize
        size.dz = self.absorberThickness
        pos = self.getPosDepth(0, 0, size, name)
        # pos.z += self.getLastObjectRightSide() - self.getLeftSideNoAbsorber()

        material = self.absorberMaterial
#        if self.firstAbsorberMaterial and isFirstAbsorber:
#            material = self.firstAbsorberMaterial
#            size.dz = self.firstAbsorberThickness
        
        self.addAbsorber(name, pos, size, material, kVisible, kGrey)

    def makeAbsorberAirGap(self):
        name = "AbsorberAirGap"
        size = Size(self.chipSize.dx, self.chipSize.dy, self.absorberAirGapThickness)
        pos = Pos(0, 0, self.objects[-1].pos.z + self.objects[-1].size.dz/2 + size.dz/2)

        self.objects.append(Box(name, "Layer", pos, size, "Air", kNotVisible, kNoColor))
        print("\033[1mAdded {} with mean pos {}\033[0m".format(self.objects[-1], self.objects[-1].pos))

    def makeKaptonGlue(self):
        name = "KaptonGlue"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kEpoxy, kVisible, kYellow)

    def makeKapton(self):
        name = "Kapton"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.kaptonThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kKapton, kVisible, kGreen)

    def makeAbsorberGlue(self):
        name = "AbsorberGlue"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kEpoxy, kVisible, kYellow)
   
    def makeSpacer(self):
        name = "Spacer"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.spacerThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kAluminium, kVisible, kYellow)

    def makeSpacerGlue(self):
        name = "SpacerGlue"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kEpoxy, kVisible, kYellow)

    def makeFDIGlue(self):
        name = "FDIGlue"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kEpoxy, kVisible, kYellow)

    def makeTopFDI_Al(self):
        name = "TopFDI_Al"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.topFDI_AlThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kAluminium, kVisible, kYellow)

    def makeTopFDI_Pi(self):
        name = "TopFDI_Pi"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.topFDI_PiThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kKapton, kVisible, kYellow)

    def makeBottomFDI_Al(self):
        name = "BottomFDI_Al"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.bottomFDI_AlThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kAluminium, kVisible, kYellow)

    def makeBottomFDI_Pi(self):
        name = "BottomFDI_Pi"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.bottomFDI_PiThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kKapton, kVisible, kYellow)

    def makeChipGlue(self):
        name = "ChipGlue"
        size = Size(self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kGlue, kVisible, kYellow)

    def makeAlCables(self):
        name = "AlCable"
        size = Size(self.chipSize.dx, self.chipSize.dy, self.AlCableThickness)
        pos = self.getPosDepth(self.chipPos.x, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kBlue)

    def makeActiveChips(self):
        name ="ActiveChip"
        size = Size(self.chipSize.dx, self.chipSize.dy, self.activeChipThickness)
        pos = self.getPosDepth(self.chipPos.x, self.chipPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kRed)

    def recalculateSize(self):
        self.thickness = 2*self.glueThickness + self.kaptonThickness + self.AlCableThickness + self.activeChipThickness + self.absorberThickness + self.absorberAirGapThickness
        self.moduleThickness = 2 * self.glueThickness + self.kaptonThickness + self.AlCableThickness + self.activeChipThickness

        self.size = Size(self.chipSize.dx, self.chipSize.dy, self.thickness)

    def addRotation(self):
        filename = "{}.placements".format(self.name)
        rotation = 180
        translationA = "0 {} {}".format(-self.chipGap.dy/2, -(self.getLastObjectRightSide() - self.getLeftSideNoAbsorber())/2.)
        translationB = "0 {} {}".format( self.chipGap.dy/2,  (self.getLastObjectRightSide() - self.getLeftSideNoAbsorber())/2.)
        self.objects.append(Rotation(self.name, filename, rotation, translationA, translationB))
        print("\033[1mAdded rotation of {} with {} degrees and {} mm +- translation\033[0m".format(self.name, rotation, translationA))

    def addCopies(self, nlayers):
        translation = "0 0 {}".format(self.size.dz)
        self.objects.append(Copy("Layer", nlayers, translation))
        print("\033[1mAdded repeater of {} with {} layers with translation {}\033[0m".format(self.name, nlayers, translation))
   
    def moveAllBoxes(self, fromZ = 0):
        leftSide = self.getLeftSide().z + fromZ
        print("Moving all boxes in {} to start from {} to {} with fromZ {}.".format(self.name, self.getLeftSide().z, leftSide, fromZ))
        print self.objects[0].pos.z - self.objects[0].size.dz/2
        for obj in self.objects:
            if obj.type == "Box":
                obj.move(Pos(0, 0, fromZ))

        print self.objects[0].pos.z - self.objects[0].size.dz/2

    def reduceAllSizes(self, reduceInUm):
        sizeDiff = Size(reduceInUm * um , reduceInUm * um, reduceInUm * um)
        for obj in self.objects:
            if obj.type == "Box" and not "Layer" in obj.name:
                if "Module" in obj.name:
                    obj.size.dz -= sizeDiff.dz/2.
                else:
                    obj.size.dz -= sizeDiff.dz


    def moveAllModuleDaughters(self):
        modulePos = 0

        for obj in self.objects:
            if obj.type == "Box":
                if "Module" in obj.name:
                    modulePos = obj.pos.z
                    print obj.name, modulePos
                    break

        for obj in self.objects:
            if obj.type == "Box":
                if "Module" in obj.mother:
                    print "Moving ", obj.name, " from ", obj.pos.z, 
                    obj.pos.z -= modulePos
                    print "to ", obj.pos.z

    def getModulePosition(self):
        modulePos = 0
        
        for obj in self.objects:
            if obj.type == "Box":
                if "Module" in obj.name:
                    modulePos = obj.pos
                    break

        return modulePos

    def getLayerPosition(self):
        layerPos = 0
        for obj in self.objects:
            if obj.type == "Box":
                if "Layer" in obj.name:
                    layerPos = obj.pos
                    break;

        return layerPos

    def getLayerSize(self):
        layerSize = 0
        for obj in self.objects:
            if obj.type == "Box":
                if "Layer" in obj.name:
                    layerSize = obj.size
                    break

        return layerSize

    def getModuleSize(self):
        moduleSize = 0
        
        for obj in self.objects:
            if obj.type == "Box":
                if "Module" in obj.name:
                    moduleSize = obj.size
                    break
        
        return moduleSize

    def moveAllLayerDaughters(self, doMoveLayer = False, moveLayerTo = 0):
        layerSize = 0
        layerPos = 0

        for obj in self.objects:
            if obj.type == "Box":
                if "Layer" in obj.name:
                    layerPos = obj.pos.z
                    if doMoveLayer:
                        print("Moving LAYER {} from {} to {}.".format(obj.name, obj.pos.z, moveLayerTo))
                        obj.pos.z = moveLayerTo
                        
                    layerSize = obj.size.dz

        for obj in self.objects:
            if obj.type == "Box":
                if "Layer" in obj.mother:
                    print "Moving ", obj.name, " from ", obj.pos.z, 
                    obj.pos.z += - layerPos
                    print "to ", obj.pos.z


    def buildModule(self):
        self.recalculateSize()
        self.makeMotherLayer()
        self.makeMotherModule()
        self.makeActiveChips()
        self.makeKapton()
        self.makeAlCables()
        self.makeKaptonGlue()
        self.makeChipGlue()
        self.makeAbsorber(1) #, kIsFirstAbsorber)
        self.makeAbsorberAirGap()
        self.recalculateSize()
        if (self.nlayers>1):
            self.addCopies(self.nlayers)
        self.moveAllBoxes(self.fromZ)

    def printGeometry(self):
        for obj in self.objects:
            obj.printBox()

    def clearFile(self):
        f = open("Module.mac", "w")
        f.write("")
        f.close()

    def correctAllNames(self):
        for obj in self.objects:
            if "FM_" in obj.name: continue

            if obj.name != "FirstModule":
                obj.name = "FM_" + obj.name
            if obj.type == "Box":
                if obj.mother == "Layer":
                    obj.mother = "FM_" + obj.mother

    def writeGeometry(self):
        for obj in self.objects:
            obj.writeBox()

    def plotGeometry(self, fig = None):
        if not fig:
            fig = plt.figure(figsize=(15,6.5))
        
        ax1 = fig.add_subplot(121)

        # (X,Y) plot
        for obj in self.objects:
            if obj.type == "Rotation" or obj.type == "Copy" or obj.name == "Absorber2": continue

            x1 = obj.pos.x - obj.size.dx/2
            y1 = obj.pos.y - obj.size.dy/2
            if obj.color:
                ax1.add_patch(patches.Rectangle((x1, y1), obj.size.dx, obj.size.dy, facecolor=obj.color))
            else:
                ax1.add_patch(patches.Rectangle((x1, y1), obj.size.dx, obj.size.dy, fill=False))

        minZ = self.getLeftSideNoAbsorber()
        maxZ = self.getLastObjectRightSide()

        print "minZ", minZ, "maxZ", maxZ

        ax2 = fig.add_subplot(122)

        # (Z,Y) plot
        for obj in self.objects:
            if obj.type == "Rotation" or obj.type == "Copy": continue
            x1 = obj.pos.x - obj.size.dx/2
            y1 = obj.pos.y - obj.size.dy/2
            z1 = obj.pos.z - obj.size.dz/2

            print("Drawing object {} from Z {} to Z {}.".format(obj.name, z1, z1+obj.size.dz))

            if obj.color:
                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, facecolor=obj.color))
            else:
                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, fill=False))

        # DRAW ROTATION
        for rot in self.objects:
            if rot.type == "Rotation":
                print("Starting rotation for {}".format(self.name))
                for obj in self.objects:
                    if obj.type == "Box" and not "Absorber" in obj.name and not "Layer" in obj.name:
                        pos = obj.pos
                        size = obj.size
                        newPos = obj.rotatePosition(rot.rotationmatrix, rot.rotation, self.chipGap.dy, -2, 2, minZ, maxZ)
                        
                        x1 = newPos.x - obj.size.dx/2.
                        y1 = newPos.y - obj.size.dy/2.
                        z1 = newPos.z - obj.size.dz/2.

                        print("Adding {} {} {} {}->{} for rotation".format(obj.name, x1, y1, z1, z1+obj.size.dz))

                        if obj.color:
                            ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, facecolor=obj.color))
                        else:
                            ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, fill=False))

        # DRAW COPIES
        for cop in self.objects:
            if cop.type == "Copy":
                for obj in self.objects:
                    if obj.type == "Box":
                        translation = float(cop.translation.split()[2])
                        for i in range(cop.nlayers):
                            pos = obj.pos
                            size = obj.size
                            newPos = Pos(pos.x, pos.y, pos.z + (i+1) * translation)
                            
                            x1 = newPos.x - obj.size.dx/2.
                            y1 = newPos.y - obj.size.dy/2.
                            z1 = newPos.z - obj.size.dz/2.

                            if obj.color:
                                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, facecolor=obj.color))
                            else:
                                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, fill=False))

        ax1.set_xlim([-30, 30])
        ax1.set_ylim([-30, 30])
        ax2.set_xlim([-5,  10])
        ax2.set_ylim([-30, 30])

        return fig

class scanner:
    def __init__(self, i=1, dx = 300*cm, dy = 300*cm, dz = 150*cm):
        self.size = Size(dx, dy, dz)
        self.name = "scanner_{}".format(i)
        self.mother = "world"
        self.pos = Pos(0, 0, 40*cm)
        self.box = Box("{}".format(self.name), self.mother, self.pos, self.size)

    def printGeometry(self):
        self.box.printBox()

    def writeGeometry(self):
        self.box.writeBox()

class World:
    def __init__(self, dx = 350*cm, dy = 350*cm, dz = 350*cm):
        self.name = "world"
        self.size = Size(dx, dy, dz)
        self.pos = Pos()
        self.box = Box(self.name, kNoMother, self.pos, self.size)

    def printGeometry(self):
        self.box.printBox()

    def writeGeometry(self):
        self.box.writeBox()

class Degrader:
   def __init__(self, dx = 150*cm, dy = 150*cm, dz = "{degraderthickness}"):
      self.name = "degrader"
      self.size = Size(dx, dy, dz)
      self.mothersize = Size(dx*1.1, dy*1.1, dz);
      self.mothername = "phantom"
      self.mother = "world"
      self.material = "Water"
      self.pos = Pos(0, 0, 0) # maximum 2*20 = 40 cm WEPL
      self.motherpos = Pos(0, 0, "{halfdegraderthickness}")
      self.motherbox = Box(self.mothername, self.mother, self.motherpos, self.mothersize)
      self.box = Box(self.name, self.mothername, self.pos, self.size, self.material)

   def printGeometry(self):
      self.motherbox.printBox()
      self.box.printBox()

   def writeGeometry(self):
      self.motherbox.writeBox()
      self.box.writeBox()


class MainMenu(Frame):
    def __init__(self, parent, firstModule, module, world, degrader, scanner1, scanner2):
        Frame.__init__(self, parent)
        self.parent = parent
        self.parent.protocol("WM_DELETE_WINDOW", self.myQuit)
        self.parent.title("GATE Geometry Builder")

        self.returnDictionary = {}
        self.firstModule = firstModule
        self.module = module
        self.world = world
        self.degrader = degrader
        self.scanner1 = scanner1
        self.scanner2 = scanner2
        self.state = None

        self.firstLayerPosition = 0

        self.topContainer = Frame(self, borderwidth=10)
        self.leftContainer = Frame(self, borderwidth=10)
        self.leftLabelContainer = Frame(self.leftContainer, borderwidth=5)
        self.leftEntryContainer = Frame(self.leftContainer, borderwidth=5)
        self.rightContainer = Frame(self, borderwidth=10)
        self.bottomContainer = Frame(self, borderwidth=10)

        self.bottomContainer.pack(side=BOTTOM)
        self.topContainer.pack()
        self.leftContainer.pack(side=LEFT)
        self.leftLabelContainer.pack(side=LEFT)
        self.leftEntryContainer.pack(side=LEFT)
        self.rightContainer.pack(side=RIGHT)

        self.var_xsize = DoubleVar()
        self.var_ysize = DoubleVar()
        self.var_xgap = DoubleVar()
        self.var_ygap = DoubleVar()
        self.var_nlayers = IntVar()
        self.var_material = StringVar()
        self.var_absorberThickness = DoubleVar()
        self.var_firstmaterial = StringVar()
        self.var_firstAbsorberThickness = DoubleVar()
        self.var_chipThickness = DoubleVar()
        self.var_useDegrader = IntVar()

        self.labeltext_xsize = StringVar()
        self.labeltext_ysize = StringVar()
        self.labeltext_xgap = StringVar()
        self.labeltext_ygap = StringVar()
        self.labeltext_nlayers = StringVar()
        self.labeltext_material = StringVar()
        self.labeltext_absorberThickness = StringVar()
        self.labeltext_firstmaterial = StringVar()
        self.labeltext_firstAbsorberThickness = StringVar()
        self.labeltext_chipThickness = StringVar()

        self.var_nlayers.set(150)
        self.var_xsize.set(270) # was 19.145
        self.var_ysize.set(135) # was 19.145
        self.var_xgap.set(0)
        self.var_ygap.set(0)
        self.var_material.set("Aluminium")
        self.var_absorberThickness.set(2)
        self.var_firstmaterial.set("Air")
        self.var_firstAbsorberThickness.set(2)
        self.var_chipThickness.set(14)
        self.var_useDegrader.set(1)
        
        self.optionsTitle = Label(self.leftContainer, text="Geometry options")

        self.labeltext_xsize.set("Chip x size [mm]")
        self.labeltext_ysize.set("Chip y size [mm]")
        self.labeltext_xgap.set("Gap between X chips [mm]")
        self.labeltext_ygap.set("Gap between Y chips [mm]")
        self.labeltext_nlayers.set("Number of layers")
        self.labeltext_material.set("Absorber material")
        self.labeltext_absorberThickness.set("Absorber thickness [mm]")
        self.labeltext_firstmaterial.set("First absorber material")
        self.labeltext_firstAbsorberThickness.set("First absorber thickness [mm]")
        self.labeltext_chipThickness.set("Sensor chip thickness [um]")

        labelHeight = 25

        self.label_xsize = Entry(self.leftLabelContainer, textvariable=self.labeltext_xsize, state=DISABLED, width=labelHeight)
        self.label_ysize = Entry(self.leftLabelContainer, textvariable=self.labeltext_ysize, state=DISABLED, width=labelHeight)
        self.label_xgap = Entry(self.leftLabelContainer, textvariable=self.labeltext_xgap, state=DISABLED, width=labelHeight)
        self.label_ygap = Entry(self.leftLabelContainer, textvariable=self.labeltext_ygap, state=DISABLED, width=labelHeight)
        self.label_nlayers = Entry(self.leftLabelContainer, textvariable=self.labeltext_nlayers, state=DISABLED, width=labelHeight)
        self.label_material = Entry(self.leftLabelContainer, textvariable=self.labeltext_material, state=DISABLED, width=labelHeight)
        self.label_absorberThickness = Entry(self.leftLabelContainer, textvariable=self.labeltext_absorberThickness, state=DISABLED, width=labelHeight)
        self.label_firstmaterial = Entry(self.leftLabelContainer, textvariable=self.labeltext_firstmaterial, state=DISABLED, width=labelHeight)
        self.label_firstAbsorberThickness = Entry(self.leftLabelContainer, textvariable=self.labeltext_firstAbsorberThickness, state=DISABLED, width=labelHeight)
        self.label_chipThickness = Entry(self.leftLabelContainer, textvariable=self.labeltext_chipThickness, state=DISABLED, width=labelHeight)
        
        self.entry_xsize = Entry(self.leftEntryContainer, textvariable=self.var_xsize, width=15)
        self.entry_ysize = Entry(self.leftEntryContainer, textvariable=self.var_ysize, width=15)
        self.entry_xgap = Entry(self.leftEntryContainer, textvariable=self.var_xgap, width=15)
        self.entry_ygap = Entry(self.leftEntryContainer, textvariable=self.var_ygap, width=15)
        self.entry_nlayers = Entry(self.leftEntryContainer, textvariable=self.var_nlayers, width=15)
        self.entry_material = Entry(self.leftEntryContainer, textvariable=self.var_material, width=15)
        self.entry_absorberThickness = Entry(self.leftEntryContainer, textvariable=self.var_absorberThickness, width=15)
        self.entry_firstmaterial = Entry(self.leftEntryContainer, textvariable=self.var_firstmaterial, width=15)
        self.entry_firstAbsorberThickness = Entry(self.leftEntryContainer, textvariable=self.var_firstAbsorberThickness, width=15)
        self.entry_chipThickness = Entry(self.leftEntryContainer, textvariable=self.var_chipThickness, width=15)
        
        self.check_useDegrader = Checkbutton(self.topContainer, text="Use water degrader phantom", variable=self.var_useDegrader)

        self.label_xsize.pack()
        self.entry_xsize.pack()
        self.label_ysize.pack()
        self.entry_ysize.pack()
        self.label_xgap.pack()
        self.entry_xgap.pack()
        self.label_ygap.pack()
        self.entry_ygap.pack()
        self.label_nlayers.pack()
        self.entry_nlayers.pack()
        self.label_material.pack()
        self.entry_material.pack()
        self.label_absorberThickness.pack()
        self.entry_absorberThickness.pack()
        self.label_firstmaterial.pack()
        self.entry_firstmaterial.pack()
        self.label_firstAbsorberThickness.pack()
        self.entry_firstAbsorberThickness.pack()
        self.label_chipThickness.pack()
        self.entry_chipThickness.pack()
        self.check_useDegrader.pack(side=LEFT)

        self.button_create = Button(self.bottomContainer, text="Show geometry", command=self.show_geometry, width=30)
        self.button_create.pack(side=LEFT)
        self.button_write = Button(self.bottomContainer, text="Write files", command=self.write_files, width=30)
        self.button_write.pack(side=LEFT)
        self.button_exit = Button(self.bottomContainer, text="QUIT", command=self.myQuit, width=30)
        self.button_exit.pack(side=LEFT)

        self.pack()

    def print_geometry(self):
        self.readyModule()
        self.firstModule.printGeometry()
        self.module.printGeometry()

    def write_files(self):
        self.firstModule.clearFile()

        self.readyModule()
       
        self.firstModule.moveAllModuleDaughters()
        self.module.moveAllModuleDaughters()
        
        self.firstModule.moveAllLayerDaughters(True, 0)
        self.module.moveAllLayerDaughters(True, self.firstLayerPosition)

        self.world.writeGeometry()
        if self.degrader and self.var_useDegrader.get():
           self.degrader.writeGeometry()

        self.scanner1.writeGeometry()
        self.scanner2.writeGeometry()

        self.firstModule.writeGeometry()
        self.module.writeGeometry()

    def show_geometry(self):
        self.readyModule()
        figure = self.firstModule.plotGeometry()
        figure = self.module.plotGeometry(figure)
        self.drawPlot(figure)

    def readyModule(self):
        self.returnDictionary = {"xsize" : self.var_xsize.get(), "ysize" : self.var_ysize.get(), "nlayers" : self.var_nlayers.get(), "xgap" : self.var_xgap.get(), "ygap" : self.var_ygap.get(), "chipthickness" : self.var_chipThickness.get(), "material" : self.var_material.get(), "absorberthickness" : self.var_absorberThickness.get(), "firstmaterial" : self.var_firstmaterial.get(), "firstabsorberthickness" : self.var_firstAbsorberThickness.get()}
        
        self.firstModule.nlayers = 1
        self.firstModule.chipSize = Size(self.returnDictionary["xsize"], self.returnDictionary["ysize"])
        self.firstModule.chipGap = Size(self.returnDictionary["xgap"], self.returnDictionary["ygap"])
        self.firstModule.absorberMaterial = self.returnDictionary["firstmaterial"]
        self.firstModule.absorberThickness = self.returnDictionary["firstabsorberthickness"]
        self.firstModule.activeChipThickness = self.returnDictionary["chipthickness"] / 1000.
#        self.firstModule.firstmaterial = self.returnDictionary["firstmaterial"]
#        self.firstModule.firstAbsorberThickness = self.returnDictionary["firstabsorberthickness"]
        self.firstModule.recalculateSize()
        self.firstModule.buildModule()
        lastPos = self.firstModule.getLayerThickness()
        self.firstModule.fromZ = -lastPos/2.
        self.firstModule.moveAllBoxes(self.firstModule.fromZ)
        self.firstModule.reduceAllSizes(0.1)
        self.firstModule.correctAllNames()

        self.module.nlayers = self.returnDictionary["nlayers"]
        self.module.chipSize = Size(self.returnDictionary["xsize"], self.returnDictionary["ysize"])
        self.module.chipGap = Size(self.returnDictionary["xgap"], self.returnDictionary["ygap"])
        self.module.absorberMaterial = self.returnDictionary["material"]
        self.module.absorberThickness = self.returnDictionary["absorberthickness"]
        self.module.activeChipThickness = self.returnDictionary["chipthickness"] / 1000.
#        self.module.firstAbsorberThickness= self.returnDictionary["firstabsorberthickness"]
        self.module.recalculateSize()
        self.module.buildModule()
        lastPos = self.module.getLayerThickness()
        self.module.fromZ = -lastPos/2 + self.firstModule.getLayerThickness()/2.
        self.module.moveAllBoxes(self.module.fromZ)
        self.module.reduceAllSizes(0.1)
        
        self.scanner1.size.dx = self.firstModule.size.dx
        self.scanner1.size.dy = self.firstModule.size.dy
        self.scanner2.size.dx = self.firstModule.size.dx
        self.scanner2.size.dy = self.firstModule.size.dy

        self.scanner1.size.dz = self.firstModule.thickness
        self.scanner1.pos.z = self.scanner1.size.dz/2
        self.scanner2.size.dz = self.module.thickness * self.module.nlayers*2
        self.scanner2.pos.z = self.scanner2.size.dz/2. + self.scanner1.pos.z + self.scanner1.size.dz/2.
        self.firstLayerPosition = -self.scanner2.size.dz/2 + self.module.thickness/2
        
        print("Scanner 1 position: {} and size: {}. Scanner 2 position: {} and size: {}. First layer position: {} and size: {}".format(self.scanner1.pos, self.scanner1.size, self.scanner2.pos, self.scanner2.size, self.firstLayerPosition, self.firstModule.thickness))

    def getState(self):
        return self.state

    def getReturnDictionary(self):
        return self.returnDictionary

    def drawPlot(self, figure):
        canvas = FigureCanvasTkAgg(figure, master=self.rightContainer)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, self.rightContainer)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

    def myQuit(self):
        self.parent.destroy()
        self.quit()

def main():
    # See if we can get the rotation right
    # fix names of mother volumes
    # Cross check geometry and Module.mac (filename should be constant)

    world = World()
    degrader = Degrader()
    scanner1 = scanner(1)
    scanner2 = scanner(2)

    firstModule = Module("FirstModule", scanner1.name, Pos(0,0,0), kAir, 1)
    module = Module("Module", scanner2.name, Pos(0,0,0), kAluminium, 20)

    root = Tk()
    mainmenu = MainMenu(root, firstModule, module, world, degrader, scanner1, scanner2)
    root.mainloop()

main()

