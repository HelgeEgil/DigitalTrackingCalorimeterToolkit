import matplotlib
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
kPCB = "PCB"
kGlue = "Glue"
kChip = "Silicon"
kAir = "Air"

kNoMaterial = False
kNoMother = False
kVisible = True
kNotVisible = False

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
        dx = self.dx + other.dx
        dy = self.dy + other.dy
        dz = self.dz + other.dz
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

    def rotatePosition(self, mat, deg, gapY, minY, maxY, minZ, maxZ, i=0, cop=0):
        if mat == "1 0 0" and deg == 180:
            mean = (maxZ + minZ)/2.
            delta = maxZ - minZ
            z = self.pos.z
            if cop:
                z += (i+1) * float(cop.translation.split()[2])
            newPos = Pos(self.pos.x, -(self.pos.y + gapY) + (maxY + minY), -z + (maxZ + minZ) + delta)
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
    def __init__(self, name, mother, pos, absorberMaterial, nlayers = 40, firstAbsorberMaterial = kNoMaterial):
        self.objects = []
        self.nlayers = 40
        self.pos = pos
        self.size = Size(6*cm, 6*cm, 3.975*mm)
        self.name = name
        self.chipSize = Size(19145 * um, 19145 * um)
        self.chipGap = Size(100. * um, -90 * um)
        self.chipLeftPos  = Pos(-self.chipSize.dx/2 - self.chipGap.dx/2, self.chipSize.dy/2)
        self.chipRightPos = Pos( self.chipSize.dx/2 + self.chipGap.dx/2, self.chipSize.dy/2)
        self.fillerPos    = Pos(0, -self.chipSize.dy/2)
        self.absorberThickness = 1.5 * mm
        self.glueThickness = 40 * um
        self.pcbThickness = 160 * um
        self.passiveChipThickness = 106 * um
        self.activeChipThickness = 14 * um
        self.airGapAfterChipThickness = 90 * um
        self.airGapAfterFillerThickness = 80 * um
        self.fillerGlueThickness = 70 * um
        self.fillerThickness = 300 * um
        self.mother = mother
        self.absorberMaterial = absorberMaterial
        self.firstAbsorberMaterial = firstAbsorberMaterial

    def getLeftSide(self):
        return self.pos - self.size / 2.

    def getRightSide(self):
        return self.pos + self.size / 2.

    def getLeftSideNoAbsorber(self):
        for obj in self.objects:
            if obj.name != "Layer" and obj.name != self.name and obj.name[:-1] != "Absorber":
                return obj.pos.z - obj.size.dz / 2

    def getLastObjectCenter(self):
        sumPoint = self.getLeftSide()
        for obj in self.objects:
            sumPoint += obj

        return sumPoint

    def getLastObjectRightSide(self):
        rightSide = -1000
        for obj in self.objects:
            if obj.name != self.name and obj.name != "Layer" and obj.type == "Box" and obj.name[:-1] != "Absorber":
                rightSide = max(rightSide, obj.pos.z + obj.size.dz/2)

        return rightSide

    def getThisObjectCenter(self, box):
        """Add objects from left to the right IFF the object overlaps in x,y with this object"""

        sumPoint = self.getLeftSide()
        sumPoint.z += box.size.dz / 2.

        print("Searching for overlaps with {}.".format(box))
        for obj in self.objects:
            if obj.type != "Box": continue
            if obj.Contains(box) and obj.name != self.name and obj.name != "Layer":
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
        self.objects.append(Box(self.name, "Layer", self.pos, self.size))
        self.objects[-1].genericRepeater = True
        print("Added {}".format(self.objects[-1]))

    def makeMotherLayer(self):
        self.objects.append(Box("Layer", self.mother, self.pos, self.size))
        print("Added {}".format(self.objects[-1]))

    def makeAbsorber(self, n):
        name = "Absorber" + str(n)
        size = Size(5 * cm, 5 * cm, 1.5 * mm)
        pos = self.getPosDepth(0, 0, size, name)
        if n == 2: 
            pos.z += self.getLastObjectRightSide() - self.getLeftSideNoAbsorber()

        material = self.absorberMaterial
        if self.firstAbsorberMaterial:
            material = self.firstAbsorberMaterial
        
        self.addAbsorber(name, pos, size, material, kVisible, kGrey)

    def makePCBGlue(self):
        name = "PCBGlue"
        size = Size(2 * self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kGlue, kVisible, kYellow)

    def makePCB(self):
        name = "PCB"
        size = Size(2 * self.chipSize.dx + self.chipGap.dx, self.chipSize.dx, self.pcbThickness)
        pos = self.getPosDepth(0, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kPCB, kVisible, kGreen)

    def makeChipGlue(self):
        name = "ChipGlue"
        size = Size(2 * self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.glueThickness)
        pos = self.getPosDepth(0, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kGlue, kVisible, kYellow)

    def makePassiveChips(self):
        name = "PassiveChipLeft"
        size = Size(self.chipSize.dx, self.chipSize.dy, self.passiveChipThickness)
        pos = self.getPosDepth(self.chipLeftPos.x, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kBlue)

        name = "PassiveChipRight"
        pos = self.getPosDepth(self.chipRightPos.x, self.chipRightPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kBlue)

    def makeActiveChips(self):
        name ="ActiveChipLeft"
        size = Size(self.chipSize.dx, self.chipSize.dy, self.activeChipThickness)
        pos = self.getPosDepth(self.chipLeftPos.x, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kRed)
        
        name = "ActiveChipRight"
        pos = self.getPosDepth(self.chipRightPos.x, self.chipRightPos.y, size, name)
        self.addObject(name, pos, size, kChip, kVisible, kRed)

    def makeAirGapAfterChip(self):
        name = "AirGapAfterChip"
        size = Size(2*self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.airGapAfterChipThickness)
        pos = self.getPosDepth(0, self.chipLeftPos.y, size, name)
        self.addObject(name, pos, size, kAir, kNotVisible, kNoColor)

    def makeFillerGlue(self):
        name = "FillerGlue"
        size = Size(2 * self.chipSize.dx + self.chipGap.dx , self.chipSize.dy, self.fillerGlueThickness)
        pos = self.getPosDepth(self.fillerPos.x, self.fillerPos.y, size, name)
        self.addObject(name, pos, size, kGlue, kVisible, kYellow)

    def makeFiller(self):
        name = "Filler"
        size = Size(2 * self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.fillerThickness)
        pos = self.getPosDepth(self.fillerPos.x, self.fillerPos.y, size, name)
        self.addObject(name, pos, size, self.absorberMaterial, kVisible, kGrey)
    
    def makeAirGapAfterFiller(self):
        name = "AirGapAfterFiller"
        size = Size(2*self.chipSize.dx + self.chipGap.dx, self.chipSize.dy, self.airGapAfterFillerThickness)
        pos = self.getPosDepth(self.fillerPos.x, self.fillerPos.y, size, name)
        self.addObject(name, pos, size, kAir, kNotVisible, kNoColor)

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
        leftSide = self.getLeftSide().z - fromZ
        print("Moving all boxes in {} to start from {} to {}.".format(self.name, self.getLeftSide().z, leftSide))
        for obj in self.objects:
            if obj.type == "Box":
                obj.move(Pos(0, 0, -leftSide))

    def buildModule(self):
        self.makeMotherLayer()
        self.makeAbsorber(1)
        self.makeMotherModule()
        self.makePCBGlue()
        self.makePCB()
        self.makeChipGlue()
        self.makePassiveChips()
        self.makeActiveChips()
        self.makeAirGapAfterChip()
        self.makeFillerGlue()
        self.makeFiller()
        self.makeAirGapAfterFiller()
        self.makeAbsorber(2)
        self.addRotation()
        self.addCopies(self.nlayers)
        self.moveAllBoxes()

    def printGeometry(self):
        for obj in self.objects:
            obj.printBox()

    def clearFile(self):
        f = open("Module.mac", "w")
        f.write("")
        f.close()

    def correctAllNames(self):
        for obj in self.objects:
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

        ax2 = fig.add_subplot(122)

        for obj in self.objects:
            if obj.type == "Rotation" or obj.type == "Copy": continue
            x1 = obj.pos.x - obj.size.dx/2
            y1 = obj.pos.y - obj.size.dy/2
            z1 = obj.pos.z - obj.size.dz/2

            if obj.color:
                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, facecolor=obj.color))
            else:
                ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, fill=False))

        for rot in self.objects:
            if rot.type == "Rotation":
                print("Starting rotation for {}".format(self.name))
                for obj in self.objects:
                    if obj.type == "Box" and obj.name[:-1] != "Absorber":
                        pos = obj.pos
                        size = obj.size
                        newPos = obj.rotatePosition(rot.rotationmatrix, rot.rotation, self.chipGap.dy, -2, 2, minZ, maxZ)
                        
                        x1 = newPos.x - obj.size.dx/2.
                        y1 = newPos.y - obj.size.dy/2.
                        z1 = newPos.z - obj.size.dz/2.

                        print("Adding {} {} {} for rotation".format(x1, y1, z1))

                        if obj.color:
                            ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, facecolor=obj.color))
                        else:
                            ax2.add_patch(patches.Rectangle((z1, y1), obj.size.dz, obj.size.dy, fill=False))

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

        ax1.set_xlim([-25, 25])
        ax1.set_ylim([-25, 25])
        ax2.set_xlim([0,  4*self.nlayers])
        ax2.set_ylim([-25, 25])

        return fig

class DigitalTrackingCalorimeter:
    def __init__(self, dx = 15*cm, dy = 15*cm, dz = 80*cm):
        self.size = Size(dx, dy, dz)
        self.name = "DigitalTrackingCalorimeter"
        self.mother = "world"
        self.pos = Pos()
        self.box = Box("DigitalTrackingCalorimeter", self.mother, self.pos, self.size)

    def printBox(self):
        self.box.printBox()

    def writeBox(self):
        self.box.writeBox()

class World:
    def __init__(self, dx = 100*cm, dy = 100*cm, dz = 100*cm):
        self.name = "world"
        self.size = Size(dz, dy, dz)
        self.pos = Pos()
        self.box = Box(self.name, kNoMother, self.pos, self.size)

    def printBox(self):
        self.box.printBox()

    def writeBox(self):
        self.box.writeBox()


class MainMenu(Frame):
    def __init__(self, parent, firstModule, module):
        Frame.__init__(self, parent)
        self.parent = parent
        self.parent.protocol("WM_DELETE_WINDOW", self.myQuit)
        self.parent.title("GATE Geometry Builder")

        self.returnDictionary = {}
        self.firstModule = firstModule
        self.module = module
        self.state = None

        self.leftContainer = Frame(self, borderwidth=10)
        self.leftLabelContainer = Frame(self.leftContainer, borderwidth=5)
        self.leftEntryContainer = Frame(self.leftContainer, borderwidth=5)
        self.rightContainer = Frame(self, borderwidth=10)
        self.bottomContainer = Frame(self, borderwidth=10)

        self.bottomContainer.pack(side=BOTTOM)
        self.leftContainer.pack(side=LEFT)
        self.leftLabelContainer.pack(side=LEFT)
        self.leftEntryContainer.pack(side=LEFT)
        self.rightContainer.pack(side=RIGHT)

        self.var_xsize = DoubleVar()
        self.var_ysize = DoubleVar()
        self.var_xgap = DoubleVar()
        self.var_ygap = DoubleVar()
        self.var_nlayers = IntVar()
        self.var_chipThickness = DoubleVar()
        self.var_material = StringVar()
        self.var_firstmaterial = StringVar()

        self.labeltext_xsize = StringVar()
        self.labeltext_ysize = StringVar()
        self.labeltext_xgap = StringVar()
        self.labeltext_ygap = StringVar()
        self.labeltext_nlayers = StringVar()
        self.labeltext_material = StringVar()
        self.labeltext_firstmaterial = StringVar()
        self.labeltext_chipThickness = StringVar()

        self.var_nlayers.set(40)
        self.var_xsize.set(19.145)
        self.var_ysize.set(19.145)
        self.var_xgap.set(0.1)
        self.var_ygap.set(-0.09)
        self.var_material.set("Tungsten")
        self.var_firstmaterial.set("Aluminium")
        self.var_chipThickness.set(0.014)
        
        self.optionsTitle = Label(self.leftContainer, text="Geometry options")

        self.labeltext_xsize.set("Chip x size [mm]")
        self.labeltext_ysize.set("Chip y size [mm]")
        self.labeltext_xgap.set("Gap between X chips [mm]")
        self.labeltext_ygap.set("Gap between Y chips [mm]")
        self.labeltext_nlayers.set("Number of layers")
        self.labeltext_material.set("Absorber material")
        self.labeltext_firstmaterial.set("First absorber material")
        self.labeltext_chipThickness.set("Sensor chip thickness [mm]")

        labelHeight = 25

        self.label_xsize = Entry(self.leftLabelContainer, textvariable=self.labeltext_xsize, state=DISABLED, width=labelHeight)
        self.label_ysize = Entry(self.leftLabelContainer, textvariable=self.labeltext_ysize, state=DISABLED, width=labelHeight)
        self.label_xgap = Entry(self.leftLabelContainer, textvariable=self.labeltext_xgap, state=DISABLED, width=labelHeight)
        self.label_ygap = Entry(self.leftLabelContainer, textvariable=self.labeltext_ygap, state=DISABLED, width=labelHeight)
        self.label_nlayers = Entry(self.leftLabelContainer, textvariable=self.labeltext_nlayers, state=DISABLED, width=labelHeight)
        self.label_material = Entry(self.leftLabelContainer, textvariable=self.labeltext_material, state=DISABLED, width=labelHeight)
        self.label_firstmaterial = Entry(self.leftLabelContainer, textvariable=self.labeltext_firstmaterial, state=DISABLED, width=labelHeight)
        self.label_chipThickness = Entry(self.leftLabelContainer, textvariable=self.labeltext_chipThickness, state=DISABLED, width=labelHeight)
        
        self.entry_xsize = Entry(self.leftEntryContainer, textvariable=self.var_xsize, width=15)
        self.entry_ysize = Entry(self.leftEntryContainer, textvariable=self.var_ysize, width=15)
        self.entry_xgap = Entry(self.leftEntryContainer, textvariable=self.var_xgap, width=15)
        self.entry_ygap = Entry(self.leftEntryContainer, textvariable=self.var_ygap, width=15)
        self.entry_nlayers = Entry(self.leftEntryContainer, textvariable=self.var_nlayers, width=15)
        self.entry_material = Entry(self.leftEntryContainer, textvariable=self.var_material, width=15)
        self.entry_firstmaterial = Entry(self.leftEntryContainer, textvariable=self.var_firstmaterial, width=15)
        self.entry_chipThickness = Entry(self.leftEntryContainer, textvariable=self.var_chipThickness, width=15)

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
        self.label_firstmaterial.pack()
        self.entry_firstmaterial.pack()
        self.label_chipThickness.pack()
        self.entry_chipThickness.pack()

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
        self.readyModule()
        self.firstModule.writeGeometry()
        self.module.writeGeometry()

    def show_geometry(self):
        self.readyModule()
        figure = self.firstModule.plotGeometry()
        figure = self.module.plotGeometry(figure)
        self.drawPlot(figure)

    def readyModule(self):
        self.returnDictionary = {"xsize" : self.var_xsize.get(), "ysize" : self.var_ysize.get(), "nlayers" : self.var_nlayers.get(), "xgap" : self.var_xgap.get(), "ygap" : self.var_ygap.get(), "chipthickness" : self.var_chipThickness.get(), "material" : self.var_material.get(), "firstmaterial" : self.var_firstmaterial.get()}
        
        self.firstModule.nlayers = 1
        self.firstModule.chipSize = Size(self.returnDictionary["xsize"], self.returnDictionary["ysize"])
        self.firstModule.chipGap = Size(self.returnDictionary["xgap"], self.returnDictionary["ygap"])
        self.firstModule.absorberMaterial = self.returnDictionary["material"]
        self.firstModule.activeChipThickness = self.returnDictionary["chipthickness"]
        self.firstModule.firstmaterial = self.returnDictionary["firstmaterial"]
        self.firstModule.buildModule()
        self.firstModule.correctAllNames()
        lastPos = self.firstModule.getLastObjectRightSide()
        print("Lastpost of firstModule is {}.".format(lastPos))
        
        self.module.nlayers = self.returnDictionary["nlayers"]
        self.module.chipSize = Size(self.returnDictionary["xsize"], self.returnDictionary["ysize"])
        self.module.chipGap = Size(self.returnDictionary["xgap"], self.returnDictionary["ygap"])
        self.module.absorberMaterial = self.returnDictionary["material"]
        self.module.activeChipThickness = self.returnDictionary["chipthickness"]
        self.module.buildModule()

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
    DTC = DigitalTrackingCalorimeter()

    firstModule = Module("FirstModule", DTC.name, Pos(0,0,0), kTungsten, 1, kAluminium)
    firstModule.clearFile()
    world.writeBox()
    DTC.writeBox()

    lastPos = firstModule.getLastObjectRightSide()
    module = Module("Module", DTC.name, Pos(0,0,0), kTungsten, 20)
    
    root = Tk()
    mainmenu = MainMenu(root, firstModule, module)
    root.mainloop()

main()

