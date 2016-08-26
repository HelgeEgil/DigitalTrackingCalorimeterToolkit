from enum import *

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

class Box:
    def __init__(self, name = False, mother = False, pos = False, size = False, material = kNoMaterial, visible = kNotVisible, color = kNoColor):
        self.name = name
        self.mother = mother
        self.pos = pos
        self.size = size
        self.material = material
        self.visible = visible
        self.type = "Box"
        self.color = color

    def printBox(self):
        if self.mother:
            print("/gate/{}/daughters/name {}".format(self.mother, self.name))
            print("/gate/{}/daughters/insert box".format(self.mother))
        print("/gate/{}/geometry/setXLength {} mm".format(self.name, self.size.dx))
        print("/gate/{}/geometry/setYLength {} mm".format(self.name, self.size.dy))
        print("/gate/{}/geometry/setZLength {} mm".format(self.name, self.size.dz))
        if self.pos.isPos:
            print("/gate/{}/placement/setTranslation {} {} {} mm".format(self.name, self.pos.x, self.pos.y, self.pos.z))
        if self.material:
            print("/gate/{}/setMaterial {}".format(self.name, self.material))
        if self.color:
            print("/gate/{}/setColor {}".format(self.name, self.color))
        if not self.visible:
            print("/gate/{}/vis/setVisible false".format(self.name))
        print("")

    def move(pos):
        self.pos += pos

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
        self.size = Size(5*cm, 5*cm, 3.975*mm)
        self.name = name
        self.chipSize = Size(1914.5 * um, 1914.5 * um)
        self.chipGap = Size(100. * um, -90 * um)
        self.chipLeftPos  = Pos(-self.chipSize.dx/2 - self.chipGap.dx/2, self.chipSize.dy/2)
        self.chipRightPos = Pos( self.chipSize.dx/2 + self.chipGap.dx/2, self.chipSize.dy/2)
        self.fillerPos    = Pos(0, -self.chipSize.dy/2)
        self.absorberThickness = 1.5 * mm
        self.glueThickness = 40 * um
        self.pcbThickness = 160 * um
        self.passiveChipThickness = 106 * um
        self.activeChipThickness = 14 * um
        self.fillerGlueThickness = 70 * um
        self.fillerThickness = 300 * um
        self.mother = mother
        self.absorberMaterial = absorberMaterial
        self.firstAbsorberMaterial = firstAbsorberMaterial

    def setChipSizes(self, dx, dy):
        self.chipSize = Size(dx, dy)
        self.chipLeftPos  = Pos(-self.chipSize.dx/2 - self.chipGap.dx/2, self.chipSize.dy/2)
        self.chipRightPos = Pos( self.chipSize.dx/2 + self.chipGap.dx/2, self.chipSize.dy/2)

    def setChipGaps(self, dx, dy):
        self.chipGap = Size(dx, dy)
        self.chipLeftPos  = Pos(-self.chipSize.dx/2 - self.chipGap.dx/2, self.chipSize.dy/2)
        self.chipRightPos = Pos( self.chipSize.dx/2 + self.chipGap.dx/2, self.chipSize.dy/2)

    def setThicknesses(self, absorber, glue, pcb, passive, active, fillerglue, filler):
        self.absorberThickness = absorber
        self.glueThickness = glue
        self.pcbThickness = pcb
        self.activeChipThickness = active
        self.passiveChipThickness = passive
        self.fillerGlueThickness = fillerglue
        self.fillerThickness = filler

    def getLeftSide(self):
        return self.pos - self.size / 2.

    def getRightSide(self):
        return self.pos + self.size / 2.

    def getLastObjectCenter(self):
        sumPoint = self.getLeftSide()
        for obj in self.objects:
            sumPoint += obj

        return sumPoint

    def getThisObjectCenter(self, box):
        """Add objects from left to the right IFF the object overlaps in x,y with this object"""

        sumPoint = self.getLeftSide()

        print("Searching for overlaps with {}.".format(box))
        for obj in self.objects:
            if obj.Contains(box) and obj.name != self.name:
                sumPoint.z = obj.pos.z + obj.size.dz/ 2 + box.size.dz/2.
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
        self.objects.append(Box(self.name, self.mother, self.pos, self.size))
        print("Added {}".format(self.objects[-1]))

    def makeMotherLayer(self):
        self.objects.append(Box("Layer", self.name, self.pos, self.size))
        print("Added {}".format(self.objects[-1]))

    def makeAbsorber(self):
        name = "Absorber"
        size = Size(5 * cm, 5 * cm, 1.5 * mm)
        pos = self.getPosDepth(0, 0, size, name)

        material = self.absorberMaterial
        if self.firstAbsorberMaterial:
            material = self.firstAbsorberMaterial
        
        self.addObject(name, pos, size, material, kVisible, kGrey)

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
    
    def buildModule(self):
        self.makeMotherModule()
        self.makePCBGlue()
        self.makePCB()
        self.makeChipGlue()
        self.makePassiveChips()
        self.makeActiveChips()
        self.makeFillerGlue()
        self.makeFiller()

    def printGeometry(self):
        for obj in self.objects:
            obj.printBox()

class DigitalTrackingCalorimeter:
    def __init__(self, dx = 15*cm, dy = 15*cm, dz = 80*cm):
        self.size = Size(dx, dy, dz)
        self.name = "DigitalTrackingCalorimeter"
        self.mother = "World"
        self.pos = Pos()
        self.box = Box("DigitalTrackingCalorimeter", self.mother, self.pos, self.size)

    def printBox(self):
        self.box.printBox()

class World:
    def __init__(self, dx = 100*cm, dy = 100*cm, dz = 100*cm):
        self.name = "World"
        self.size = Size(dz, dy, dz)
        self.pos = Pos()
        self.box = Box("World", kNoMother, self.pos, self.size)

    def printBox(self):
        self.box.printBox()

def main():
    world = World()
    DTC = DigitalTrackingCalorimeter()

    module = Module("Module", DTC.name, Pos(0,0,0), kTungsten, 39)
    module.buildModule()
    module.printGeometry()

main()
