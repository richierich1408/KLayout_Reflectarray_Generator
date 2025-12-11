`# ----------------------------------------------
# Script Recorded by Ansys Electronics Desktop Version 2024.1.0
# 17:39:47  Mar 19, 2024
# ----------------------------------------------
import ScriptEnv
import math
PI = math.pi
import csv

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.GetActiveProject()
oDesign = oProject.GetActiveDesign()
oEditor = oDesign.SetActiveEditor("3D Modeler")


# Star: Editable variables:
# Do not edit anything outside this box.
# ----------------------------------------------

# Path for the lookup table
dataPath = r'/home/davidhardy/Documents/Transparent Reflectarray/Large Reflectarray Test/Python/LookupTable_SolidRings_SolidGP.csv'

lambda0 = 10.71                         # [mm]
k0 = 2*math.pi / lambda0                # free-space wavenumber
h = 0.5                                 # [mm] - height of the substrate
uc_p = lambda0/4                        # [mm] - period of the unit cell in the x-direction
uc_q = uc_p / math.cos(math.pi/6)       # [mm] - period of the unit cell in the y-direction

# Dimensions of array
N = 19*2+1                                  # ODD number of columns, should equal ~ boundaryRadius*2 / uc_p
M = 21*2+1                                  # ODD number of rows, should equal ~ boundaryRadius*2 / (uc_q + 2*uc_p/sin(pi/6))
boundaryRadius = 25.4*2                   # [mm]
cushion = uc_q/2 + 1                             # [mm] 


# Initialize hexagonal ring parameters        
ringDiagonalLength = 1                  # [mm]
ringWidth = 0.2                         # [mm]
ringInnerDiagonalLength = (2/math.sqrt(3))*(math.sqrt(3)/2 * ringDiagonalLength - 2*ringWidth)

# Coordinates of the horn antenna's apex centre 
x0 = 0        # [mm]
y0 = 0        # [mm]
F = 75.31        # [mm] - F is just the z-coordinate

# Beam-steering 
theta_b = 30 # [degrees]
phi_b = 30   # [degrees]

theta_b = math.radians(theta_b)  # (convert to radians)
phi_b = math.radians(phi_b)  # (convert to radians)

# Lookup Table 
diagonal_vect = []
width_vect = []
phase_vect = []

with open(dataPath, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader)
    for row in reader:
        diagonal_vect = diagonal_vect + [float(row[0])]
        width_vect = width_vect + [float(row[1])]
        phase_vect = phase_vect + [float(row[2])]

# End: editable variables
# ----------------------------------------------


# These equations were for RA with a 1-D phase distribution
# referenceX = -(N-1)/2 * uc_p
# phaseGrad = 2*PI/lambda0 * uc_p * math.sin(theta_bs * PI/180)
# phaseGrad = math.degrees(phaseGrad)

ringNames = "" # this is used to collect all the ring names (E.g. "ring1, ring2") for the unite operation at the end of the script

for n in range(int(N)):
    for m in range(int(M)):
        centerX = (n - (N-1) / 2) * uc_p + (m%2)/2*uc_p
        centerY_p = (m - (M-1) / 2) * (3/4)
        centerY = centerY_p*uc_q
        
        # At the center X coordinate set the phase
        # Equation below was used for 1-D phase distribution array
        # idealPhase = ( phaseGrad * (centerX - referenceX)/uc_p ) % 360

        d = math.sqrt((x0 - centerX) ** 2 + (y0 - centerY) ** 2 + F ** 2) 
        idealPhase = math.degrees(k0 * (d - (centerX * math.cos(phi_b) + centerY * math.sin(phi_b)) * math.sin(theta_b)))
        #idealPhase -= -163.79
        idealPhase = idealPhase % 360 

        # look up the iteration with the closest phase; assume that the csv is organized in order of increasing phase
        matchIndex = 0
        for i in range(len(phase_vect)):
            if phase_vect[i] >= idealPhase:
                matchIndex = i 
                break

        ringDiagonalLength = diagonal_vect[matchIndex]
        ringWidth = width_vect[matchIndex]
        ringInnerDiagonalLength = (2/math.sqrt(3))*(math.sqrt(3)/2 * ringDiagonalLength - 2*ringWidth)
        
        if (centerX**2 + centerY**2) < (boundaryRadius - cushion)**2:
            
            # Create Ring object, then inner ring, and then perform a substraction

            # Outer hexagon
            oEditor.CreateRegularPolygon(
                [
                    "NAME:PolygonParameters",
                    "XCenter:="        , "{} mm".format(0),
                    "YCenter:="        , "{}".format(0),
                    "ZCenter:="        , "{} mm".format(h/2),
                    "XStart:="        , "{} mm".format(0),
                    "YStart:="        , "{} mm".format(ringDiagonalLength/2),
                    "ZStart:="        , "{} mm".format(h/2),
                    "NumSides:="        , "6",
                    "WhichAxis:="        , "Z"
                ],
                [
                    "NAME:Attributes",
                    "Name:="        , "Ring_{}_{}".format(n,m),
                    "Color:="        , "(0 255 0)",
                    "PartCoordinateSystem:=", "Global",
                ])

            # Inner hexagon 
            oEditor.CreateRegularPolygon(
                [
                    "NAME:PolygonParameters",
                    "XCenter:="        , "{} mm".format(0),
                    "YCenter:="        , "{}".format(0),
                    "ZCenter:="        , "{} mm".format(h/2),
                    "XStart:="        , "{} mm".format(0),
                    "YStart:="        , "{} mm".format(ringInnerDiagonalLength/2),
                    "ZStart:="        , "{} mm".format(h/2),
                    "NumSides:="        , "6",
                    "WhichAxis:="        , "Z"
                ],
                [
                    "NAME:Attributes",
                    "Name:="        , "RingInner_{}_{}".format(n,m),
                    "Color:="        , "(0 255 0)",
                    "PartCoordinateSystem:=", "Global",
                ])

            # Subtract the two to make a ring
            oEditor.Subtract(
                [
                    "NAME:Selections"    ,
                    "Blank Parts:="        , "Ring_{}_{}".format(n,m),
                    "Tool Parts:="        , "RingInner_{}_{}".format(n,m)
                ],
                [
                    "NAME:SubtractParameters",
                    "KeepOriginals:=", False 
                ]
            )
                
            # Move the ring to the correct x/y coordinate
            oEditor.Move(
            	[
            		"NAME:Selections",
            		"Selections:="		, "Ring_{}_{}".format(n,m),
            		"NewPartsModelFlag:="	, "Model"
            	], 
            	[
            		"NAME:TranslateParameters",
            		"TranslateVectorX:="	, "{}mm".format(centerX),
            		"TranslateVectorY:="	, "{}mm".format(centerY),
            		"TranslateVectorZ:="	, "0mm"
            	])
            
               
            ringNames += "Ring_{}_{}, ".format(n,m) # not necessary unless doing unite operation at the end


 
