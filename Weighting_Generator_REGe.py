from string import Template
import pathlib
import os
import time
import matplotlib.pyplot as plt
import numpy as np
import csv
import shutil
import concurrent.futures
import pandas as pd
import atexit
import imageio





inlet=3 #State 3 in the concentration file which is the inlet to the medical room/pump
outlet=4 #State 4 in the concentration file which is the outlet from the medical room/pump

current_path = pathlib.Path().resolve()



perBatch=100000
batches=10
detector_loc=inlet  #inlet or outlet
files = os.listdir(current_path)
detector_rad=3.8
detector_area=np.pi*(detector_rad**2)
column_length=240
column_thickness=detector_rad*2
shield_thickness=1.2

def linear_space_string(start, end, num_points):
    # Generate linearly spaced numbers
    numbers = np.linspace(start, end, num_points)
    
    # Convert each number to a string with desired format (optional: format to a certain precision)
    number_strings = [f"{num:.2f}" for num in numbers]  # Format to 2 decimal places
    
    # Join the strings with a space
    result_string = " ".join(number_strings)
    
    return result_string


MAVRIC_input = '''=mavric
Medical Room Detector Response:  Use CE transport
v7.1-28n19g

'-------------------------------------------------------------------------------
' Cross Section Information
' Load the Lead, Air, Copper and Stainless Steel cross secions for use in the geometry
'-------------------------------------------------------------------------------
read comp
    lead 1 end
    dry-air 2 end
    copper 3 end
    ss304 4 end
end comp

'-------------------------------------------------------------------------------
'Geometry Block - SCALE standard geometry package (SGGP)
' A SS pipe filled with the FLiBe source, next to a Lead collimator filled with Air,
' leading to a point detector on the other side. The right side of the collimator is covered in copper
' cuboid collimator because stacking lead blocks is easier
'-------------------------------------------------------------------------------
read geometry
    global unit 1
        cuboid 10   -10  ''' + str(-10+column_length+0.5) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2)  + " " + str(column_thickness/2) + " " + str(-column_thickness/2) + '''
        cuboid 20   -10  ''' + str(-10+column_length) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        cuboid 30   -10   ''' + str(-10+column_length+0.5) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        ycylinder 40    0.635 '''  + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + ''' origin x=-15
        ycylinder 60    0.501 '''  + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + ''' origin x=-15
        cuboid 50  -25.0 ''' + str(-10+column_length + 20) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        xcylinder 70   ''' + str(detector_rad)  + " " + str(-10+column_length+9.5) + " " + str(-10+column_length+10.5) + '''
        



        media 2 1  10
        media 1 1  20 -10
        media 3 1  30 -10 -20
        media 4 1  40 -60
        media 0 1  60 
        media 0 1  70 
        media 2 1  50 -10 -20 -30 -40 -70 -60
    boundary 50
end geometry

'-------------------------------------------------------------------------------
' Definitions Block
'-------------------------------------------------------------------------------
read definitions
'Create Importance Map
    gridGeometry 100
        title="Denovo Grid"
        xLinear 5  -25.0 ''' + str(-10+column_length + 20) + '''
        yLinear 5  ''' + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + '''
        zLinear 5  ''' + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + '''
        xPlanes '''+ linear_space_string(-10,-10 +column_length, 20) +''' end
        yPlanes '''+ linear_space_string(-column_thickness/2,column_thickness/2, 5) +''' end
        zPlanes '''+ linear_space_string(-column_thickness/2,column_thickness/2, 5) +''' end
    end gridGeometry
    
'Position of Point Detector
    location 1  
        position ''' + str(-10 + column_length + 10) + ''' 0.0 0.0   
    end location
'Generate Dose Response as well
    response 5
        title="ANSI standard (1977) photon flux-to-dose-rate factors"
        doseData=9504
    end response


'Use the f71 file created by origen.
'Use designated state (correspondes to location in medical room) for parameter 5 (photon flux)
'*** IMPORTANT. filename needs to be changed to where the .f71 file was saved to***

    distribution 1 
        special="origensBinaryConcentrationFile" 
        filename="''' + str(current_path) + '/Weighting_Generator.f71"' + '''
        parameters 3 5 end 
    end distribution

'Create 1400 Discrete Energy bins from 10 keV to 3000 keV
    energyBounds 1
        linear 1400 3e6 1e4 
    end energyBounds

'The results should be indicative of a 1200 second COunt
    timeBounds 1
        title="linear command"
        linear 1 0.0 1200
    end timeBounds
end definitions

read importancemap
    gridgeometryid=100
    adjointSource 1
        boundingbox -25.0 30.0  10.0 -10.0  10.0 -10.0
        responseid=5
    end adjointSource
'    macromaterial
'        mmsubcell=9
'        mmtolerance=0.001
'        mmRTSpeed
'    end macromaterial
    quadrature=6
    legendre=5
    respweighting
'    xblocks=8
'    yblocks=8
end importancemap



'-------------------------------------------------------------------------------
' Sources Block
'-------------------------------------------------------------------------------
read sources
'Photon source
    src 1
        title="irradiated flibe"
        useNormConst
        ycylinder 0.5 ''' + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        origin x=-15.0 y=0.0 z=0.0
        photons
        eDistributionID=1
    end src
end sources

'-------------------------------------------------------------------------------
' Tallies Block
'-------------------------------------------------------------------------------
read tallies

'Calculte Photon response at desired locaton
    pointDetector 15  
        photon 
        locationID=1 
        responseID=5  
        energyBoundsID=1
        timeBoundsID=1
    end pointDetector
end tallies

'-------------------------------------------------------------------------------
' Parameters Block
'-------------------------------------------------------------------------------
read parameters
    ceLibrary="ce_v7.1_endf.xml"
    randomSeed=0
    photons
    fissionMult=0  secondaryMult=0 
    perBatch=''' + str(perBatch) + ''' batches=''' + str(batches) + '''
end parameters



end data
end'''
    
m = open("Flux_Weighting.inp", "w")
m.write(MAVRIC_input)
m.close()

os.system('scalerte "' + str(current_path) + '/Flux_Weighting.inp"')
#for file in files:
#    if (file.startswith("Flux_Weighting")) and os.path.isfile(file) and file.find("adjoint") == -1:
#        os.unlink(file)