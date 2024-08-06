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





plt.close('all')




inlet=3 #State 3 in the concentration file which is the inlet to the medical room/pump
outlet=4 #State 4 in the concentration file which is the outlet from the medical room/pump

current_path = pathlib.Path().resolve()

detector_rad=3.8
detector_distance=2.5
detector_area=np.pi*(detector_rad**2)
abs_efficiency = pd.read_excel(str(current_path) + '/efficiency.xlsx', sheet_name='Sheet1')
efficiency_energy = abs_efficiency['Energy (keV)']
geom_efficiency=0.5*(1-np.cos(np.arctan(detector_rad/detector_distance)))
efficiency=abs_efficiency['efficiency(%)']/geom_efficiency/100

fuel_salt=0 #1 for yes, 0 for no
num_loops=11
current_loop=0
skip_loops=10
flow_rate=10 #cm/s (Range from 10 to 200)
rad_time=10/flow_rate
looptime=rad_time*1200
count_time=100000 #Seconds
perBatch=10000
batches=1
detector_loc=inlet  #inlet or outlet
precision = ""
num_procs=110
half_procs=30
files = os.listdir(current_path)
column_length=240
column_thickness=detector_rad*2
shield_thickness=1.2

if detector_loc==inlet:
    loc = "inlet"
elif detector_loc==outlet:
    loc = "outlet"

if fuel_salt:
    salt_type="Fuel_Salt"
else:
    salt_type="Clean_Salt"



if perBatch*batches*num_procs>999999:
    csv_path = str(current_path) + "/" + salt_type + "_" + loc + "_" + str(flow_rate) + "cm_per_s__" + str(num_loops) + 'loops_' + str(column_length) + 'cm_column_length.csv'
    
data_vals=[]
# Usage
image_folder = str(current_path)  # Folder where images are stored
output_path = str(current_path) + "/" + salt_type + "_" + loc + "_" + str(flow_rate) + "cm_per_s__" + str(num_loops) + 'loops_' + str(column_length) + 'cm_column_length.mp4'     # Path where the video will be saved
fps = 6  # Frames per second

def cleanup():
    print("Saving data to CSV...")
    with open(csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data_vals)
    print("Data successfully saved to CSV.")
    # Get list of image files, ensuring they are in the correct order
    image_files = [os.path.join(image_folder, img) for img in sorted(os.listdir(image_folder)) if img.endswith(".png")]

    # Create a video writer object
    writer = imageio.get_writer(output_path, fps=fps)

    # Append images to the video writer
    for image_file in image_files:
        image = imageio.imread(image_file)
        writer.append_data(image)
        print(f"Added {image_file} to video")

    writer.close()  # Close the writer to finalize the video
    print("Video created successfully.")

    for file in files:
        if (file.startswith("MAVRIC")) and os.path.isfile(file):
            os.unlink(file)


# Register the cleanup function
atexit.register(cleanup)



# Define unique markers and styles for each element, with their x positions
element_markers = {
    511: ('o', 'Annihilation', 'green'),  # Circle marker, green
    850: ('^', 'Mn-56', 'red'),  # Triangle upwards marker, red
    1368: ('s', 'Na-24', 'blue'),  # Square marker, blue
    1635: ('*', 'F-20', 'purple'),  # Star marker, purple
    1778: ('D', 'Al-28', 'orange'),  # Diamond marker, orange
    1813: ('^', 'Mn-56', 'red'),  # Triangle upwards marker, red
    2114: ('^', 'Mn-56', 'red'),  # Triangle upwards marker, red
    2167: ('p', 'Cl-38', 'cyan'),  # Pentagon marker, cyan
    2754: ('s', 'Na-24', 'blue')  # Square marker, blue
}

if fuel_salt:
    # Adding new elements if Fuel_salt equals 1
    element_markers.update({
        100: ('>', 'La-145', 'teal'),  # Diamond marker, orange
        399: ('v', 'La-144', 'aquamarine'),  # Square marker, blue
        544: ('8', 'La-144 & Kr-90', 'navy'),  # Octagon marker, blue
        600: ('<', 'Cs-140', 'yellow'),  # Square marker, blue
        845: ('v', 'La-144', 'aquamarine'),  # Square marker, blue
        1118: ('h', 'Kr-90', 'maroon'),  # Triangle upwards marker, red
        1313: ('p', 'I-136', 'deeppink'),  # Triangle upwards marker, red
        1428: ('H', 'Sr-94', 'black'),  # Hexagon marker, black
        1565: ('x', 'Br-86', 'fuchsia'),  # Hexagon marker, black
    })

xs_lib = str(current_path) + '/cx.f33"'

for file in files:
    if (file.endswith(".png") or file.startswith("MAVRIC")) and os.path.isfile(file): #or file.startswith(".nfs"):
        os.unlink(file)

seed=0
def run_command(x):
    if x>0:
        shutil.copyfile(str(current_path) + '/MAVRIC_SECTION_0.inp', str(current_path) + '/MAVRIC_SECTION_' + str(x) + '.inp')
        with open('MAVRIC_SECTION_' + str(x) + '.inp', 'r', encoding='utf-8') as file:
            data = file.readlines()

        seed_loc = 0
        bounds_loc = 0

        for line_num, line in enumerate(data):
            if "randomSeed" in line:
                seed_loc = line_num
            if "energyBounds 1" in line:
                bounds_loc = line_num + 1

        data[seed_loc] = "    randomseed=" + str(x+seed) + "\n"
        if x>=half_procs:
            data[bounds_loc]="        linear 638 3e6 1.63955e6 \n"
        with open('MAVRIC_SECTION_' + str(x) + '.inp','w',encoding='utf-8') as file:
            file.writelines(data)
    command = 'scalerte "' + str(current_path) + '/MAVRIC_SECTION_' + str(x) + '.inp"'
    os.system(command)




# Calculate rolling averages with a window of 5
def rolling_average(data, window=5):
    pad_width = window // 2
    padded_data = np.pad(data, pad_width, mode='reflect')
    rolling_avg = np.convolve(padded_data, np.ones(window)/window, mode='valid')
    return rolling_avg

def linear_space_string(start, end, num_points):
    # Generate linearly spaced numbers
    numbers = np.linspace(start, end, num_points)
    
    # Convert each number to a string with desired format (optional: format to a certain precision)
    number_strings = [f"{num:.2f}" for num in numbers]  # Format to 2 decimal places
    
    # Join the strings with a space
    result_string = " ".join(number_strings)
    
    return result_string

if os.path.exists("ORIGEN_SECTION.f71"):
  os.remove("ORIGEN_SECTION.f71")

for x in range(num_loops):
    if x<current_loop:
        continue
    seed = x
    print("loop" + str(x+1))
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
        cuboid 10   -10  ''' + str(-10+column_length+0.5) + " " + str(column_thickness/2) + " " + str(-column_thickness/2)  + " " + str(column_thickness/2) + " " + str(-column_thickness/2) + '''
        cuboid 20   -10  ''' + str(-10+column_length) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        cuboid 30   -10   ''' + str(-10+column_length+0.5) + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + " " + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness) + '''
        ycylinder 40    0.635 '''  + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + ''' origin x=-15
        ycylinder 60    0.510 '''  + str(column_thickness/2 + shield_thickness) + " " + str(-column_thickness/2 - shield_thickness)  + ''' origin x=-15
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
        '''
            
    t = Template('''filename="$filename"''')
    concentration_file = str(current_path) + '/ORIGEN_SECTION.f71'
    filename = t.substitute({'filename' : concentration_file})
     
    MAVRIC_final  =     '''
        parameters ''' + str(detector_loc) + ''' 5 end 
    end distribution

'Create 1400 Discrete Energy bins from 10 keV to 3000 keV
    energyBounds 1
        linear 762 1.63741e6 1e4 
    end energyBounds

'The results should be indicative of a 1200 second COunt
    timeBounds 1
        title="linear command"
        linear 1 0.0 ''' + str(looptime) + '''
    end timeBounds
end definitions

read importancemap
    gridgeometryid=100
    adjointFluxes="''' + str(current_path) + '/Flux_Weighting.adjoint.dff"' + '''
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
    randomSeed=''' + str(seed) + '''
    photons
    fissionMult=0  secondaryMult=0 


    perBatch=''' + str(perBatch) + ''' batches=''' + str(batches) + '''
end parameters



end data
end'''
    
    MAVRIC_input = "".join([MAVRIC_input, filename,MAVRIC_final])
    
    m = open("MAVRIC_SECTION_0.inp", "w")
    m.write(MAVRIC_input)
    
    m.close()
    
    
    
    
    ORIGEN_input = '''=origen

options{
    print_xs=no             %print cross sections
    digits=6                %digits=6 is high-precision
    fixed_fission_energy=no %set to yes to use 200 MeV/fission
}

bounds{gamma = [1400I 3e6 1e4]
        } %1000 Linearly spaced bins}

case(irrad){
    title="Single fluid FLiBe depletion calculation"

    lib{
       file="''' + xs_lib + '''
       pos=1
    }
    
    mat{
       %29 kg of FLiBe with 4 9s Lithium in the whole loop
       units=GRAMS
       '''
    if os.path.exists("ORIGEN_SECTION.f71"):
        iso = 'load{ file="' + str(current_path) + '/ORIGEN_SECTION_PREVIOUS.f71" pos=5 }'
        print("Used existing Concentration File")
    elif fuel_salt == 0: #initial Concentration
        iso = ''' iso=[f=12.19	li=2.24	be=1.44	%FliBe
            ca=0.00153784890292959 k=0.00104094968077379	ba=0.00104094961215541	na=0.000895869565217391	 %Trace Elements
            sr=0.000281818074443671 cl=0.000110623576123982	al=8.48718535469107E-05	rb=6.70995448937407E-06	
            cr=5.9120133212123E-06 mn=4.89645308924485E-06	ni=3.36223117739475E-06	zn=2.01661332279501E-06	
            cs=1.55235697940503E-06 sc=2.75652173913043E-07	mo=6.70995378547183E-08	co=6.70995372148073E-08]
        '''
        print("Created Concentration File")
    elif fuel_salt == 1: #initial Concentration
        iso = ''' iso=[f=760 li7=140 li6=0.01 be=90 %FLiBe
             u238=95 u235=5 %Uranium
             na=0.056 al=0.0053 cl=0.0069 k=0.065 ca=0.096 %Trace Elements
             sc=0.0000172 cr=0.00037 mn=0.000306 co=0.0000042 
             ni=0.00021 zn=0.000126 rb=0.00042 sr=0.0176 
             mo=0.0000042 cs=0.000097 ba=0.065]
        '''
        print("Created Concentration File")

    ORIGEN_final = '''
    feed=[Cr=4.13861E-10 Mn=2.60644E-10 Fe=2.3775E-10 Sr=2.3775E-11
          Na=-1.23278E-10 Al=-2.20139E-11 K=-4.84306E-11 Sc=-3.12597E-13 Ni=-1.585E-11 Mo=-3.12597E-13 Cs=-4.71097E-13 Ba=-1.98125E-10] %g/s

    }

    time{
        units=SECONDS 
        t=[ ''' + str(rad_time) + " " + str(rad_time*40) + " " + str(rad_time*1160) + " " +str(looptime) + ''' ]
    }

    flux= [ 1e13 0 0 0] %1E13 n/cm^2-s
    % change cutoff to absolute curies by element,
% in step of interest (7), but print GRAMS
    print{
        cutoff_step = 1 % default -1 for average
        rel_cutoff = no % default is yes for cutoff in percent
        % only print above 1e-3 curies
        cutoffs[ GRAMS=1e-20 ]   % default is 1e-6 percent
        nuc{
         total=yes
         units=GRAMS
        }
    }
    gamma{brem_medium=H2O}     
    save = yes %only save begin and end
} %end case



end'''
    
    
    ORIGEN_input = "".join([ORIGEN_input, iso,ORIGEN_final])
    
    
    
    
    o = open("ORIGEN_SECTION.inp", "w")
    o.write(ORIGEN_input)
    
    o.close()
    
    

    # If the output file exists, rename it.
    if os.path.exists("ORIGEN_SECTION.f71"):
      os.replace("ORIGEN_SECTION.f71", "ORIGEN_SECTION_PREVIOUS.f71")
#    
#    


    os.system('scalerte "' + str(current_path) + '/ORIGEN_SECTION.inp"')
    
    
    
    
    
    
    
    # print(MAVRIC_input)
    time_to_wait = 1000
    time_counter = 0
    while not os.path.exists("ORIGEN_SECTION.f71"):
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            print("Cocentration File Not Created")
            break
    # print(time_counter)
    if time_counter < time_to_wait:
        print("ORIGEN Run Successfully")
    
    # If file exists, delete it.
    
    # if os.path.exists("MAVRIC_SECTION.pd15.txt"):
    #   for file in files:
    #     if file.endswith(".pd15.txt"):
    #         os.remove(file)
    #   print("Deleted Old Spectra")

    if x%skip_loops == 0: 
        if os.path.exists("ORIGEN_SECTION.f71"):
            print("MAVRIC Running")
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # Launch all the tasks and wait for them to complete
                executor.map(run_command, range(num_procs))
            
        
            
        time_counter = 0
        while 1==1:
            files = os.listdir(current_path)
            num_files=0
            for file in files:
                if file.endswith(".pd15.txt"):
                    num_files+=1
            time.sleep(1)
            print(num_files)
            time_counter += 1
            if time_counter > time_to_wait:
                print("Flux Output Not Created")
                break
            if num_files==num_procs:
                break   
        

        energy = []
        flux = []
        uncertainty = []
        bad = []
        first_fail=0
        second_fail=0
        if time_counter < time_to_wait:
            print("MAVRIC Run Successfully")
            for z in range(num_procs): #loop through all monte carlo simulation results
                bins=0
                found=0
                if z==half_procs:
                    bad=[]
                    split=len(flux)
                with open("MAVRIC_SECTION_" + str(z) + ".pd15.txt", "r") as pd:
                    # Read all lines of the file into a list
                    values = pd.readlines()
                    if (len(values) > 150):
                        for lines in values:
                            if found == 2:
                                first_fail=0
                                second_fail=0
                                break
                            if lines.find("--") != -1:
                                found += 1
                                continue
                            if found == 1: #read the lines containing flux data
                                try:
                                    if bins in bad:
                                        #print("skipping index " + str(bins))
                                        continue
                                    if z==half_procs or second_fail==1: #append the higher energy data to the lists
                                        flux.append(float(lines[66:77])/(num_procs-half_procs))
                                        energy.append(float(lines[17:28])/1000)
                                        uncertainty.append((float(lines[79:90])/(num_procs-half_procs))**2)
                                    elif z>half_procs: #add the averaged values to the proper indeces
                                        #print(bad)
                                        #print(len(flux))
                                        #print(split+bins)
                                        flux[split+bins]+=float(lines[66:77])/(num_procs-half_procs)
                                        uncertainty[split+bins]+=(float(lines[79:90])/(num_procs-half_procs))**2
                                        bins+=1
                                    elif z==0 or first_fail==1: #append the lower energy data to the list
                                        flux.append(float(lines[66:77])/half_procs) 
                                        #print(lines[11:14]) 
                                        #print(len(flux))                               
                                        energy.append(float(lines[17:28])/1000)
                                        uncertainty.append((float(lines[79:90])/half_procs)**2)
                                    elif z<half_procs: #add the averaged values to the proper indeces
                                        #print(bad)
                                        #print(len(flux))
                                        #print(bins)
                                        flux[bins]+=float(lines[66:77])/half_procs
                                        #print("MAVRIC_SECTION_" + str(z))
                                        #print(lines)
                                        #print(bins)
                                        uncertainty[bins]+=(float(lines[79:90])/half_procs)**2
                                        bins+=1
        
                                except ValueError:
                                    if z==0 or z==half_procs:
                                        bad.append(len(flux))
                                        if len(flux) > len(uncertainty):
                                            flux.pop()
                                            energy.pop()
                                    continue
                    else:
                        if z==0:
                            first_fail=1
                        if z==half_procs:
                            second_fail=1
            

            # Filter the lists
            uncertainty=np.sqrt(uncertainty)
            mask = efficiency_energy.apply(lambda x: np.any(np.isclose(x, energy, atol=1.3)))
            filtered_efficiency = efficiency[mask]
            flux=[g for _, g in sorted(zip(energy, flux))]
            uncertainty=[g for _, g in sorted(zip(energy, uncertainty))]
            counts=flux*(filtered_efficiency)*detector_area*count_time
            uncertainty=uncertainty*(filtered_efficiency)*detector_area*count_time
            energy.sort()
            filtered_energy, filtered_counts, filtered_uncertainty, filtered_flux = zip(*[(e, c, u, f) for e, c, u, f in zip(energy, counts, uncertainty, flux) if f>1E-30])#(c+u) >= 1])
            # Convert the tuples numpy arrays for plotting
            filtered_energy = np.array(filtered_energy)
            filtered_counts = np.array(filtered_counts)
            filtered_flux = np.array(filtered_flux)
            filtered_uncertainty = np.array(filtered_uncertainty)
            print("Counts per second = " + str(sum(filtered_counts)/count_time))
        

        # Plot the data
        plt.figure(figsize=(12, 6))
        # Create a mask to exclude specific x positions from the initial plot
        mask = np.ones(len(filtered_energy), dtype=bool)
        for position in element_markers.keys():
            # Exclude points within ±2 units of the position
            mask &= ~(np.abs(filtered_energy - position) <= 4)

        # Plotting the general data without special points
        plt.errorbar(filtered_energy[mask], filtered_counts[mask], yerr=filtered_uncertainty[mask], fmt='o',
                    color='gray', ecolor='lightgray', elinewidth=1, capsize=2,
                    markeredgewidth=1, markeredgecolor='gray', markerfacecolor='none')

        # Dictionary to collect legend handles and labels to avoid duplicates
        legend_handles = {}
        # Compute the rolling averages for counts
        local_avg_counts = rolling_average(filtered_counts)
        # Plot highlighted data points with error bars
        for position, (marker, label, color) in element_markers.items():
            highlight_mask = np.abs(filtered_energy - position) <= 2
            condition_mask = filtered_counts > 1.1 * local_avg_counts
            combined_mask = highlight_mask & condition_mask
            if any(combined_mask):
                line, caplines, (bars,) = plt.errorbar(filtered_energy[combined_mask], filtered_counts[combined_mask], 
                yerr=filtered_uncertainty[combined_mask],
                fmt=marker, markersize=10, color="gray", ecolor="gray",
                elinewidth=1.5, capsize=3, markeredgewidth=2,
                markeredgecolor=color, markerfacecolor='none', label=label)
                # Store the handle with the label in a dictionary to ensure uniqueness
                if label not in legend_handles:
                    legend_handles[label] = line

        # Creating legend without repeating labels
        plt.legend(legend_handles.values(), legend_handles.keys(), loc='upper left', bbox_to_anchor=(1, 1))

        time_since=rad_time*40*(detector_loc-4)**2+rad_time*1160*(detector_loc-3)**2+looptime*x #time since loop startup in seconds
        plt.title("loop " + str(x+1) + " (" + str(int((time_since/3600))) + " Hours, "+ str(int((time_since/60)%60)) + " Minutes "+ str(int(time_since%60)) + " Seconds Post Loop Startup)")
        # Set the y-axis to linear scale
        plt.yscale('log')
        plt.grid()
        plt.xlabel("Gamma Energy (keV)")
        plt.ylabel("Gamma Photon Counts")
        ax = plt.gca()
        if fuel_salt:
            ax.set_ylim([1E3, 1E10])
        else:
            ax.set_ylim([1E-1, 1E7])
        ax.set_xlim([0, 3000])
        plt.tight_layout()
        # Fill in the area under the curve
        #plt.fill_between(energy, flux, color='blue')
    
        # save plot
        
        plt.savefig("loop_" + str(x+1).zfill(4))
        plt.close('all')

        filtered_energy= filtered_energy.tolist()
        filtered_counts=filtered_counts.tolist()
        filtered_uncertainty=filtered_uncertainty.tolist()
        filtered_flux=filtered_flux.tolist()
        
        filtered_energy[0:0] = ["Loop " + str(x+1), "Gamma Energy (keV)"]
        filtered_counts[0:0] = [" ", "Detector Counts"]
        filtered_uncertainty[0:0] = [" ", "Uncertainty (Counts)"]
        filtered_flux[0:0] = [" ", "Flux (n/cm^2-s)"]
        data_vals.append(filtered_energy)
        data_vals.append(filtered_counts)
        data_vals.append(filtered_uncertainty)
        data_vals.append(filtered_flux)
        data_vals.append([" "])


        with open('ORIGEN_SECTION.inp','r',encoding='utf-8') as file:
            data = file.readlines()

        time_loc = 0
        flux_loc = 0
        for line_num, line in enumerate(data):
            if "t=[" in line:
                time_loc=line_num
            if "flux= [" in line:
                flux_loc=line_num
       
        inner_time=str(rad_time) + " " + str(looptime) + " " 
        data[time_loc]="        t=[" + inner_time*(skip_loops-1) + "]\n"
        inner_flux="1e13 0 "
        data[flux_loc]="   flux= ["+inner_flux*(skip_loops-1) +"]\n"
        with open('ORIGEN_SECTION.inp','w',encoding='utf-8') as file:
            file.writelines(data)
        os.system('scalerte "' + str(current_path) + '/ORIGEN_SECTION.inp"')
        current_loop=x+(skip_loops)
        os.unlink(str(current_path) + '/ORIGEN_SECTION.inp')