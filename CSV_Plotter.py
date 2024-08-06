import pathlib
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import imageio
inlet=3 #State 3 in the concentration file which is the inlet to the medical room/pump
outlet=4 #State 4 in the concentration file which is the outlet from the medical room/pump

current_path = pathlib.Path().resolve()
files = os.listdir(current_path)

for file in files:
    if (file.endswith(".png") or file.startswith("MAVRIC")) and os.path.isfile(file): 
        os.unlink(file)

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
column_length=240
shield_thickness=1.2
plot_flux=1

if detector_loc==inlet:
    loc = "inlet"
elif detector_loc==outlet:
    loc = "outlet"

if fuel_salt:
    salt_type="Fuel_Salt"
else:
    salt_type="Clean_Salt"

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

current_path = pathlib.Path().resolve()

output_path = str(current_path) + '/Clean_Salt_inlet_10cm_per_s__11loops_240cm_column_length'     # Path where the video will be saved
with open(output_path + '.csv', 'r') as file:
    lines = file.readlines()

def rolling_average(data, window=5):
    pad_width = window // 2
    padded_data = np.pad(data, pad_width, mode='reflect')
    rolling_avg = np.convolve(padded_data, np.ones(window)/window, mode='valid')
    return rolling_avg

x=0
y=0
for row in lines:
    y+=1
    if row.find("Loop") != -1:
        x = int(row[4:row.find(',')])
        print(x)
    if row.find("Gamma") != -1:
        # Split the string by commas and strip any whitespace
        components = [item.strip() for item in row.split(',')]
        # Filter out the non-numeric components and convert to floats
        gamma = np.array([float(item) for item in components if item.replace('.', '', 1).isdigit()])
    if row.find("Detector") != -1:
        # Split the string by commas and strip any whitespace
        components = [item.strip() for item in row.split(',')]
        # Filter out the non-numeric components and convert to floats
        # Filter out the non-numeric components and convert to floats
        filtered_energy=[]
        filtered_counts=[]
        for energy, counts in zip(gamma,components[2:]):
            try:
                # Try to convert the flux value to a float
                counts_float = float(counts)
                # If successful, add both energy and flux values to the filtered lists
                filtered_energy.append(float(energy))
                filtered_counts.append(counts_float)
            except ValueError:
                # If conversion fails, print a message and skip the value
                print(f"Skipping invalid flux value: {counts}")
    if row.find("Uncertainty") != -1:
        # Split the string by commas and strip any whitespace
        components = [item.strip() for item in row.split(',')]
        # Filter out the non-numeric components and convert to floats
        filtered_energy=[]
        filtered_uncertainty=[]
        for energy, uncertainty in zip(gamma,components[2:]):
            try:
                # Try to convert the flux value to a float
                uncertainty_float = float(uncertainty)
                # If successful, add both energy and flux values to the filtered lists
                filtered_energy.append(float(energy))
                filtered_uncertainty.append(uncertainty_float)
            except ValueError:
                # If conversion fails, print a message and skip the value
                print(f"Skipping invalid flux value: {uncertainty}")
        if plot_flux == 0:
            filtered_energy = np.array(filtered_energy)
            filtered_counts = np.array(filtered_counts)
            filtered_uncertainty = np.array(filtered_uncertainty)
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

            time_since=rad_time*40*(detector_loc-4)**2+rad_time*1160*(detector_loc-3)**2+looptime*(x-1) #time since loop startup in seconds
            plt.title("loop " + str(x) + " (" + str(int((time_since/3600))) + " Hours, "+ str(int((time_since/60)%60)) + " Minutes "+ str(int(time_since%60)) + " Seconds Post Loop Startup)")
            # Set the y-axis to linear scale
            plt.yscale('log')
            plt.grid()
            plt.xlabel("Gamma Energy (keV)")
            plt.ylabel("Gamma Photon Counts")
            
            ax = plt.gca()
            if fuel_salt:
                ax.set_ylim([1E3, 1E10])
            else:
                ax.set_ylim([1E-1, 1E8])
            ax.set_xlim([0, 3000])
            plt.tight_layout()
            # Fill in the area under the curve
            #plt.fill_between(energy, flux, color='blue')
        
            # save plot
            
            plt.savefig("loop_" + str(x).zfill(4))
            plt.close('all')
    if row.find("Flux") != -1:
        # Split the string by commas and strip any whitespace
        components = [item.strip() for item in row.split(',')]
        # Filter out the non-numeric components and convert to floats
        #filtered_flux = np.array([float(item) for item in components if item.replace('.', '', 1).isdigit()])
        filtered_energy=[]
        filtered_flux=[]
        for energy, flux in zip(gamma,components[2:]):
            try:
                # Try to convert the flux value to a float
                flux_float = float(flux)
                # If successful, add both energy and flux values to the filtered lists
                filtered_energy.append(float(energy))
                filtered_flux.append(flux_float)
            except ValueError:
                # If conversion fails, print a message and skip the value
                print(f"Skipping invalid flux value: {flux}")
        if plot_flux:
            # Plot the data
            plt.figure(figsize=(12, 6))
            plt.plot(filtered_energy,filtered_flux)

            time_since=rad_time*40*(detector_loc-4)**2+rad_time*1160*(detector_loc-3)**2+looptime*(x-1) #time since loop startup in seconds
            plt.title("loop " + str(x) + " (" + str(int((time_since/3600))) + " Hours, "+ str(int((time_since/60)%60)) + " Minutes "+ str(int(time_since%60)) + " Seconds Post Loop Startup)")
            # Set the y-axis to linear scale
            plt.yscale('log')
            plt.grid()
            plt.xlabel("Gamma Energy (keV)")
            plt.ylabel("Gamma Flux (n/cm^2-s)")
            ax = plt.gca()
            if fuel_salt:
                ax.set_ylim([1E3, 1E10])
            elif plot_flux:
                ax.set_ylim([1E-20, 1E3])
            else:
                ax.set_ylim([1E-1, 1E8])
            ax.set_xlim([0, 3000])

            plt.tight_layout()
            # Fill in the area under the curve
            #plt.fill_between(energy, flux, color='blue')
        
            # save plot
            
            plt.savefig("loop_" + str(x).zfill(4))
            plt.close('all')
    
        

image_files = [os.path.join(str(current_path), img) for img in sorted(os.listdir(str(current_path))) if img.endswith(".png")]

fps = 6  # Frames per second
# Create a video writer object
writer = imageio.get_writer(output_path + '.mp4', fps=fps)

# Append images to the video writer
for image_file in image_files:
    image = imageio.imread(image_file)
    writer.append_data(image)
    print(f"Added {image_file} to video")

writer.close()  # Close the writer to finalize the video
print("Video created successfully.")