
# Enter your Python code here

import pya
import os # Import the os module

# Create a new layout
layout = pya.Layout()
top = layout.create_cell("Gratings_Flat")
layer1 = layout.layer(1, 0) # GDS Layer 1, Datatype 0

# Grating parameters
bar_width = 0.5 # in layout units (microns if not specified)
pitch = 1.0     # in layout units (pitch = bar_width + spacing)
grating_length = 10.0
num_bars = 50

# Generate the bars
for i in range(num_bars):
    # Calculate the x-coordinate of the current bar
    x_start = i * pitch
    x_end = x_start + bar_width
    
    # Create a box (rectangle) for the bar
    # pya.DBox uses double precision coordinates
    box = pya.DBox(x_start, 0, x_end, grating_length)
    
    # Insert the box into the top cell on the specified layer
    top.shapes(layer1).insert(box)

# Define the directory path and filename
directory = "/Users/hrishitdas/KLayout/Example Python Scripts"
filename = "grating.gds"

# Use os.path.join to create the full, correct path
full_path = os.path.join(directory, filename)

# Write the layout using the combined path
layout.write(full_path)

print(f"File saved to: {full_path}")

# Optional: Load the generated layout in the KLayout GUI if running as a macro
#if pya.Application.instance():
#    pya.MainWindow.instance().load_layout("grating.gds", 1)
    

