
# Enter your Python code here

import pya
import math
import csv
import os

# --- Editable Variables ---

# Path for the lookup table (MODIFY THIS PATH FOR YOUR SYSTEM)
# The CSV should have columns: diagonal_length_um, width_um, phase_deg
dataPath = r"/Users/hrishitdas/KLayout/Reflectarray/LUT.csv"

# Output GDS file path (MODIFY THIS PATH)
output_gds_path = os.path.expanduser("/Users/hrishitdas/KLayout/Reflectarray/Reflectarray_KLayout.gds")

# Frequency and Wavelength
frequency_THz = 30.0
lambda0_um = 300000.0 / frequency_THz # Speed of light (um/s) / frequency (Hz)
k0 = 2 * math.pi / lambda0_um

# Array Dimensions (User defined)
# User requested 1000 x 1000 array instances.
N_instances = 10 
M_instances = 10

# Beam-steering angles (MODIFY THESE INPUTS)
# e.g., theta_b = 0 for broadside
theta_b_deg = 30  # [degrees]
phi_b_deg = 30    # [degrees]

theta_b = math.radians(theta_b_deg)
phi_b = math.radians(phi_b_deg)

# --- KLayout Setup ---

layout = pya.Layout()
# Set Database Unit (DBU) to 1 nanometer (0.001 um) for precision
layout.dbu = 0.001 
top_cell = layout.create_cell("TOP_ARRAY")
layer1 = layout.layer(1, 0) # Layer for the metal rings

# Use DBU (double precision) for all geometry definitions
dbu = layout.dbu 

# --- Helper Functions ---

def create_hexagonal_ring_cell(lyt, diagonal_length_um, width_um, layer_index):
    """Creates a single hexagonal ring cell (unit cell master)."""
    cell_name = f"Ring_D{diagonal_length_um:.2f}_W{width_um:.2f}um"
    uc_cell = lyt.create_cell(cell_name)
    
    # Calculate inner diagonal length
    inner_diag_length_um = (2 / math.sqrt(3)) * ((math.sqrt(3) / 2) * diagonal_length_um - 2 * width_um)
    
    # Define points for outer and inner hexagons (using DPoint for double precision)
    num_sides = 6
    outer_points = []
    inner_points = []
    
    for i in range(num_sides):
        angle = 2 * math.pi * i / num_sides
        outer_points.append(pya.DPoint(math.cos(angle) * diagonal_length_um / 2, 
                                       math.sin(angle) * diagonal_length_um / 2))
        inner_points.append(pya.DPoint(math.cos(angle) * inner_diag_length_um / 2, 
                                       math.sin(angle) * inner_diag_length_um / 2))

    outer_poly = pya.DPolygon(outer_points)
    inner_poly = pya.DPolygon(inner_points)
    
    # Perform boolean subtraction using the Region class for precise shapes
    outer_region = pya.DRC.drc_info_get_layout_region(lyt, uc_cell, layer_index)
    outer_region.insert(outer_poly)

    inner_region = pya.DRC.drc_info_get_layout_region(lyt, uc_cell, layer_index)
    inner_region.insert(inner_poly)
    
    # Subtract inner from outer
    ring_region = outer_region - inner_region
    
    # Insert resulting shapes into the cell
    uc_cell.shapes(layer_index).insert(ring_region)
    
    return uc_cell

def get_ring_params_from_table(idealPhase, phase_vect, diagonal_vect, width_vect):
    """Performs lookup in the CSV data to find matching parameters."""
    # Assumes table is sorted by increasing phase
    matchIndex = 0
    for i in range(len(phase_vect)):
        if phase_vect[i] >= idealPhase:
            matchIndex = i 
            break
    
    return diagonal_vect[matchIndex], width_vect[matchIndex]

# --- Main Generation Logic ---

print("Loading lookup table...")

diagonal_vect = []
width_vect = []
phase_vect = []

try:
    # Use standard Python open() as pya environment supports it
    with open(dataPath, 'r') as csvfile: 
        reader = csv.reader(csvfile, delimiter=',')
        next(reader) # Skip header row
        for row in reader:
            diagonal_vect.append(float(row[2])) # [um] 3rd column in csv
            width_vect.append(float(row[1]))   # [um] 2nd column in csv
            phase_vect.append(float(row[0]))   # [deg] 1st column in csv
except Exception as e:
    print(f"Error reading CSV file: {e}")
    # Stop script execution if lookup fails
    raise SystemExit(f"Could not read lookup table at {dataPath}")

print(f"Loaded {len(phase_vect)} entries from lookup table.")

# Unit cell periods for hexagonal grid (based on your input scaling)
# Assuming a nominal size of lambda0/4 from your script
uc_p_um = 5   # Period x [um]
uc_q_um = uc_p_um / math.cos(math.pi/6) # Period y [um]

F_um = 100000.0 # Focal length scaled to micrometers (1m)
x0, y0 = 0.0, 0.0 # Horn center

# Create a dictionary to store unique unit cells to reuse them
created_cells = {}

print("Generating array instances...")

# Iterate over N x M instances to calculate required phase at each location
for n in range(N_instances):
    for m in range(M_instances):
        # Calculate center coordinates for the hexagonal grid
        centerX = (n - (N_instances - 1) / 2) * uc_p_um + (m % 2) / 2 * uc_p_um
        # centerY calculation from HFSS script seems different for true hex grid Y offset
        # A common q calculation is based on uc_p/sqrt(3) vertical step
        centerY = (m - (M_instances - 1) / 2) * uc_q_um * 0.75 # Adjusted Y spacing approximation from HFSS script

        # Calculate the required ideal phase at this location for beam steering
        d = math.sqrt((x0 - centerX)**2 + (y0 - centerY)**2 + F_um**2) 
        idealPhase_rad = k0 * (d - (centerX * math.cos(phi_b) + centerY * math.sin(phi_b)) * math.sin(theta_b))
        idealPhase_deg = math.degrees(idealPhase_rad) % 360 

        # Look up ring dimensions based on the required phase
        ring_d, ring_w = get_ring_params_from_table(idealPhase_deg, phase_vect, diagonal_vect, width_vect)
        
        # Use cell hierarchy: create a master cell for this specific ring type if it doesn't exist
        cell_key = f"{ring_d:.2f}_{ring_w:.2f}"
        if cell_key not in created_cells:
            uc_master_cell = create_hexagonal_ring_cell(layout, ring_d, ring_w, layer1)
            created_cells[cell_key] = uc_master_cell
        else:
            uc_master_cell = created_cells[cell_key]
        
        # Place an instance of the master cell at the correct location
        # Use DTrans (double precision transformation) for placement
        trans = pya.DTrans(centerX, centerY)
        top_cell.insert(pya.DCellInst(uc_master_cell.cell_index(), trans))


print(f"Total unique cell masters created: {len(created_cells)}")
print(f"Total instances placed: {N_instances * M_instances}")
print(f"Writing GDS file to: {output_gds_path}")

# Write the layout to GDS file using the specified absolute path
layout.write(output_gds_path)
print("GDS file generation complete.")
