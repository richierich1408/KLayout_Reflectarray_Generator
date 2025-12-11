# Enter your Python code here
import os
import pya

layout = pya.Layout()
top = layout.create_cell("TOP")
l1 = layout.layer(1, 0)
top.shapes(l1).insert(pya.Box(0, 0, 1000, 2000))

# Define the output path in your user's home directory (e.g., on the Desktop)
output_directory = os.path.expanduser("~/Desktop")
output_file = os.path.join(output_directory, "t.gds")

# Write the layout to the specified full path
layout.write(output_file)

# Optional: Print the location to the KLayout Console for verification
print(f"GDS file saved to: {output_file}")