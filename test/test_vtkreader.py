import sys
sys.path.append('../vector_field_plot_tool')
from vtk_reader import VTKReader

vr = VTKReader('sphere000000.vtu')
vr.extract_data()

# Now we have self.vector_field and self.coordinates defined
