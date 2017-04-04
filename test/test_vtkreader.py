import sys
sys.path.append('../vector_field_plot_tool')
from vtk_reader import VTKReader

vr = VTKReader('sphere000000.vtu')
vr.extract_data()

# Now we have self.vector_field and self.coordinates defined

a, b, c = vr.interpolate_data(-3, 3, -3, 3,
                              nx=10, ny=10,
                              _filter=None)

# print(c)

from plot_vector_field import plot_vector_field
import matplotlib
matplotlib.pyplot.rcParams['backend'] = 'TkAgg'

plot_vector_field(vr, -3, 3, -3, 3, colorbar=True, savefig='test.pdf')
