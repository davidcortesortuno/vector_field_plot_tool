import sys
sys.path.append('../vector_field_plot_tool')
from vtk_reader import VTKReader
from plot_vector_field import plot_vector_field, plot_scalar_field
import matplotlib.pyplot as plt

vr = VTKReader('sphere000000.vtu')
vr.extract_data()
# Now we have self.vector_field and self.coordinates defined

# We can interpolate some data:
vr.interpolator_method = 'cubic'
a, b, c = vr.interpolate_data(-3, 3, -3, 3,
                              nx=10, ny=10,
                              _filter=None)

# print(c)

# import matplotlib
# matplotlib.pyplot.rcParams['backend'] = 'TkAgg'

ax = plot_scalar_field(vr, -3, 3, -3, 3, hsv_map=True)
plot_vector_field(vr, -3, 3, -3, 3, ax=ax)
plt.savefig('test.pdf', bbox_inches='tight')
