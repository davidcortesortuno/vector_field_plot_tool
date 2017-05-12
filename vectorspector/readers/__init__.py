from . import array_reader
from . import vector_field_function_reader
try:
    from . import vtk_reader
except ImportError as err:
    print('WARNING: Cannot load the VTK reader since a library was not found')
    print('Exception error message:')
    print(err)
