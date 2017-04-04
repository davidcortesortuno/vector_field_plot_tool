from vector_field_reader import VectorFieldReader
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


# -----------------------------------------------------------------------------


class VTKReader(VectorFieldReader):
    """

    """

    def __init__(self, vtk_file, vtkfiletype='XMLUnstructuredGrid'):

        super().__init__()

        self.vtk_file = vtk_file
        self.vtkfiletype = vtkfiletype

        grid_reader_function = 'vtk' + vtkfiletype + 'Reader'
        self.reader_function = getattr(vtk, grid_reader_function)()

    def extract_data(self, rotate=None):
        """
        Extract the data from a VTK file using the python-vtk library
        The data is passed to a numpy array for better manipulation
        """

        # Check the type of file to load it according to its grid structure
        # which was specified when calling the class These dictionary entries
        # return the corresponding VTK Reader functions, so they can be easily
        # called according to the class argument
        #
        # Rotate argument is a rotation in the x-y plane. To use, set rotate
        # equal to an angle in radians.

        self.reader = self.reader_function
        self.reader.SetFileName(self.vtk_file)

        if not self.vtkfiletype.startswith('XML'):
            # For Non XML vtk files:
            self.reader.ReadAllVectorsOn()
            self.reader.ReadAllScalarsOn()

        self.reader.Update()

        # Get the coordinates of the nodes in the mesh
        nodes_vtk_array = self.reader.GetOutput().GetPoints().GetData()

        # Get The vector field (data of every node)
        vf_vtk_array = self.reader.GetOutput().GetPointData().GetArray(0)

        # Transform the coordinates of the nodes into a Numpy array and
        # save them to the corresponding class objects
        nodes = vtk_to_numpy(nodes_vtk_array)
        self.coordinates = nodes
        if rotate:
            self.x = nodes[:, 0] * np.cos(rotate) - nodes[:, 1] * np.sin(rotate)
            self.y = nodes[:, 0] * np.sin(rotate) + nodes[:, 1] * np.cos(rotate)
            self.z = nodes[:, 2]
        else:
            self.x, self.y, self.z = (nodes[:, 0],
                                      nodes[:, 1],
                                      nodes[:, 2]
                                      )

        # Transform the vector field data into a Numpy array and save
        if rotate:
            vf = vtk_to_numpy(vf_vtk_array)
            vfx = vf[:, 0] * np.cos(rotate) - vf[:, 1] * np.sin(rotate)
            vfy = vf[:, 0] * np.sin(rotate) + vf[:, 1] * np.cos(rotate)
            vfz = vf[:, 2]
            self.vector_field = np.zeros_like(vf)
            self.vector_field[:, 0] = vfx
            self.vector_field[:, 1] = vfy
            self.vector_field[:, 2] = vfz
        else:
            self.vector_field = vtk_to_numpy(vf_vtk_array)

