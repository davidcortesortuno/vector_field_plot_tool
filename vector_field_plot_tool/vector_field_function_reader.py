from vector_field_reader import VectorFieldReader
import numpy as np


# -----------------------------------------------------------------------------


class VectorFieldFunctionReader(VectorFieldReader):
    """
    Reader for a function that returns a vector field as a function
    of coordinates
    """

    def __init__(self, vf_function):

        super(VectorFieldReader, self).__init__()

        self.vf_function = vf_function
