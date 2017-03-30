from vector_field_reader import VectorFieldReader
import numpy as np


# -----------------------------------------------------------------------------


class VectorFieldFunctionReader(VectorFieldReader):
    """

    """

    def __init__(self, vf_function):

        super(VectorFieldReader, self).__init__()

        self.vf_function = vf_function
