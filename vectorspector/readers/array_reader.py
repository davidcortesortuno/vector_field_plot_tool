from ..vector_field_reader import VectorFieldReader
import numpy as np


# -----------------------------------------------------------------------------


class ArrayReader(VectorFieldReader):
    """
    Reader for data obtained directly form two arrays with:
    coordinates, vector field data
    """

    def __init__(self, coordinates, vector_field):

        super().__init__()

        self.vector_field = vector_field
        self.coordinates = coordinates
