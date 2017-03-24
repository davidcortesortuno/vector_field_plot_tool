import numpy as np


# -----------------------------------------------------------------------------

class VectorFieldReader(object):
    """
    Base class for the readers
    """

    def __init__(self):

        self.field_data = None
