import numpy as np
import scipy.interpolate


# -----------------------------------------------------------------------------

class VectorFieldReader(object):
    """
    Base class for the readers
    """

    def __init__(self):

        self.vector_field = None
        self.coordinates = None

    def interpolate_data(self,
                         x_min, x_max, y_min, y_max,
                         nx=20, ny=20,
                         interpolator='scipy',
                         interpolator_method='linear',
                         _filter=None
                         ):
        """

        x_min, y_min, etc   :: Range of interpolation for the vector field

        """

        x, y = self.coordinates[:, 0], self.coordinates[:, 1]

        if _filter is None:
            _filter = x < 2 * np.max(x)

        # Coordinates for the discretised quiver/field plot
        xi = np.linspace(x_min, x_max, nx)
        yi = np.linspace(y_min, y_max, ny)

        interp_data = [None] * 3

        if interpolator == 'scipy':

            xi, yi = np.meshgrid(xi, yi)

            # Fill values are NaN to avoid plotting them
            for i in range(3):
                interp_data[i] = scipy.interpolate.griddata(
                    (x[_filter], y[_filter]),
                    self.vector_field[:, i][_filter],
                    (xi, yi),
                    method='linear',
                    fill_value=np.nan
                )

        return xi, yi, interp_data
