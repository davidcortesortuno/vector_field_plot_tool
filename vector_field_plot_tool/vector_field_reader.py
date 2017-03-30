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
                         nx_q=20, ny_q=20,
                         interpolator='scipy',
                         interpolator_method='linear',
                         _filter=None
                         ):
        """

        x_min, y_min, etc   :: Range of interpolation for the vector field

        """
        if _filter is None:
            _filter = self.coordinates[:, 0] < 2 * np.max(self.coordinates[:, 0])

        x, y = self.coordinates[:, 0], self.coordinates[:, 1]

        # Coordinates for the discretised quiver/field plot
        xi = np.linspace(x_min, x_max, nx_q)
        yi = np.linspace(y_min, y_max, ny_q)

        interp_data = np.zeros((len(xi), 3))

        if interpolator == 'scipy':

            xi, yi = np.meshgrid(xi, yi)

            # Fill values are NaN to avoid plotting them
            for i in range(3):
                interp_data[:, i] = scipy.interpolate.griddata(
                    (x[_filter], y[_filter]),
                    self.vector_field[:, i][_filter],
                    (xi, yi),
                    method='linear',
                    fill_value=np.nan
                )

        return xi, yi, interp_data
