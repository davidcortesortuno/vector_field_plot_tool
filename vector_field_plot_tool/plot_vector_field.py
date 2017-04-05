import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


# -----------------------------------------------------------------------------


def plot_vector_field(data_reader,  # one of our reader classes
                      x_min, x_max,
                      y_min, y_max,
                      vf_component='z',
                      normalize_data=True,
                      figsize=(8., 8.),
                      cmap='viridis',
                      hsv_map=False,
                      cmap_alpha=1.,
                      quiver_map=None,
                      quiver_type='interpolated_cmap',
                      quiver_color='k',
                      pivot='mid',
                      nx=20,
                      ny=20,
                      ax=None,
                      clim=[-1, 1],
                      frame=True,
                      **quiver_args
                      ):

    """
    Make a 2D plot from the data extracted of the VTK file, using
    a colormap with interpolated values.
    IT IS NECESSARY to run the extract_data() function before.
    If a new range of z_values is required, simply reassign the self.z_min
    and self.z_max attributes
    When setting quiver_type as interpolated, the numbers of arrows can be
    controled specifying the nx and ny parameeters, which are the
    number of entities along x and y respectively.
    OPTIONS:
    x_min, x_max        :: Range of spatial x values to be used in the 2D
                           plot to interpolate the data for the colormap
    y_min, y_max        :: Range of spatial y values to be used in the 2D
                           plot to interpolate the data for the colormap
    vf_component        :: Component of the vector field that is going
                           to be shown as the magnitude of every entity
                           in a colormap. By default, it is plotted
                           the z component magnitude of the vectors.
                           Options:
                                'x', 'y', 'z'
    normalize_data      :: Set False if the colorbar ticks values are in the
                           range of the real data. By default, the colormap is
                           normalised from -1 to 1
    nx, ny              :: Resolution in the x and y directions
                           for the interpolations using the data points,
                           i.e. the number of divisions
                           between x_min and x_max; y_min and y_max
    figsize             :: Dimensions of the plot as a tuple,
                           (8, 8) by default
    cmap                :: Palette of the colourmap (not considered if
                           using the hsv_map option)
    hsv_map             :: With this option the colormap is going to use
                           the HSV palette, where the x and y
                           components of the
                           vectors are mapped into the Hue values and the z
                           component is done in the S and V, so that the
                           maximum z values are shown in white and the
                           minimum in black. For 2 dimensional vector
                           fields, this makes the colormap to be only black
                           (since all the z components are set to zero),
                           thus this option can be passed as:
                                '2d' or '3d'.
                           The 2d option set S and V as 1, so the plot
                           shows the full color. The 3d option makes the
                           mapping black or white according to the z
                           component of the field.
                           This mapping is useful for showing a vector
                           field without a quiver plot.
    cmap_alpha          :: Transparency value of the colourmap
    quiver_map          :: Colour palette of the arrows of the vector
                           field. By default it is the inverted
                           palette of cmap
    quiver_type         :: By default the quiver plot is not interpolated,
                           it shows all the data points in the specified
                           spatial ranges (raw data), and it is shown with
                           a colormap.  This option lets the user choose to
                           interpolate the vector field and if a colormap
                      cmap=plt.get_cmap(cmap),
                           or a single color is used. The options are:
                                'interpolated_cmap', 'interpolated_color',
                                'raw_cmap', 'raw_color'
    quiver_color        :: Arrow color if one of the 'color' options was
                           specified in the quiver_type argument
    pivot               :: By default we make the arrows to be drawn at the
                           center of the grid nodes. This option is from
                           the matplotlib quiver function
    nx, ny              :: Resolution in the x and y directions for the
                           arrows in the quiver plot if one of the
                           interpolated quiver_type options are passed
                           (number of divisions between x_min and x_max;
                           y_min and y_max). By default: 20 x 20 arrows are
                           drawn
    frame               :: Frame of the plot
    ax                  :: Can be a predefined matplotlib axis object to
                           show the plot on it. This is useful to make
                           a grid of plots
    **quiver_args       :: Any extra keyword arguments for the quiver plot

    # DEPRECATED::
    interpolator_hsv_method     :: Method for scipy, for the HSV mapping.
                                   Default: 'linear'
    interpolator_quiver_method :: Method for scipy or natgrid when
                                  interpolating the quiver plot, default:
                                  'linear' or 'nn'
    TODO:
        Add polar components
        Add titles
    """

    # Vector field components
    cs = {'x': 0, 'y': 1, 'z': 2}

    # Extract data from the data reader
    # _filter = np.logical_and(data_reader.coordinates[:, 2] > 2.5,
    #                          data_reader.coordinates[:, 2] < 4)

    (xi, yi,
     quiv_xyz) = data_reader.interpolate_data(x_min, x_max,
                                              y_min, y_max,
                                              nx=nx, ny=ny,
                                              )

    # ---------------------------------------------------------------------
    # Now plot in matplotlib ----------------------------------------------
    # ---------------------------------------------------------------------

    # Use a predefined axis if possible
    if not ax:
        fig = plt.figure(figsize=figsize, frameon=frame)
        ax = fig.add_subplot(111)

    if not hsv_map:
        pass
        # Plot the colour map with the interpolated values of v_i
        # ax.pcolormesh(xi, yi, quiv_xyz[cs[vf_component]],
        #               cmap=plt.get_cmap(cmap),
        #               vmin=-1, vmax=1,
        #               alpha=cmap_alpha)
    else:
        pass
        # Plot the colour map with the HSV colours
        # ax.imshow(zi, interpolation='None',
        #           extent=[np.min(xi), np.max(xi),
        #                   np.min(yi), np.max(yi)],
        #           vmin=-1, vmax=1,
        #           origin='lower'
        #           )

    if (quiver_type == 'interpolated_cmap'
            or quiver_type == 'raw_cmap'):
        ax.quiver(xi,
                  yi,
                  quiv_xyz[0],
                  quiv_xyz[1],
                  # paint the vectors according to the
                  # component of m_component
                  quiv_xyz[cs[vf_component]],
                  cmap=cmap,
                  pivot=pivot,
                  clim=clim,
                  **quiver_args
                  )

    elif (quiver_type == 'interpolated_color'
            or quiver_type == 'raw_colour'):
        ax.quiver(xi,
                  yi,
                  quiv_xyz[0],
                  quiv_xyz[1],
                  color=quiver_color,
                  pivot=pivot,
                  clim=clim,
                  **quiver_args
                  )
    # elif not quiver_type:
    #     pass
    else:
        print('Specify a valid option for the quiver plot')
        return

    # if not quiver_map:
    #     quiver_map = cmap + '_r'

    # # Use whole data if the vector field is not inerpolated
    # if (quiver_type == 'raw_cmap'
    #         or quiver_type == 'raw_colour'):
    #     quiv['vx'] = self.vf[:, 0][self.data_filter],
    #     quiv['vy'] = self.vf[:, 1][self.data_filter]
    #     quiv['vz'] = self.vf[:, 2][self.data_filter]

    #     xi_q, yi_q = x, y

    return ax


def plot_scalar_field(data_reader,  # one of our reader classes
                      x_min, x_max,
                      y_min, y_max,
                      vf_component='z',
                      normalize_data=True,
                      nx=20, ny=20,
                      xlim=None,
                      ylim=None,
                      figsize=(8., 8.),
                      cmap='gist_earth',
                      hsv_map=False,
                      alpha=1.,
                      frame=True,
                      ax=None,
                      clim=[-1, 1],
                      x_label=r'$x$',
                      y_label=r'$y$',
                      ):

    """
    This function plots the scalar field as a background
    """

    # Vector field components
    cs = {'x': 0, 'y': 1, 'z': 2}

    (xi, yi,
     scalar_xyz) = data_reader.interpolate_data(x_min, x_max,
                                                y_min, y_max,
                                                nx=nx, ny=ny
                                                )
    dx = (xi[0, 1] - xi[0, 0]) * 0.5
    dy = (yi[1, 0] - yi[0, 0]) * 0.5

    # ---------------------------------------------------------------------
    # Now plot in matplotlib ----------------------------------------------
    # ---------------------------------------------------------------------
    # Use a predefined axis if possible
    if not ax:
        fig = plt.figure(figsize=figsize, frameon=frame)
        ax = fig.add_subplot(111)

    if hsv_map:
        field_data = np.vstack((scalar_xyz[0],
                                scalar_xyz[1],
                                scalar_xyz[2]
                                )).reshape(3, -1).T
        rgb_data = generate_rgbs(field_data)
        rgb_data.shape = (scalar_xyz[0].shape[0], scalar_xyz[0].shape[0], 3)
        # rgb_data = rgb_data.T.reshape(3, -1)
        # grid_nx = scalar_xyz[0].shape[0]
        # rgbs = [rgb_data[0].reshape(grid_nx, -1),
        #         rgb_data[1].reshape(grid_nx, -1),
        #         rgb_data[2].reshape(grid_nx, -1)
        #         ]

        # Plot the colour map with the interpolated values of v_i
        ax.imshow(rgb_data,
                  extent=[np.min(xi) - dx, np.max(xi) + dx,
                          np.min(yi) - dy, np.max(yi) + dy],
                  vmin=-1, vmax=1,
                  alpha=alpha)
    else:
        # Plot the colour map with the HSV colours
        ax.imshow(scalar_xyz[cs[vf_component]],
                  interpolation='None',
                  extent=[np.min(xi) - dx, np.max(xi) + dx,
                          np.min(yi) - dy, np.max(yi) + dy],
                  vmin=clim[0], vmax=clim[1],
                  origin='lower',
                  alpha=alpha,
                  cmap=plt.get_cmap(cmap),
                  )

    return ax


def plot_colorbar(ax,
                  cmap='viridis',
                  hsv_map=False,
                  normalize_data=True,
                  colorbar_label='',
                  clim=[-1, 1]
                  ):
    if hsv_map:
        cmap_cb = matplotlib.cm.get_cmap(name='hsv')
    else:
        cmap_cb = matplotlib.cm.get_cmap(name=cmap)

    if normalize_data or hsv_map:
        norm = matplotlib.colors.Normalize(-1, 1)
    else:
        norm = matplotlib.colors.Normalize(vmin=clim[0],
                                           vmax=clim[1]
                                           )

    # Add axes for the colorbar with respect to the top image
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)

    # Colorbar
    cbar = matplotlib.colorbar.ColorbarBase(cax,
                                            cmap=cmap_cb,
                                            norm=norm,
                                            # ticks=[-1, 0, 1],
                                            orientation='vertical',
                                            )

    cbar.set_label(colorbar_label, rotation=270)

    # # Label HSV colorbar accordingly
    # if hsv_map:
    #     cbar.set_ticks([1, 0, -1])
    #     cbar.set_ticklabels([r'$2\pi$', r'$\pi$', r'$0$'])
    #     # cbar.update_ticks()


# TOOLS -----------------------------------------------------------------------
import colorsys


def convert_to_rgb(hls_color):
    return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                        hls_color[1],
                                        hls_color[2]
                                        ))


def generate_rgbs(field_data):
    """
    field_data      ::  (n, 3) array
    """
    hls = np.ones_like(field_data)
    hls[:, 0] = np.arctan2(field_data[:, 1],
                           field_data[:, 0]
                           )
    hls[:, 0][hls[:, 0] < 0] = hls[:, 0][hls[:, 0] < 0] + 2 * np.pi

    hls[:, 1] = 0.5 * (field_data[:, 2] + 1)
    # self.hls[:, 2] = 0.5 * (self.sim.spin.reshape(-1, 3)[:, 2] + 1)

    rgbs = np.apply_along_axis(convert_to_rgb, 1, hls)

    # Some RGB values can get very small magnitudes, like 1e-10:
    rgbs[rgbs < 0] += 2 * np.pi

    return rgbs
