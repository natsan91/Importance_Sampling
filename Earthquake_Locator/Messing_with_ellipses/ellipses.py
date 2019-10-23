from __future__ import division
import pylab
import numpy

import os
os.environ["PROJ_LIB"] = "C:\\Users\\natha\\Anaconda3\\Library\\share"

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

class Basemap(Basemap):
    def ellipse(self, x0, y0, a, b, n, ax=None, **kwargs):
        """
        Draws a polygon centered at ``x0, y0``. The polygon approximates an
        ellipse on the surface of the Earth with semi-major-axis ``a`` and 
        semi-minor axis ``b`` degrees longitude and latitude, made up of 
        ``n`` vertices.

        For a description of the properties of ellipsis, please refer to [1].

        The polygon is based upon code written do plot Tissot's indicatrix
        found on the matplotlib mailing list at [2].

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.patches.Polygon

        RETURNS
            poly : a maptplotlib.patches.Polygon object.

        REFERENCES
            [1] : http://en.wikipedia.org/wiki/Ellipse


        """
        ax = kwargs.pop('ax', None) or self._check_ax()
        g = pyproj.Geod(a=self.rmajor, b=self.rminor)
        # Gets forward and back azimuths, plus distances between initial
        # points (x0, y0)
        azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0+a, x0], [y0, y0+b])
        tsid = dist[0] * dist[1] # a * b

        # Initializes list of segments, calculates \del azimuth, and goes on 
        # for every vertex
        seg = [self(x0+a, y0)]
        AZ = numpy.linspace(azf[0], 360. + azf[0], n)
        for i, az in enumerate(AZ):
            # Skips segments along equator (Geod can't handle equatorial arcs).
            if numpy.allclose(0., y0) and (numpy.allclose(90., az) or
                numpy.allclose(270., az)):
                continue

            # In polar coordinates, with the origin at the center of the 
            # ellipse and with the angular coordinate ``az`` measured from the
            # major axis, the ellipse's equation  is [1]:
            #
            #                           a * b
            # r(az) = ------------------------------------------
            #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
            #
            # Azymuth angle in radial coordinates and corrected for reference
            # angle.
            azr = 2. * numpy.pi / 360. * (az + 90.)
            A = dist[0] * numpy.sin(azr)
            B = dist[1] * numpy.cos(azr)
            r = tsid / (B**2. + A**2.)**0.5
            lon, lat, azb = g.fwd(x0, y0, az, r)
            x, y = self(lon, lat)

            # Add segment if it is in the map projection region.
            if x < 1e20 and y < 1e20:
                seg.append((x, y))

        poly = Polygon(seg, **kwargs)
        ax.add_patch(poly)

        # Set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

        return poly

pylab.close('all')
pylab.ion()

m = Basemap(width=12000000, height=8000000, resolution='l', projection='stere',
            lat_ts=50, lat_0=50, lon_0=-107.)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(numpy.arange(-80.,81.,20.))
m.drawmeridians(numpy.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
# draw ellipses
ax = pylab.gca()
for y in numpy.linspace(m.ymax/20, 19*m.ymax/20, 9):
    for x in numpy.linspace(m.xmax/20, 19*m.xmax/20, 12):
        lon, lat = m(x, y, inverse=True)
        poly = m.ellipse(lon, lat, 3, 1.5, 100, facecolor='green', zorder=10,
            alpha=0.5)
pylab.title("Ellipses on stereographic projection")
#input()