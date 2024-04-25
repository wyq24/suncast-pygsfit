from astropy import constants as const, units as u
import numpy as np
from astropy.coordinates import SkyCoord

def extent_convertor(cheader):
    # Replace direct assignments with lookups
    CRVAL1 = cheader['CRVAL1']
    CDELT1 = cheader['CDELT1']
    CRPIX1 = cheader['CRPIX1']
    NAXIS1 = cheader['NAXIS1']

    CRVAL2 = cheader['CRVAL2']  # Assuming similar keys for Y-axis
    CDELT2 = cheader['CDELT2']
    CRPIX2 = cheader['CRPIX2']
    NAXIS2 = cheader['NAXIS2']

    # Calculate extent
    x_extent_left = CRVAL1 - (CRPIX1 - 1) * CDELT1
    x_extent_right = x_extent_left + NAXIS1 * CDELT1
    y_extent_bottom = CRVAL2 - (CRPIX2 - 1) * CDELT2
    y_extent_top = y_extent_bottom + NAXIS2 * CDELT2
    extent = [x_extent_left, x_extent_right, y_extent_bottom, y_extent_top]
    return extent



def sfu2tb_2d(frequency, flux_2d, area=None, square=True, reverse=False, verbose=False):
    """
    Convert a 2D array of flux densities (in sfu) to brightness temperatures (in K), or vice versa,
    across a given frequency, considering the source's area or assuming a square/elliptical shape.

    Parameters:
        frequency: Frequency in Hz (assumed same for the entire 2D array)
        flux_2d: 2D array of flux densities in sfu (if not reverse) or brightness temperatures in K (if reverse)
        area: Source area in arcsec^2 (if known)
        square: Assume square-shaped source if True; elliptical if False (ignored if area is specified)
        reverse: Convert brightness temperature to flux density if True
        verbose: Print detailed operations if True
    """

    c = const.c.cgs
    k_B = const.k_B.cgs
    sfu = u.Jy * 1e4  # sfu to Jansky conversion

    # Frequency handling
    frequency = frequency * u.Hz if not hasattr(frequency, 'unit') else frequency

    # Flux handling
    flux_unit = u.K if reverse else sfu
    flux_2d = flux_2d * flux_unit if not hasattr(flux_2d, 'unit') else flux_2d

    if area is None:
        raise ValueError('Area must be provided or calculable from source size.')

    area = area * (u.arcsec ** 2) if not hasattr(area, 'unit') else area

    sr = area.to(u.radian ** 2)
    factor = c ** 2 / (2 * k_B * frequency ** 2 * sr)

    if reverse:
        if verbose: print('Converting input brightness temperatures in K to flux densities in sfu for each point.')
        result = (flux_2d / factor).to(sfu, equivalencies=u.dimensionless_angles())
    else:
        if verbose: print('Converting input flux densities in sfu to brightness temperatures in K for each point.')
        result = (flux_2d * factor).to(u.K, equivalencies=u.dimensionless_angles())

    return result


def world2pix(xy_coord, refmap):
    solar_x = xy_coord[0] * u.arcsec
    solar_y = xy_coord[1] * u.arcsec
    coord = SkyCoord(solar_x, solar_y, frame=refmap.coordinate_frame)
    return refmap.world_to_pixel(coord)






