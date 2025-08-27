# Standard library
from astropy.time import Time
import astropy.units as u
import astropy.coordinates as acoords
import h5py
import sunpy.coordinates as coords
import numpy as np


def get_naif_body_codes_dict():
    """ Return a dictionary with the names and naif codes of bodies that can be looked up in Horizons for use in
    generating an offline ephemeris for HUXt."""

    bodies = {'PSP': -96,
              'SOLO': -144,
              'STA': -234,
              'STB': -235,
              'MERCURY': 199,
              'VENUS': 299,
              'EARTH': 399,
              'MARS': 499,
              'JUPITER': 599,
              'SATURN': 699}
    return bodies


def zerototwopi(angles):
    """
    Function to constrain angles to the 0 - 2pi domain.
    Args:
        angles: a numpy array of angles
    Returns:
        angles_out: a numpy array of angles in the 0 - 2pi domain.
    """
    twopi = 2.0 * np.pi
    angles_out = angles
    a = -np.floor_divide(angles_out, twopi)
    angles_out = angles_out + (a * twopi)
    return angles_out


def make_ephemeris():
    ephem = h5py.File('ephemeris_new.hdf5', 'w')

    bodies_dict = get_naif_body_codes_dict()
    for body, naif_code in bodies_dict.items():

        body_group = ephem.create_group(body)

        # Get start and stop times for the ephemeris, based on the body
        t_start = Time('1963-01-01T00:00:00').to_datetime()
        t_stop = Time('2029-01-01T00:00:00').to_datetime()
        if body in ['STA', 'STB']:
            t_start = Time('2007-01-01T00:00:00')
            if body == 'STA':
                t_stop = Time('2025-11-10T00:00:00').to_datetime()
            elif body == 'STB':
                t_stop = Time('2024-10-25T00:00:00').to_datetime()
        elif body == 'PSP':
            t_start = Time('2018-08-13T00:00:00')
        elif body == 'SOLO':
            t_start = Time('2020-02-11T00:00:00')

        if body in ['STA', 'STB', 'PSP', 'SOLO']:
            time_lookup = {'start': t_start, 'stop': t_stop, 'step': '6H'}
        else:
            time_lookup = {'start': t_start, 'stop': t_stop, 'step': '24H'}

        body_coords = coords.get_horizons_coord(naif_code, time_lookup)

        for coord_sys in ['CARR', 'HEEQ', 'HAE']:
            coord_group = body_group.create_group(coord_sys)
            coord_group.create_dataset('time', data=body_coords.obstime.jd)

            if coord_sys == 'CARR':
                this_coord = body_coords.transform_to(coords.HeliographicCarrington(observer="self"))
            elif coord_sys == 'HEEQ':
                this_coord = body_coords.transform_to(coords.HeliographicStonyhurst())
            elif coord_sys == 'HAE':
                this_coord = body_coords.transform_to(acoords.HeliocentricMeanEcliptic())

            if coord_sys == 'HAE':
                rad = coord_group.create_dataset('radius', data=this_coord.distance.to(u.km).value)
            else:
                rad = coord_group.create_dataset('radius', data=this_coord.radius.to(u.km).value)

            rad.attrs['unit'] = u.km.to_string()
            lon = coord_group.create_dataset('longitude', data=np.rad2deg(zerototwopi(this_coord.lon.radian)))
            lon.attrs['unit'] = u.deg.to_string()
            lat = coord_group.create_dataset('latitude', data=np.rad2deg(this_coord.lat.radian))
            lat.attrs['unit'] = u.deg.to_string()

            ephem.flush()

    ephem.close()
    return


if __name__ == "__main__":

    make_ephemeris()
