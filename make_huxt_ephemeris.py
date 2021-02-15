# Standard library
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
import h5py
# Our own library for using spice with STEREO (https://github.com/LukeBarnard/stereo_spice)
from stereo_spice.coordinates import StereoSpice
spice = StereoSpice()

def make_ephemeris():
    ephem = h5py.File('ephemeris.hdf5', 'w')
    stereo_begin = Time('2007-01-01T00:00:00')

    bodies = ['EARTH', 'VENUS', 'MERCURY', 'STA', 'STB']
    for body in bodies:

        body_group = ephem.create_group(body)

        if (body == 'STA') | (body=='STB'):
            t_start = '2007-01-01T00:00:00'
        else:
            t_start = '1963-01-01T00:00:00'

        t_stop = '2029-01-01T00:00:00'
        times = pd.date_range(t_start, t_stop, freq='4H')
        times = Time(times.to_pydatetime(), format='datetime')

        for system in ['CARR', 'HEEQ', 'HAE']:

            coord_group = body_group.create_group(system)

            coords = np.zeros((times.size, 3))*np.NaN

            for i, time in enumerate(times):

                if (body in ['STA', 'STB']) & (time < stereo_begin):
                    continue

                coords[i, :] = spice.get_lonlat(time, body, system, degrees=True)

            coord_group.create_dataset('time', data=times.jd)
            rad = coord_group.create_dataset('radius', data=coords[:, 0])
            rad.attrs['unit'] = u.km.to_string()
            lon = coord_group.create_dataset('longitude', data=coords[:, 1])
            lon.attrs['unit'] = u.deg.to_string()
            lat = coord_group.create_dataset('latitude', data=coords[:, 2])
            lat.attrs['unit'] = u.deg.to_string()

            ephem.flush()

    ephem.close()
    return

if __name__=="__main__":
    make_ephemeris()