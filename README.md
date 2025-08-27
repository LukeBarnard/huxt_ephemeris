# huxt_ephemeris

This script produces the ephemeris file used by HUXt. It produces a HDF5 file of the coordinates of Mercury, Venus,
Earth, Mars, Jupiter, and Saturn, as well as  STEREO-A and STEREO-B, for the period 1963-01-01 until 2029-01-01, at 
4 hour resolution. Coordintaes are provided in the Carrington, HEEQ, and HAE systems. Coordinates for STEREO-A and
STEREO-B after 2020 are based on a predicted ephemeris, rather than the definitive ephemeris.

This script depends on the [stereo_spice](https://github.com/LukeBarnard/stereo_spice) package. 

