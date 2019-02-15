import numpy as np
from numpy import pi, sin, cos, exp, log, sqrt
import astropy.units as u
from astropy.coordinates import SkyCoord

def predict_pos(pulsar, epoch, format='mjd'):
    '''
    Predict the position of a pulsar based on its proper motion.
    
    Input
    -----
    `pulsar`: Pulsar for which to estimate distance
    `epoch`: Time at which to estimate position (integer MJD)
    
    Output
    ------
    `pos`: Predicted position (`SkyCoord` object)
    '''
    pulsar_ra = pulsar.equatorial().ra.to(u.deg).value
    pulsar_dec = pulsar.equatorial().dec.to(u.deg).value
    position_cov = pulsar.equatorial_cov()
    title = pulsar.name
    try:
        julian_years_elapsed = (epoch - pulsar.posepoch)/365.25
        ra_offset = pulsar.pm_equatorial()[0]/cos(pulsar_dec*pi/180)*julian_years_elapsed/1000
        dec_offset = pulsar.pm_equatorial()[1]*julian_years_elapsed/1000
        pm_offset_cov = pulsar.pm_equatorial_cov()*julian_years_elapsed**2/1000**2
        cos_dec = cos(pulsar_dec*pi/180)
        pm_offset_cov *= np.array([[1/cos_dec**2, 1/cos_dec], [1/cos_dec, 1]])
        position_cov += pm_offset_cov
    except TypeError:
        ra_offset = 0
        dec_offset = 0
    
    return SkyCoord(pulsar_ra + ra_offset/3600, pulsar_dec + dec_offset/3600,
                    unit='deg', frame='icrs')
