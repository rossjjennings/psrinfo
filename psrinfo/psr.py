import numpy as np
from numpy import pi, sin, cos, exp, log, sqrt
import astropy.units as u
from astropy.units import pc, yr, km, s, deg, radian, mas, hourangle
from astropy.coordinates import Angle, Galactic, ICRS, BarycentricTrueEcliptic
from functools import wraps

def cached(func):
    name = func.__name__
    
    @wraps(func)
    def wrapper(self):
        try:
            return self._cache[name]
        except KeyError:
            self._cache[name] = func(self)
            return self._cache[name]
        except AttributeError:
            self._cache = {name: func(self)}
            return self._cache[name]
    
    return wrapper

class pulsar:
    def __init__(self, name, coords, **other_data):
        self.name = name
        self.coords = coords
        self.__dict__.update(other_data)
    
    @staticmethod
    def from_equatorial(rec):
        other_data = dict()
        try:
            try:
                coords = ICRS(ra = Angle(rec.ra, unit=hourangle),
                              dec = Angle(rec.dec, unit=deg),
                              pm_ra_cosdec = float(rec.pmra) * mas/yr,
                              pm_dec = float(rec.pmdec) * mas/yr)
                coord_cov = np.array([[rec.ra_err**2*(180/12)**2, 0],
                                      [0, rec.dec_err**2]])
                pm_cov = np.array([[rec.pmra_err**2, 0],
                                   [0, rec.pmdec_err**2]])
                other_data.update(coord_cov=coord_cov, pm_cov=pm_cov)
            except:
                coords = ICRS(ra = Angle(rec.ra, unit=hourangle),
                              dec = Angle(rec.dec, unit=deg))
                coord_cov = np.array([[rec.ra_err**2*(180/12)**2, 0],
                                      [0, rec.dec_err**2]])
                other_data.update(coord_cov=coord_cov)
        except AttributeError:
            try:
                coords = ICRS(ra = Angle(rec.lon, unit=hourangle),
                              dec = Angle(rec.lat, unit=deg),
                              pm_ra_cosdec = float(rec.pmlon) * mas/yr,
                              pm_dec = float(rec.pmlat) * mas/yr)
                coord_cov = np.array([[rec.lon_err**2*(180/12)**2, 0],
                                      [0, rec.lat_err**2]])
                pm_cov = np.array([[rec.pmlon_err**2, 0],
                                   [0, rec.pmlat_err**2]])
                other_data.update(coord_cov=coord_cov, pm_cov=pm_cov)
            except ValueError:
                coords = ICRS(ra = Angle(rec.lon, unit=hourangle),
                              dec = Angle(rec.lat, unit=deg))
                coord_cov = np.array([[rec.lon_err**2*(180/12)**2, 0],
                                      [0, rec.lat_err**2]])
                other_data.update(coord_cov=coord_cov)
        
        exclude = ['name', 'ra', 'dec', 'pmra', 'pmdec',
                   'lat', 'lon', 'pmlat', 'pmlon']
        other_data.update({name: item for name, item in rec._asdict().items()
                           if name not in exclude})
        return pulsar(rec.name, coords, **other_data)
    
    @staticmethod
    def from_ecliptic(rec):
        other_data = dict()
        try:
            try:
                coords = BarycentricTrueEcliptic(lon = Angle(rec.elon, unit=deg),
                                                 lat = Angle(rec.elat, unit=deg),
                                                 pm_lon_coslat = float(rec.pmelon) * mas/yr,
                                                 pm_lat = float(rec.pmelat) * mas/yr)
                coord_cov = np.array([[rec.elon_err**2, 0],
                                      [0, rec.elat_err**2]])*3600**2
                pm_cov = np.array([[rec.pmelon_err**2, 0],
                                   [0, rec.pmelat_err**2]])
                other_data.update(coord_cov=coord_cov, pm_cov=pm_cov)
            except:
                coords = BarycentricTrueEcliptic(lon = Angle(rec.elon, unit=deg),
                                                 lat = Angle(rec.elat, unit=deg))
                coord_cov = np.array([[rec.elon_err**2, 0],
                                      [0, rec.elat_err**2]])*3600**2
                other_data.update(coord_cov=coord_cov)
        except AttributeError:
            try:
                coords = BarycentricTrueEcliptic(lon = Angle(rec.lon, unit=deg),
                                                 lat = Angle(rec.lat, unit=deg),
                                                 pm_lon_coslat = float(rec.pmlon) * mas/yr,
                                                 pm_lat = float(rec.pmlat) * mas/yr)
                coord_cov = np.array([[rec.lon_err**2, 0],
                                      [0, rec.lat_err**2]])*3600**2
                pm_cov = np.array([[rec.pmlon_err**2, 0],
                                   [0, rec.pmlat_err**2]])
                other_data.update(coord_cov=coord_cov, pm_cov=pm_cov)
            except ValueError:
                coords = BarycentricTrueEcliptic(lon = Angle(rec.lon, unit=deg),
                                                 lat = Angle(rec.lat, unit=deg))
                coord_cov = np.array([[rec.lon_err**2, 0],
                                      [0, rec.lat_err**2]])*3600**2
                other_data.update(coord_cov=coord_cov)
        
        exclude = ['name', 'elon', 'elat', 'pmelon', 'pmelat',
                   'lat', 'lon', 'pmlat', 'pmlon']
        other_data.update({name: item for name, item in rec._asdict().items()
                           if name not in exclude})
        return pulsar(rec.name, coords, **other_data)
    
    @cached
    def equatorial(self):
        return self.coords.transform_to(ICRS)
    
    @cached
    def ecliptic(self):
        return self.coords.transform_to(BarycentricTrueEcliptic)
    
    @cached
    def galactic(self):
        return self.coords.transform_to(Galactic)
    
    @cached
    def galactic_rad(self):
        l = self.galactic().l.to(radian).value
        b = self.galactic().b.to(radian).value
        return np.array([l, b])
    
    @cached
    def pm_equatorial(self):
        pm_ra = self.equatorial().pm_ra_cosdec
        pm_dec = self.equatorial().pm_dec
        return np.array([pm_ra.value, pm_dec.value])
        
    @cached
    def pm_ecliptic(self):
        pm_lon = self.ecliptic().pm_lon_coslat
        pm_lat = self.ecliptic().pm_lat
        return np.array([pm_lon.value, pm_lat.value])
    
    @cached
    def pm_galactic(self):
        pm_l = self.galactic().pm_l_cosb
        pm_b = self.galactic().pm_b
        return np.array([pm_l.value, pm_b.value])
    
    @cached
    def equatorial_to_ecliptic(self):
        c_ra = ICRS(self.equatorial().ra, self.equatorial().dec,
                    pm_ra_cosdec = 1 * mas/yr, pm_dec = 0 * mas/yr)
        c_dec = ICRS(self.equatorial().ra, self.equatorial().dec,
                    pm_ra_cosdec = 0 * mas/yr, pm_dec = 1 * mas/yr)
        e_ra = c_ra.transform_to(BarycentricTrueEcliptic)
        e_dec = c_dec.transform_to(BarycentricTrueEcliptic)
        return np.array([[e_ra.pm_lon_coslat.value, e_dec.pm_lon_coslat.value],
                         [e_ra.pm_lat.value, e_dec.pm_lat.value]])
    
    @cached
    def equatorial_to_galactic(self):
        c_ra = ICRS(self.equatorial().ra, self.equatorial().dec,
                    pm_ra_cosdec = 1 * mas/yr, pm_dec = 0 * mas/yr)
        c_dec = ICRS(self.equatorial().ra, self.equatorial().dec,
                    pm_ra_cosdec = 0 * mas/yr, pm_dec = 1 * mas/yr)
        g_ra = c_ra.transform_to(Galactic)
        g_dec = c_dec.transform_to(Galactic)
        return np.array([[g_ra.pm_l_cosb.value, g_dec.pm_l_cosb.value],
                         [g_ra.pm_b.value, g_dec.pm_b.value]])
    
    @cached
    def ecliptic_to_galactic(self):
        e_lon = BarycentricTrueEcliptic(self.ecliptic().lon, self.ecliptic().lat,
                    pm_lon_coslat = 1 * mas/yr, pm_lat = 0 * mas/yr)
        e_lat = BarycentricTrueEcliptic(self.ecliptic().lon, self.ecliptic().lat,
                    pm_lon_coslat = 0 * mas/yr, pm_lat = 1 * mas/yr)
        g_lon = e_lon.transform_to(Galactic)
        g_lat = e_lat.transform_to(Galactic)
        return np.array([[g_lon.pm_l_cosb.value, g_lat.pm_l_cosb.value],
                         [g_lon.pm_b.value, g_lat.pm_b.value]])
    
    @cached
    def equatorial_cov(self):
        if type(self.coords) is ICRS:
            return self.coord_cov
        elif type(self.coords) is BarycentricTrueEcliptic:
            transform = self.equatorial_to_ecliptic()
            return np.matmul(transform, np.matmul(self.coord_cov, transform.T))
    
    @cached
    def ecliptic_cov(self):
        if type(self.coords) is ICRS:
            transform = np.linalg.inv(self.equatorial_to_ecliptic())
            return np.matmul(transform, np.matmul(self.coord_cov, transform.T))
        elif type(self.coords) is BarycentricTrueEcliptic:
            return self.coord_cov
    
    @cached
    def galactic_cov(self):
        if type(self.coords) is ICRS:
            transform = self.equatorial_to_galactic()
            return np.matmul(transform, np.matmul(self.coord_cov, transform.T))
        elif type(self.coords) is BarycentricTrueEcliptic:
            transform = self.ecliptic_to_galactic()
            return np.matmul(transform, np.matmul(self.coord_cov, transform.T))
    
    @cached
    def pm_equatorial_cov(self):
        if type(self.coords) is ICRS:
            return self.pm_cov
        elif type(self.coords) is BarycentricTrueEcliptic:
            transform = self.equatorial_to_ecliptic()
            return np.matmul(transform, np.matmul(self.pm_cov, transform.T))
    
    @cached
    def pm_ecliptic_cov(self):
        if type(self.coords) is ICRS:
            transform = np.linalg.inv(self.equatorial_to_ecliptic())
            return np.matmul(transform, np.matmul(self.pm_cov, transform.T))
        elif type(self.coords) is BarycentricTrueEcliptic:
            return self.pm_cov
    
    @cached
    def pm_galactic_cov(self):
        if type(self.coords) is ICRS:
            transform = self.equatorial_to_galactic()
            return np.matmul(transform, np.matmul(self.pm_cov, transform.T))
        elif type(self.coords) is BarycentricTrueEcliptic:
            transform = self.ecliptic_to_galactic()
            return np.matmul(transform, np.matmul(self.pm_cov, transform.T))
    
    def set_pm_equatorial(self, pm_ra_cosdec, pm_dec, pm_ra_cosdec_err,
                          pm_dec_err, pm_ra_dec_cov = 0):
        self.coords = ICRS(self.equatorial().ra, self.equatorial().dec,
                           pm_ra_cosdec = pm_ra_cosdec * mas/yr, pm_dec = pm_dec * mas/yr)
        self.pm_cov = np.array([[pm_ra_cosdec_err**2, pm_ra_dec_cov],
                                [pm_ra_dec_cov, pm_dec_err**2]])
        del self._cache
   
    def set_pm_ecliptic(self, pm_lon_coslat, pm_lat, pm_lon_coslat_err,
                        pm_lat_err, pm_lon_lat_cov = 0):
        self.coords = BarycentricTrueEcliptic(self.ecliptic().lon, self.ecliptic().lat,
                                              pm_lon_coslat = pm_lon_coslat * mas/yr,
                                              pm_lat = pm_lat * mas/yr)
        self.pm_cov = np.array([[pm_lon_coslat_err**2, pm_lon_lat_cov],
                                [pm_lon_lat_cov, pm_lat_err**2]])
        del self._cache
    
    def set_pm_galactic(self, pm_l_cosb, pm_b, pm_l_cosb_err,
                        pm_b_err, pm_l_b_cov = 0):
        self.coords = Galactic(self.galactic().l, self.galactic().b,
                               pm_l_cosb = pm_l_cosb * mas/yr, pm_b = pm_b * mas/yr)
        self.pm_cov = np.array([[pm_l_cosb_err**2, pm_l_b_cov],
                                [pm_l_b_cov, pm_b_err**2]])
        del self._cache
    
    def __repr__(self):
        coordtype = repr(self.coords).split()[0][1:]
        return "<PSR {} ({})>".format(self.name, coordtype)
