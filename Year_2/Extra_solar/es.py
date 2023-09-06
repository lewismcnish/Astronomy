import numpy as np     # always import numpy at the start
from scipy.interpolate import interp1d      # you will need this scipy function 

def shiftSpectrum(wavelength, intensity, v):
    c = 3e8
    shifted_wavelength = wavelength * (1.0 + v/c)
    shifted_intensity = intensity / (1.0 + v/c)
    interpolated_intensity = interp1d(shifted_wavelength,shifted_intensity,bounds_error=False,fill_value=(intensity[0], intensity[-1]))(wavelength)
    return interpolated_intensity
