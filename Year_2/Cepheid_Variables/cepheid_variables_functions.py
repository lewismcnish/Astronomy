import numpy as np     # always import numpy at the start

def powerSpectrum(relT,mag,periods):
    N = len(periods)
    power = np.zeros(N)
    for i,P in enumerate(periods):
        phi = 2.0*np.pi*relT/P
        power[i] = (np.sum(mag*np.cos(phi)))**2 + (np.sum(mag*np.sin(phi)))**2 
    return power/(N**2)
