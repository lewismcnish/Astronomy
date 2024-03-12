import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import pandas as pd


def plot(juliandate: np.ndarray, magnitude: np.ndarray, mag_err: np.ndarray) -> None:
    plt.style.use('report.mplstyle')
    plt.figure(figsize=(10, 5))
    plt.scatter(juliandate, magnitude )
    plt.errorbar(juliandate, magnitude, yerr=mag_err, fmt='o')
    plt.xlabel('Time (JD)')
    plt.ylabel('Magnitude')
    plt.title('Light Curve')
    plt.show()    

def plot_lc(file: str) -> None:
    data = pd.read_csv(file)
    jd = data['JD']
    mag = data['Magnitude']
    try:
        mag_err = data['Mag_err']
        
    except:
        mag_err = None
    plot(jd, mag, mag_err)
