import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from a345_utilities import print_header, mem_usage     
from matplotlib.patches import Rectangle as rect 
import numpy as np
import os, time
from photutils import DAOStarFinder
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
from a345_utilities import print_header 
from photutils import DAOStarFinder
from astropy.stats import mad_std
from astropy.coordinates import SkyCoord
from photutils import CircularAnnulus
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from astropy.time import Time
import re
import warnings
warnings.filterwarnings('ignore')
from astropy.stats import sigma_clipped_stats
import pandas as pd



from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError




RADIUS = 19

def cal_star_mags(star:str,band:str,exposure:str, radius: int = RADIUS, solver: bool = True, Threshold: int = 4) -> None:
    path_cal = '/Volumes/external_2T/calibration/2023-10/neg10c/master'
    path_data =  '/Volumes/external_2T'
    target_jd = []
    target_mags = []
    calibration_mags = []

    '''
    Function to calibrate raw .fits files, plate solve and extract magnitudes \n
    Star is a string \n
    Band is the colour band \n 
    Exposure is the exposure time in seconds (time with s)
    '''
    
    data_cal = np.transpose(np.loadtxt(path_data + '/'+ 'cal_stars/' + star + '_calibration_stars.txt', skiprows=1, delimiter=","))

    mag_g_cal=data_cal[3]
    mag_g_err = data_cal[4]
    mag_i_cal=data_cal[7]
    mag_i_err = data_cal[8]
    mag_r_cal=data_cal[5]
    mag_r_err = data_cal[6]
    ra_cal=data_cal[1]
    dec_cal=data_cal[2]
    
    if solver == True:
        ast = AstrometryNet()
        ast.API_URL = 'http://nova.astro.gla.ac.uk/api' # local server
        ast.api_key = 'XXXXXXXX'
        ast.URL = 'http://nova.astro.gla.ac.uk'
    else:
        ast = AstrometryNet()
        ast.API_URL = 'http://nova.astrometry.net/api' # local server
        ast.URL = 'http://nova.astrometry.net'
        #API key for Giles Hammond
        ast.api_key = 'cyilczrjxdbmnhum'
            
    counter = 0
    with fits.open(path_cal+'/dark_flat_gr_3s_master.fits') as hdu:
        flatdark_data = hdu[0].data
    
    with fits.open(path_cal+'/dark_'+exposure+'_master.fits') as hdu:
        dark_data = hdu[0].data

    if re.search('G',band):
        mag_cal1 = mag_g_cal
        mag_err = mag_g_err
        with fits.open(path_cal+'/flat_g_master.fits') as hdu:
            flat = hdu[0].data      
    
    if re.search('I',band):
        mag_cal1 = mag_i_cal
        mag_err = mag_i_err
        
        with fits.open(path_cal+'/flat_i_master.fits') as hdu:
            flat = hdu[0].data 
    
    if re.search('R',band):
        mag_cal1 = mag_r_cal
        mag_err = mag_r_err
        
        with fits.open(path_cal+'/flat_r_master.fits') as hdu:
            flat = hdu[0].data 
    
    star_list = os.listdir(path_data+'/'+star+'/'+band)       

    for img in star_list:
        print(img)

            
        if img.endswith('.fits'):          
            with fits.open(path_data+'/'+star+'/'+band+'/'+img) as hdu:
                img_header = hdu[0].header
                img_data = hdu[0].data
            if re.search('_E_',img):
                img_data = np.rot90(img_data, k=2)

            img_c = (img_data-dark_data)*(np.mean(flat-flatdark_data))/(flat-flatdark_data)
            mean, median, std = sigma_clipped_stats(img_c)
            daofind = DAOStarFinder(fwhm=4, threshold = Threshold*std)  
            sources = daofind(img_c)

            for col in sources.colnames:  
                sources[col].info.format = '%.8g'
            sources.sort('flux')
            sources.reverse()
            if counter == 0:
                print(len(sources))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++PLATE SOLVING++++++++++++++++++++++++++++++++++++++++++++++++++++
            try:
                wcs_header = None 
                wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'],
                                                        img_c.shape[1], 
                                                        img_c.shape[0],
                                                        #scale_est =  0.724,
                                                        #scale_units = 'arcsecperpix',
                                                        solve_timeout=300)
                if wcs_header:
                    print('Success')
                else:
                    print('Failed') 
                
                img_header.update(wcs_header)
                if counter == 0:
                    
                    print(img_header['OBJCTDEC'],img_header['OBJCTRA'])
                hdu = fits.PrimaryHDU()
                hdu.header.update(img_header)       
                hdu.data = img_c
                # hdu.writeto(filename+'corrected'+img, overwrite=True)
                wcs = WCS(img_header)
            
            
            except TimeoutError:       
                print('\n -> ##FAIL: Timeout while solving, try a longer timeout, optmise number of sources (200-800 seems about right)')
                wcs = WCS(img_header)
            # with fits.open(filename+'corrected'+img) as hdu:
            #     reloaded_header = hdu[0].header
            #     reloaded_data   = hdu[0].data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++PLATE SOLVING++++++++++++++++++++++++++++++++++++++++++++++++++++
            # number of sources from DAO Starfinder to plot
            N_source=20
            # number of calibration stars to plot




            # plot yellow circles around the sources found by DAO starfinder
            positions_dao = np.transpose((sources['xcentroid'], sources['ycentroid']))  
            apertures_define = CircularAperture(positions_dao, r=25)
            apertures_dao=apertures_define[0:N_source+1]
            # apertures_dao.plot(color='yellow', lw=1, alpha=1)

            r1=radius
            r2=r1+2
            r3=r2+4
            
            
            
            # plot blue circles around the calibration stars from VizieR
            for i in range(len(ra_cal)):
                
                mag_cal1 = np.sort(np.array(mag_cal1))
                sourcex , sourcey = wcs.wcs_world2pix(ra_cal[i],dec_cal[i],1)
                source = np.transpose((sourcex, sourcey))
                source_aperture = CircularAperture(source, r1)
                source_annulus = CircularAnnulus(source, r2, r3)
                source_phot = [source_aperture, source_annulus]
                # source_aperture.plot(color='blue', lw=2, alpha=1) 
                # source_annulus.plot(color='deepskyblue', lw=2, alpha=1)
                phot_table_source = aperture_photometry(img_c, source_phot)
                for col in phot_table_source.colnames:
                    phot_table_source[col].info.format = '%.8g'  # for consistent table output
                bkg_mean_cal = float(phot_table_source[0]['aperture_sum_1'] / source_annulus.area)
                bcal = bkg_mean_cal * source_aperture.area
                cal_flux=float(phot_table_source1[0]['aperture_sum_0'] - bcal)
                mag_cal_arr = []
                for j in range(len(ra_cal)):
                    if i == j:
                        continue
                    sourcex , sourcey = wcs.wcs_world2pix(ra_cal[i],dec_cal[i],1)
                    source = np.transpose((sourcex, sourcey))
                    source_aperture = CircularAperture(source, r1)
                    source_annulus = CircularAnnulus(source, r2, r3)
                    source_phot = [source_aperture, source_annulus]
                    # source_aperture.plot(color='blue', lw=2, alpha=1) 
                    # source_annulus.plot(color='deepskyblue', lw=2, alpha=1)
                    phot_table_source = aperture_photometry(img_c, source_phot)
                    for col in phot_table_source.colnames:
                        phot_table_source[col].info.format = '%.8g'  # for consistent table output
                    bkg_mean_cal = float(phot_table_source[0]['aperture_sum_1'] / source_annulus.area)
                    bcal = bkg_mean_cal * source_aperture.area
                    cal_flux2=float(phot_table_source1[0]['aperture_sum_0'] - bcal)
                    mag_cal=mag_cal1[j] + 2.5*np.log10(cal_flux/cal_flux2)

                    mag_cal_arr.append(mag_cal)
                    
                mag_cal_final = np.mean(mag_cal_arr)
                calibration_mags.append(mag_cal_final)
                

            
            
            
            source1_x, source1_y= wcs.wcs_world2pix(ra_cal,dec_cal,1)
            source1 = np.transpose((source1_x, source1_y))
            source1_aperture = CircularAperture(source1, r1)  
            source1_annulus = CircularAnnulus(source1, r2, r3)

            source1_phot = [source1_aperture, source1_annulus]
            # source1_aperture.plot(color='blue', lw=2, alpha=1)
            # source1_annulus.plot(color='deepskyblue', lw=2, alpha=1)


            source2_x, source2_y =(sources['xcentroid'] , sources['ycentroid'] )

            Coor_x, Coor_y = wcs.wcs_pix2world(source2_x, source2_y, 1) 
            source2=np.transpose((source2_x, source2_y))
            source2_aperture = CircularAperture(source2, r1)  
            # source2_aperture.plot(color='red', lw=2, alpha=1)

            source2 = np.transpose((source2_x, source2_y))
            source2_aperture = CircularAperture(source2, r1)  
            source2_annulus = CircularAnnulus(source2, r2, r3)

            source2_phot = [source2_aperture, source2_annulus]
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++APERTURE PHOTOMETRY+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            phot_table_source1 = aperture_photometry(img_c, source1_phot)
            phot_table_source2 = aperture_photometry(img_c, source2_phot)

            for col in phot_table_source1.colnames:
                phot_table_source1[col].info.format = '%.8g'  # for consistent table output

            bkg_mean_cal1 = float(phot_table_source1[0]['aperture_sum_1'] / source1_annulus.area)
            bcal1 = bkg_mean_cal1 * source1_aperture.area

            cal1_flux=float(phot_table_source1[0]['aperture_sum_0'] - bcal1)


            for col in phot_table_source2.colnames:
                phot_table_source2[col].info.format = '%.8g'  # for consistent table output

            bkg_mean_targ = float(phot_table_source2[0]['aperture_sum_1'] / source2_annulus.area)

            targcal = bkg_mean_targ * source2_aperture.area

            r = np.linspace(1, 30, 30)
            test_r2 = r+2
            test_r3 = r+4
            if counter == 0:
                flux = []
                plt.figure()
                for i in np.arange(1,31,1):
                    
                    source2 = np.transpose((source2_x, source2_y))
                    source2_aperture = CircularAperture(source2, i)  
                    source2_annulus = CircularAnnulus(source2, test_r2[i-1], test_r3[i-1])

                    source2_phot = [source2_aperture, source2_annulus]
                    phot_table_source2 = aperture_photometry(img_c, source2_phot)

                    bkg_mean_targ = float(phot_table_source2[0]['aperture_sum_1'] / source2_annulus.area)

                    targcal = bkg_mean_targ * source2_aperture.area
                    targ_flux=float(phot_table_source2[0]['aperture_sum_0'] - targcal)
                    flux.append(targ_flux)
                    
                plt.plot(r,flux)
                plt.show()
            
            targ_flux=float(phot_table_source2[0]['aperture_sum_0'] - targcal)
            mag_targ=mag_cal1[0] + 2.5*np.log10(cal1_flux/targ_flux)

            t_fits=img_header['DATE-OBS']
            t = Time(t_fits, format='isot', scale='utc')
            t_jd=t.jd 
            target_mags.append(mag_targ)
            target_jd.append(t_jd)
            print(mag_targ)
            
            if counter == 0:
                break
            counter +=1
    g_inds = mag_g_cal.argsort()
    g_err = mag_g_err[g_inds]
    g_cal = mag_g_cal[g_inds]
    i_inds = mag_i_cal.argsort()
    i_err = mag_i_err[i_inds]
    i_cal = mag_i_cal[i_inds]
    r_inds = mag_r_cal.argsort()
    r_err = mag_r_err[r_inds]
    r_cal = mag_r_cal[r_inds]
    
    data = np.transpose(np.array([calibration_mags, g_cal, g_err, i_cal, i_err, r_cal, r_err]))

    np.savetxt(star + 'err.txt',data,delimiter=',',header='mags cali stars, g mag, g err, i mag, i err, r mag, r err')


    
    return calibration_mags, g_cal, g_err, i_cal, i_err, r_cal, r_err


