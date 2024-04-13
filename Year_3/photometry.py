import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from a345_utilities import print_header    
from matplotlib.patches import Rectangle as rect 
import numpy as np
import os, time
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
from photutils.aperture import aperture_photometry, CircularAperture
from photutils.aperture import CircularAnnulus
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')
from astropy.stats import sigma_clipped_stats
import pandas as pd
import re
from astropy.time import Time
from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError


# =================================Calibration============================================
def calibration(star, band , exposure:str = '60s', six_x_six: bool = True):
    
    '''This function calibrates the images of a star in a given band and exposure time.
    
    
    Parameters:
    -----------
    star : str
        Name of the star.
    band : str
        Band of the images.
    exposure : str
        Exposure time of the images.
    six_x_six : bool
        If True, the function will use the 6x6 images, otherwise it will use the 1x1 images.
        '''
    counter = 0
    
    if six_x_six == True:
        path_cal = '/Volumes/external_2T/6x6_cal/master'
        path_data =  '/Volumes/external_2T'
    else:
        path_cal = '/Volumes/external_2T/calibration/2023-10/neg10c/master'
        path_data =  '/Volumes/external_2T'
    
    
    with fits.open(path_cal+'/dark_flat_gr_3s_master.fits') as hdu:
        flatdark_data = hdu[0].data
    
    with fits.open(path_cal+'/dark_'+exposure+'_master.fits') as hdu:
        dark_data = hdu[0].data

    if re.search('G',band):
        with fits.open(path_cal+'/flat_g_master.fits') as hdu:
            flat = hdu[0].data      
    
    if re.search('I',band):
        with fits.open(path_cal+'/flat_i_master.fits') as hdu:
            flat = hdu[0].data 
    
    if re.search('R',band):
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
            if counter == 0:
                plt.figure(figsize=(10,10))
                plt.imshow(img_c, cmap='gray', origin='lower')
                plt.colorbar()
                plt.title('Corrected image')
                plt.show()
                counter += 1
            hdu = fits.PrimaryHDU(img_c)
            hdu.header = img_header
            if six_x_six == True:
                if not os.path.exists(path_data+'/corrected/6x6/'+star+'/'+band):
                    os.makedirs(path_data+'/corrected/6x6/'+star+'/'+band)
                hdu.writeto(path_data+'/corrected/6x6/'+star+'/'+band+'/'+img+'_c.fits',overwrite=True)
            else:
                if not os.path.exists(path_data+'/corrected/'+star+'/'+band):
                    os.makedirs(path_data+'/corrected/'+star+'/'+band)
                hdu.writeto(path_data+'/corrected/'+star+'/'+band+'/'+img+'_c.fits',overwrite=True)
            
            


# ===================================Platesolve============================================
def platesolve(star,band,solver:bool = True, six_x_six: bool = True):
    
    '''This function platesolves the images of a star in a given band.
    
    
    Parameters:
    -----------
    star : str
        Name of the star.
    band : str
        Band of the images.
    solver : bool
        If True, the function uses the online server. If False, the function uses the Giles server.'''
    plt.style.use('report.mplstyle')
    if six_x_six == True:
        
        path_data =  '/Volumes/external_2T/corrected/6x6'
    else:
        path_data =  '/Volumes/external_2T/corrected'
    
    counter = 0 
    
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
        
    star_list = os.listdir(path_data+'/'+star+'/'+band)
    for img in star_list:
        print(img)

            
        if img.endswith('.fits'):          
            with fits.open(path_data+'/'+star+'/'+band+'/'+img) as hdu:
                img_header = hdu[0].header
                img_data = hdu[0].data
            mean, median, std = sigma_clipped_stats(img_data)
            if six_x_six == True:
                
                daofind = DAOStarFinder(fwhm=12, threshold=5*std) 
            else:
                daofind = DAOStarFinder(fwhm=4.2, threshold=4*std)
            sources = daofind(img_data)

            for col in sources.colnames:  
                sources[col].info.format = '%.8g'
            sources.sort('flux')
            sources.reverse()
            # if counter == 0:
            #     plt.figure(figsize=(10,10))
            #     plt.imshow(img_data, cmap='gray', origin='lower')
            #     plt.colorbar()
            #     plt.title('Corrected image')
            #     plt.scatter(sources['xcentroid'], sources['ycentroid'], s=100, edgecolor='red', facecolor='none')
            #     plt.show()
            print(len(sources))
            #     counter += 1
            try:
                wcs_header = None 
                wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'],
                                                        img_data.shape[1], 
                                                        img_data.shape[0],
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
                hdu = fits.PrimaryHDU(img_data)
                hdu.header.update(img_header) 
                if six_x_six == True:
                    if not os.path.exists(path_data+'/wcs/6x6/'+star+'/'+band):
                        os.makedirs(path_data+'/wcs/6x6/'+star+'/'+band)
                    hdu.writeto(path_data+'/wcs/6x6/'+star+'/'+band+'/'+img,overwrite=True)
                else:
                    if not os.path.exists(path_data+'/wcs/'+star+'/'+band):
                        os.makedirs(path_data+'/wcs/'+star+'/'+band)
                    hdu.writeto(path_data+'/wcs/'+star+'/'+band+'/'+img,overwrite=True)
            
            except TimeoutError:       
                print('\n -> ##FAIL: Timeout while solving, try a longer timeout, optmise number of sources (200-800 seems about right)')
    
    
    
    

# ===================================Photometry============================================

def photometry(star: str, band: str, radius: int,cal_index: int = 4, six_x_six: bool = True):
    
    '''
    Completes photometry on pre calibrated images and platesolved images.
    
    Parameters:
    star (str): Name of the star
    band (str): Filter band
    radius (int): Radius of the aperture for photometry
    cal_index (int): Index of the calibration star (default value is 4 the target star)
    
    
    '''
    
    
    
    counter = 0
    
    if six_x_six == True:
        cal_star_path = '/Volumes/external_2T'
        path_data =  '/Volumes/external_2T/corrected/wcs/6x6'
    else:
        cal_star_path = '/Volumes/external_2T'
        path_data =  '/Volumes/external_2T/corrected/wcs'
        
    
    target_jd = []
    cal_star_mags = []
    air_mass = []
    
    if six_x_six == True:
        data_cal = np.transpose(np.loadtxt(cal_star_path + '/'+ 'cal_stars/' + star[0:10] + '_calibration_stars.txt', skiprows=1, delimiter=","))
    else:
        data_cal = np.transpose(np.loadtxt(cal_star_path + '/'+ 'cal_stars/' + star + '_calibration_stars.txt', skiprows=1, delimiter=","))
    mag_g_cal=data_cal[3]
    mag_g_err = data_cal[4]
    mag_i_cal=data_cal[7]
    mag_i_err = data_cal[8]
    mag_r_cal=data_cal[5]
    mag_r_err = data_cal[6]
    ra_cal=data_cal[1]
    dec_cal=data_cal[2]
    
    if re.search('G',band):
        mag_cal1 = mag_g_cal
        mag_err = mag_g_err
        
    if re.search('I',band):
        mag_cal1 = mag_i_cal
        mag_err = mag_i_err
        
    if re.search('R',band):
        mag_cal1 = mag_r_cal
        mag_err = mag_r_err
        
    star_list = os.listdir(path_data+'/'+star+'/'+band)     
    
    for img in star_list:
        print(img) 
        
        if img.endswith('.fits'):          
            with fits.open(path_data+'/'+star+'/'+band+'/'+img) as hdu:
                img_header = hdu[0].header
                img_data = hdu[0].data
                
            mean, median, std = sigma_clipped_stats(img_data)
            
            if six_x_six == True:
                
                daofind = DAOStarFinder(fwhm=12, threshold=5*std) 
                sources = daofind(img_data)
            else:
                daofind = DAOStarFinder(fwhm=4.2, threshold=4*std)
                sources = daofind(img_data)
            for col in sources.colnames:  
                sources[col].info.format = '%.8g'
            sources.sort('flux')
            sources.reverse()
            if counter == 0:
                print(img_header['OBJCTDEC'],img_header['OBJCTRA'])
                print((sources))
            wcs = WCS(img_header)
        
            # number of calibration stars to plot
            # print(wcs.wcs_pix2world(sources['xcentroid'], sources['ycentroid'], 1))


            # plot yellow circles around the sources found by DAO starfinder
            positions_dao = np.transpose((sources['xcentroid'], sources['ycentroid']))  


            r1=radius
            r2=r1+2
            r3=r2+4

            # plot blue circles around the calibration stars from VizieR
            source1_x, source1_y= wcs.wcs_world2pix(ra_cal,dec_cal,1)
            source1 = np.transpose((source1_x, source1_y))
            source1_aperture = CircularAperture(source1, r1)  
            source1_annulus = CircularAnnulus(source1, r2, r3)

            source1_phot = [source1_aperture, source1_annulus]
            # source1_aperture.plot(color='blue', lw=2, alpha=1)
            # source1_annulus.plot(color='deepskyblue', lw=2, alpha=1)

            source2_x, source2_y =(sources['xcentroid'] , sources['ycentroid'] )

            # Coor_x, Coor_y = wcs.wcs_pix2world(source2_x, source2_y, 1) 
            source2=np.transpose((source2_x, source2_y))
            source2_aperture = CircularAperture(source2, r1)  
            # source2_aperture.plot(color='red', lw=2, alpha=1)

            source2 = np.transpose((source2_x, source2_y))
            source2_aperture = CircularAperture(source2, r1)  
            source2_annulus = CircularAnnulus(source2, r2, r3)

            source2_phot = [source2_aperture, source2_annulus]
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++APERTURE PHOTOMETRY+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            phot_table_source1 = aperture_photometry(img_data, source1_phot)
            phot_table_source2 = aperture_photometry(img_data, source2_phot)

            for col in phot_table_source1.colnames:
                phot_table_source1[col].info.format = '%.8g'  # for consistent table output

            bkg_mean_cal1 = float(phot_table_source1[6]['aperture_sum_1'] / source1_annulus.area)
            bcal1 = bkg_mean_cal1 * source1_aperture.area

            cal1_flux=float(phot_table_source1[6]['aperture_sum_0'] - bcal1)


            for col in phot_table_source2.colnames:
                phot_table_source2[col].info.format = '%.8g'  # for consistent table output


            bkg_mean_targ = float(phot_table_source2[cal_index]['aperture_sum_1'] / source2_annulus.area)

            targcal = bkg_mean_targ * source2_aperture.area

            r = np.linspace(1, 30, 30)
            test_r2 = r+2
            test_r3 = r+4
    # =======================================air mass calculation===============================================================
            alt = img_header['OBJCTALT'][0:2]+'.'+img_header['OBJCTALT'][3:5]+img_header['OBJCTALT'][6:8]
            
            zenith_angle = 90 - float(alt)
            
            airmass = 1/np.cos(np.radians(zenith_angle))
            air_mass.append(airmass)

    # =======================================air mass calculation===============================================================
    # ---------------------------------------Growth Curve-------------------------------------------------------------------
            if counter == 0:
                flux = []
                plt.figure(figsize=(10,10))
                plt.ylabel('Flux e$^-$ pixel$^{-2}$')
                plt.xlabel('Radius (pixels)')
                
                for i in np.arange(1,31,1):
                    
                    source2 = np.transpose((source2_x, source2_y))
                    source2_aperture = CircularAperture(source2, i)  
                    source2_annulus = CircularAnnulus(source2, test_r2[i-1], test_r3[i-1])

                    source2_phot = [source2_aperture, source2_annulus]
                    phot_table_source2 = aperture_photometry(img_data, source2_phot)

                    bkg_mean_targ = float(phot_table_source2[cal_index]['aperture_sum_1'] / source2_annulus.area)

                    targcal = bkg_mean_targ * source2_aperture.area
                    targ_flux=float(phot_table_source2[cal_index]['aperture_sum_0'] - targcal)
                    flux.append(targ_flux)
                
                plt.plot(r,flux, color = 'blue', label = 'Curve of Growth')
                plt.vlines(9, min(flux), max(flux), color = 'red', linestyle = '--', label = 'Aperture Radius = 9')
                plt.legend()
                plt.savefig('graphs/growth_curve.pdf', bbox_inches='tight', pad_inches=0.1)
                plt.show()
                
                
    # ----------------------------------------------------------------------------------------------------------------------------
            targ_flux=float(phot_table_source2[cal_index]['aperture_sum_0'] - targcal)
            mag_targ=mag_cal1[4] + 2.5*np.log10(cal1_flux/targ_flux)

            t_fits=img_header['DATE-OBS']
            t = Time(t_fits, format='isot', scale='utc')
            t_jd=t.jd 
            cal_star_mags.append(mag_targ)
            target_jd.append(t_jd)
            # print(mag_targ)
            counter +=1
        
    dict = {'JD': target_jd, 'Mag': cal_star_mags, 'Air Mass': air_mass}
    df = pd.DataFrame(dict)
    # if six_x_six == True:
    #     df.to_csv('airmass/6x6/'+band[4]+'/cal_star_'+str(cal_index)+'_'+band[4]+'_airmass.csv')
    # else:
    #     df.to_csv('airmass/'+star+'/'+band[4]+'/cal_star_'+str(cal_index)+'_'+band[4]+'_airmass.csv')
    return target_jd, cal_star_mags, air_mass
    