import numpy as np
from bbconv import prism
from bbconv import convfwhm
from bbconv import convr

#**************************************************************************

# This hides the fortran convolution code so it works nicely with rest of
# code


#**************************************************************************

def prism_non_uniform(obspec,modspec,resel):


    fluxout = prism(np.asfortranarray(obspec),np.asfortranarray(modspec),resel)[0:obspec[0,:].size]

    return fluxout


def conv_uniform_FWHM(obspec,modspec,fwhm):

    fluxout = convfwhm(np.asfortranarray(obspec),np.asfortranarray(modspec),fwhm)[0:obspec[0,:].size]
    
    return fluxout


        

def conv_uniform_R(obspec,modspec,R):

    fluxout = convr(np.asfortranarray(obspec),np.asfortranarray(modspec),R)[0:obspec[0,:].size]
    

    return fluxout
    
    
    
    
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################




#### CONVOLVING THE MODEL SPECTRA WITH THE NON-UNIFORM RESOLVING POWER (NIRSPEC + MIRI)

#wl =  obs_data[:, 0]  #np.arange(1,16,1)   lambda in microns



def JWST_R(wl):
    
    '''Function that provides the Resolving Power (R) as a function of wavelength (wl) in micron,
       depending on jwst instrument'''
    
    if slit_limited == 1.0 and source_limited == 0.0:
        
        R = wl/(d_wl * resel_slit_lim)
        
        return R, wl
    
    
    elif source_limited == 1.0 and slit_limited == 0.0:
        
        theta_lambda = (648000/np.pi) * 1.028 * (wl/D) #FWHM in arcsec
        
        #resel = ((648000/np.pi) * 1.028 * (wl/D)) / theta_pix
        
        R = wl / (d_wl * (theta_lambda / theta_pix))
        
        return R, wl
    
    
    elif source_limited == 1.0 and slit_limited == 1.0:
        
        ### nirspec
        
        nirspec_wl = wl[(wl >= nirspec_wl_range[0]) & (wl <= nirspec_wl_range[1])]
        
        nirspec_d_wl = d_wl[(wl >= nirspec_wl_range[0]) & (wl <= nirspec_wl_range[1])] 
        
        R_nirspec = nirspec_wl/(nirspec_d_wl * resel_slit_lim)
        
        ### miri
        
        miri_wl = wl[(wl >= miri_wl_range[0]) & (wl <= miri_wl_range[1])]
        
        miri_d_wl = d_wl[(wl >= miri_wl_range[0]) & (wl <= miri_wl_range[1])]
        
        #resel = ((648000/np.pi) * 1.028 * (miri_wl/D)) / theta_pix
        
        theta_lambda = (648000/np.pi) * 1.028 * (miri_wl/D)
        
        R_miri = miri_wl / (miri_d_wl * (theta_lambda / theta_pix))
        
        #plt.plot(nirspec_wl,R_nirspec)
        #plt.plot(miri_wl,R_miri)
        
        return np.concatenate((R_nirspec, R_miri)), np.concatenate((nirspec_wl, miri_wl))
    
    else:
        
        print('Error: please check source_limited and slit_limited.')
        



def conv_non_uniform_R(model_flux, model_wl, R, obs_wl):
    """
    Convolve a model spectrum with a wavelength-dependent resolving power 
    onto the observed wavelength grid ???

    Parameters:
    - model_flux: 1D array of model flux values.
    - model_wl: 1D array of model wl values.
    - obs_wl: 1D array of observed wl values.
    - R: 1D array of resolving power values (for the obs_wl grid.)

    Returns:
    - convolved_flux: 1D array of convolved flux values on the obs_wl grid.
    """
    # create the array for the convolved flux
    convolved_flux = np.zeros_like(obs_wl)

    for i, wl_center in enumerate(obs_wl): 
        
        # compute FWHM and sigma for each wl
        # print('wl_center', wl_center)
        # print('R[i]', R[i])
        
        fwhm = wl_center / R[i]
        # print('fwhm', fwhm)
        sigma = fwhm / 2.355


        # compute the Gaussian kernel for the current wl
       
        gaussian_kernel = np.exp(-((model_wl-wl_center) ** 2) / (2 * sigma **2))
        #print('gaussian_kernel before normalisation', gaussian_kernel)

        # normalisation
        gaussian_kernel /= np.sum(gaussian_kernel)
        # print('gaussian_kernel after normalisation', gaussian_kernel)



        # apply the kernel to the flux
        convolved_flux[i] = np.sum(model_flux * gaussian_kernel)
    
    return convolved_flux
