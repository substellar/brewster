import numpy as np
from bbconv import spex
from bbconv import convfwhm
from bbconv import convr

#**************************************************************************

# This hides the fortran convolution code so it works nicely with rest of
# code


#**************************************************************************

def spex_non_uniform(obspec,modspec):


    fluxout = spex(np.asfortranarray(obspec),np.asfortranarray(modspec))[0:obspec[0,:].size]

    return fluxout


def conv_uniform_FWHM(obspec,modspec,fwhm):

    fluxout = confwhm(np.asfortranarray(obspec),np.asfortranarray(modspec),fwhm)[0:obspec[0,:].size]
    
    return fluxout


        

def conv_uniform_R(obspec,modspec,R):

    fluxout = convr(np.asfortranarray(obspec),np.asfortranarray(modspec),R)[0:obspec[0,:].size]
    

    return fluxout
