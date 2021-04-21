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
