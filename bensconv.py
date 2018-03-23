import numpy as np

#**************************************************************************

# FILE:instrument_non_uniform.py

#

# DESCRIPTION: This function takes in a temperature (T) and a

# wavelength grid (wl) and returns a blackbody flux grid.

#

# USAGE: 1. Import function ---- '>>> from blackbody import blackbody'

#        2. Call blackbody.py -- '>>> B = blackbody(T,wl)'

#**************************************************************************

def spex_non_uniform(wlgrid,wno, Fp):

    #delta=wlgrid-np.roll(wlgrid,1)

    #delta[0]=delta[1]

    wlgrid = np.reshape(wlgrid, (wlgrid.size))
    wno = np.reshape(wno,(wno.size))
    Fp = np.reshape(Fp,(Fp.size))

    ###more accurate approximation

    delta1=wlgrid-np.roll(wlgrid,1)

    delta2=np.abs(wlgrid-np.roll(wlgrid,-1))

    delta=0.5*(delta1+delta2)

    
    delta[0]=delta[1]

    delta[-1]=delta[-2]

    szobs=wlgrid.shape[0]

    Fratio_int=np.zeros(szobs)

    wl=1E4/wno[::-1]

    Fp=Fp[::-1]

    for i in range(szobs):

        sigma=delta[i] * 3.3 / 2.355

        gauss=np.exp(-(wl-wlgrid[i])**2/(2*sigma**2))

        gauss=gauss/np.sum(gauss)

        Fratio_int[i]=np.sum(gauss*Fp)



    return Fratio_int


def conv_uniform_FWHM(modspec,obspec,fwhm):

    wlmod = np.reshape(modspec[0,:],modspec[0,:].size)
    Fp = np.reshape(modspec[1,:],modspec[1,:].size)
    wlobs = np.reshape(obspec[0,:],obspec[0,:].size)
    
    szobs=wlobs.shape[0]

    Fratio_int=np.zeros(szobs)


    for i in range(szobs):

        sigma = fwhm / 2.355

        gauss = np.exp(-(wlmod-wlobs[i])**2/(2*sigma**2))

        gauss = gauss/np.sum(gauss)

        Fratio_int[i]=np.sum(gauss*Fp)

    return Fratio_int


        

def conv_uniform_R(modspec,obspec,R):

    wlmod = np.reshape(modspec[0,:],modspec[0,:].size)
    Fp = np.reshape(modspec[1,:],modspec[1,:].size)
    wlobs = np.reshape(obspec[0,:],obspec[0,:].size)
    
    szobs=wlobs.shape[0]
    
    Fratio_int=np.zeros(szobs)


    for i in range(szobs):
        # sigma is FWHM / 2.355
        sigma = (wlmod / R) / 2.355

        gauss = np.exp(-(wlmod-wlobs[i])**2/(2*sigma**2))

        gauss = gauss/np.sum(gauss)

        Fratio_int[i]=np.sum(gauss*Fp)



    return Fratio_int
