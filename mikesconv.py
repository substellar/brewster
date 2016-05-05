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

def instrument_non_uniform(wlgrid,wno, Fp):

    #delta=wlgrid-np.roll(wlgrid,1)

    #delta[0]=delta[1]



    ###more accurate approximation

    delta1=wlgrid-np.roll(wlgrid,1)

    delta2=np.abs(wlgrid-np.roll(wlgrid,-1))

    delta=0.5*(delta1+delta2)

    delta[0]=delta[1]

    delta[-1]=delta[-2]



    szmod=wlgrid.shape[0]

    Fratio_int=np.zeros(szmod)

    wl=1E4/wno[::-1]

    Fp=Fp[::-1]

    for i in range(wlgrid.shape[0]):

        sigma=delta[i]/2.355

        gauss=np.exp(-(wl-wlgrid[i])**2/(2*sigma**2))

        gauss=gauss/np.sum(gauss)

        Fratio_int[i]=np.sum(gauss*Fp)



    return Fratio_int
