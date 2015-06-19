import numpy as np
import scipy as sp
import forwardmodel

def run(intemp,invmr,pcover,cloudparams,r2d2,grav,dlam,do_clouds,gasnum,cloudnum):
    # get the ngas and ncloud bits
    ngas = invmr.shape[0]
    inlayer = np.arange(0,15.25,1)
    layer = np.arange(0,15.25,0.25)
    # Hard code nlayers
    nlayers = 61
    # spline fit with no smoothing)
    tfit = np.interpolate.splrep(inlayer,temp,s=0)
    temp = np.interpolate.splev(layer,tfit, der=0)

    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile    
    VMR = np.empty((ngas,nlayers),dtype='d')
    if invmr.size > invmr.shape[0]:
        for i in range(0,ngas):
            vfit = np.interpolate.splrep(inlayer,invmr[i,:],s=0)
            VMR[i,:] = np.interpolate.splev(layer,vfit,der=0)
    else:
        for i in range(0,ngas):
            VMR[i,:] = invmr[i]

    # now need to translate cloudparams in to cloud profile
    if do_clouds == 1:
        npatch = cloudparams.shape[0]
        ncloud = cloudparams.shape[1]
        # 3rd dimension of cloudparams for simple slab model are:
        # 0) number density
        # 1) top layer id (or pressure)
        # 2) base ID (these are both in 61 layers)
        # 3) rg
        # 4) rsig
        cloudrad = np.empty((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.empty_like(cloudrad)
        cloudprof = np.zeros_like(cloudrad)
        ndens= np.reshape(cloudparams['f0'],(npatch,ncloud))
        c1 = np.reshape(cloudparams['f1'],(npatch,ncloud))
        c2 = np.reshape(cloudparams['f2'],(npatch,ncloud))
        rad = np.reshape(cloudparams['f3'],(npatch,ncloud))
        sig = np.reshape(cloudparams['f4'],(npatch,ncloud))
        for i in range(0, npatch):
            for j in range(0, ncloud):
                b1 = c1[i,j] - 1
                b2 = c2[i,j] -1 
                cloudprof[i,b1:b2+1,j] = ndens[i,j]
                cloudrad[i,:,j] = rad[i,j]
                cloudsig[i,:,j] = sig[i,j]        
        
    # now we can call the forward model

    spec = forwardmodel.marv(w1,w2,temp,logg,R2D2,gasnum,VMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof)


    # outspec now 2d wave,flux
    # now shift wavelen by delta_lambda and regrid the flux

    
    # convolve with instrumental profile


    # rebin to observed dispersion
