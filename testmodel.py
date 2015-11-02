
# coding: utf-8

# In[1]:

import scipy as sp
import numpy as np
import pickle
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.io.idl import readsav
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from pysynphot import observation
from pysynphot import spectrum
import forwardmodel
import ciamod


# In[2]:

# set up coarse pressure grid and fine pressure grid
logcoarsePress = np.arange(-4.0, 2.5, 0.5)
coarsePress = pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.5, 0.1)
finePress = pow(10,logfinePress)
print finePress.size
print coarsePress.size
print finePress


# In[3]:

array = pickle.load(open("test_H2H2_H2He_CIA_H2O.pic", "rb")) 
leveltemp = array[0]
levelpress = array[1]
mikespec = np.array([array[2],array[3]],dtype='f')
mikespec[0] = 10000.0 / mikespec[0]
print levelpress.size
print levelpress


# In[4]:

mikepress = np.empty(levelpress.size - 1,dtype='float64')
miketemp = np.empty(leveltemp.size -1, dtype='float64')
for i in range(0,mikepress.size):
    mikepress[i] = np.sqrt(levelpress[i] * levelpress[i+1])
mtfit = interp1d(np.log10(levelpress),leveltemp)
miketemp = mtfit(np.log10(mikepress)) 
tfit = interp1d(np.log10(mikepress),miketemp,bounds_error=False,fill_value=miketemp[miketemp.size-1])
temp = tfit(np.log10(finePress))
temp[0] = temp[1] - 3.0
print temp
#temp[:] = 1000.
#print temp


# In[5]:

# Get the linelists
ngas = 5
gasnum = np.asfortranarray(np.array([1,2,20,4,5],dtype='i'))
lists = ["../Linelists/xsecarrH2O_1wno_500_10000.save","../Linelists/xsecarrCH4_1wno_500_10000.save","../Linelists/xsecarrK_new_1wno_500_10000_02.save","../Linelists/xsecarrCO_1wno_500_10000_02.save","../Linelists/xsecarrCO2_1wno_500_10000_02.save" ]

#lists = ["/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrCH4_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrK_new_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrCO_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrCO2_1wno_500_10000_02.save" ]

# In[6]:
#x=readsav('/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save')
x=readsav('../Linelists/xsecarrH2O_1wno_500_10000.save')
inlinelist=x.xsecarr  #3D array with Nwavenubmers x Ntemps x Npressure
inlinetemps=np.asfortranarray(x.t,dtype='float64')
inpress=x.p
inwavenum=x.wno
ntemps = inlinetemps.size
npress= finePress.size
nwave = inwavenum.size
#logpress = np.arange(-5.,2.5,0.125)
#press = 10.**logpress
#print press
# Here we are interpolating the linelist onto our fine pressure scale. 
linelist = (np.ones([ngas,npress,ntemps,nwave],order='F')).astype('float64', order='F')
for gas in range (0,ngas):
    inlinelist=readsav(lists[gas]).xsecarr
    for i in range (0,ntemps):
        for j in range (0,nwave):
            pfit = interp1d(np.log10(inpress),np.log10(inlinelist[:,i,j]))
            linelist[gas,:,i,j] = np.asfortranarray(pfit(np.log10(finePress)))
#print linelist.shape
#print np.result_type(linelist)
press = finePress*1000.


# In[7]:

#intemp = np.loadtxt("16temps.dat",dtype='float32')
#intemp = np.full(16,1000.,dtype='f')
#inlayer = np.arange(0,15.25,1)  
#layer = np.arange(0,15.00,0.25)
#print layer.size
#tfit = interpolate.splrep(inlayer,intemp,s=0)
#temp = np.asfortranarray(interpolate.splev(layer,tfit, der=0),dtype='float32')
w1 = 1.05
w2 = 5.0
logg = 4.5
R2D2 = 1.0
use_disort = 0
#print inpress
#print mikepress


# In[8]:

VMR1 = np.full((npress,),(-3.5)) # water
vmr2 = np.full((npress,),(-3.4)) # ch4
vmr3 = np.full((npress,),(-8.0))  # K
vmr4 = np.full((npress,),(-7.5)) # CO
vmr5 = np.full((npress,),(-8.2)) # CO2
#print VMR1.shape
logVMR = np.asfortranarray(np.reshape((VMR1,vmr2,vmr3,vmr4,vmr5),(ngas,npress)),dtype='float64')
print logVMR.shape


# In[9]:

pcover = 1.0
do_clouds = 0
cloudnum = np.array([1],dtype='i')
#cloudname = np.reshape((cname),(1,1))
cloudrad = np.full((1,npress,1),1e-4)
cloudsig = np.full((1,npress,1),1e-5)
cloudprof = np.full((1,npress,1),0.0)


# In[10]:

cia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(cia, dtype='float32')
ciatemps = np.asfortranarray(ciatemps, dtype='float32')
print cia.shape
print ciatemps.shape
print cia.dtype


# In[11]:

outspec = forwardmodel.marv(w1,w2,temp,logg,R2D2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort)


