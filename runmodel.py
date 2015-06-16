
# coding: utf-8

# In[1]:

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import scipy
import numpy as np
#from IPython.display import display
#get_ipython().magic(u'matplotlib inline')
#%config InlineBackend.close_figures = False
from scipy.io.idl import readsav
import forwardmodel


# In[2]:
intemp = np.loadtxt("temps.dat",dtype='f')
temp = np.asfortranarray(intemp)
print temp

w1 = 0.8
w2 = 2.5
logg = 4.5
R2D2 = 1
gasnum = np.asfortranarray(np.array([1],dtype='i'))



# In[3]:

VMR1 = np.full((61,),8e-4)
#vmr2 = np.loadtxt("h2o.dat",unpack=True)
#vmr3 = np.loadtxt("K.dat",unpack=True)
#print VMR1.shape
VMR = np.reshape((VMR1),(1,61))
print VMR.shape


# In[4]:

pcover = 1.0
do_clouds = 0
cloudnum = np.array([1],dtype='i')
#cloudname = np.reshape((cname),(1,1))
cloudrad = np.full((1,61,1),1e-4)
cloudsig = np.full((1,61,1),1e-5)
cloudprof = np.full((1,61,1),0.0)


# In[ ]:
test_spec = forwardmodel.marv(w1,w2,temp,logg,R2D2,gasnum,VMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof)


# In[ ]:



