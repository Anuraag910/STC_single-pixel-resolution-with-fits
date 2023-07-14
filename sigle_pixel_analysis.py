#!/usr/bin/env python
# coding: utf-8

# In[148]:


from astropy.io import fits
import numpy as np
import astropy
from astropy.table import Table
import matplotlib.pyplot as plt


# In[161]:


path = '/home/suman/tifr/'
file = path + '202307011_1506_Am241_on_det1_30000pkts.fits'

hdul = fits.open(file)[1]


# In[164]:


data = hdul.data
tab = Table(data)
tab


# In[165]:


data['pixid']


# In[41]:


#mask = 0 < np.any(data['pixid']) < 3


# In[42]:


#plt.plot(data["pixid"][mask],data["pha"][mask])


# In[166]:


detid = data[data['detid']== 1]
data_ = detid[detid['pixid'] == 100]

np.shape(data_)
data_


# In[113]:


tab = Table(data)
#tab


# In[167]:


tab_ = Table(data_)
tab_


# In[168]:


len(np.unique(data_['time']))


# In[169]:


import matplotlib.pyplot as plt


# In[170]:


data_['pha']


# In[172]:


plt.plot((data_['time']),(data_['pha']),".")
plt.xlabel("Aqusiton time")
plt.ylabel("PHA")
plt.xscale("log")
plt.yscale("log")
#plt.xlim(0, 100000000)


# In[187]:


N,bins,_ = plt.hist(data_['pha'], bins=500)
plt.xlabel("PHA")
plt.ylabel("No. of Counts")
plt.show()


# In[174]:


len(N)


# In[175]:


bin_centers = (bins[:-1] + bins[1:]) / 2


# In[176]:


def gauss(x,amp,mean,stdev):
    return amp*np.exp(-(x-mean)**2/(2*stdev**2))


# In[177]:


gauss(1,10,9.3,0.3)


# In[178]:


p0 = [np.max(N),np.mean(data_['pha']),np.std(data_['pha'])]


# In[179]:


p0


# In[180]:


from scipy.optimize import curve_fit


# In[181]:


params, pcov = curve_fit(gauss,bin_centers,N,p0)
errors = np.sqrt(np.diag(pcov))


# In[182]:


print("Fit results")
for p,e in zip(params, errors):
    print(f"{p:0.1f} +- {e:0.1f}")
print("Approximate resolution : {:0.1f}%".format(100* 2.35 * params[2] / params[1]))


# In[183]:


fit_curve = gauss(bin_centers,*params)


# In[186]:


plt.plot(bin_centers, fit_curve, 'r-', label='Gaussian Fit')
N,bins,_ = plt.hist(data_['pha'], bins= len(bins)-1)
plt.xlabel("PHA")
plt.ylabel("No. of Counts")
plt.show()


# In[ ]:




