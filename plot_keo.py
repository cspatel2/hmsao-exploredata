#%%
from matplotlib import colors
import numpy as np 
import matplotlib.pyplot as plt
import xarray as xr
import os
from glob import glob
from datetime import datetime,timedelta
# %%
fdir = '../hmsao-processed'
files = glob(os.path.join(fdir, '*.nc'))
len(files)

# %%
def get_date_from_fn(fn:str):
    if isinstance(fn,list):
        return [get_date_from_fn(f)  for f in fn]
    return os.path.basename(fn).split('_')[1]

# def get_data_from_fn(fn:str):
#     bn = os.path.basename(fn)
#     words = bn.split('_')
#     w = words[-1].strip('.nc').split('[')
#     wl = w[0]
#     idx = w[-1].strip(']')
#     return wl,idx
#%%

window = '5577'
date = '*'
wlfiles = glob(os.path.join(fdir,f'{date}*{window}*.nc'))
print(len(wlfiles))


# %%
dates = get_date_from_fn(wlfiles)
dates = np.unique(dates)

#%%
date = dates[3]
files = [f for f in files if date in f]
print(len(files))
# %%
ds = xr.open_mfdataset(files)
#%%
time = [datetime.fromtimestamp(t) for t in ds['tstamp'].values]
#%%
ds = ds.assign_coords(time = ('tstamp',time))
#%%


data = ds.intensity.sel(gamma = slice(40,52)).sum(dim='wavelength').diff(dim = 'timstamp')
#%%
vmin = np.nanpercentile(data,0.1)
vmax = np.nanpercentile(data,99.9)
#%%
data.plot(vmin = 1e-3, y = 'gamma', x = 'time', norm=colors.LogNorm())



# %%

# %%
ds
# %%
t_start = datetime(2025,1,17,9,0).timestamp()
t_end = t_start + timedelta(hours=9).timestamp()



# %%
