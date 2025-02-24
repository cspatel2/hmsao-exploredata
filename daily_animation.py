# %%
from __future__ import annotations
import pytz
import astropy.io.fits as fits
import os
from glob import glob
from matplotlib import pyplot as plt
from dateutil.parser import parse
from typing import SupportsFloat as Numeric
import subprocess
from typing import Iterable
from tqdm import tqdm
import xarray as xr
import numpy as np
from datetime import datetime, timezone, timedelta
# import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
# from uncertainties import ufloat, unumpy as unp

# %%
# fdir = '../hmsao-processed'
# files = glob(os.path.join(fdir, '*.nc'))
# len(files)

# %%


# def get_date_from_fn(fn: str):
#     if isinstance(fn, list):
#         return [get_date_from_fn(f) for f in fn]
#     return os.path.basename(fn).split('_')[1]

# def get_data_from_fn(fn:str):
#     bn = os.path.basename(fn)
#     words = bn.split('_')
#     w = words[-1].strip('.nc').split('[')
#     wl = w[0]
#     idx = w[-1].strip(']')
#     return wl,idx


# %%
# sweden = pytz.timezone('Europe/Stockholm')
# fdir = '../hmsao-processed'
fdir = '/mnt/d/Projects/hms-ao/hmsao-processed'
window = '5577'
date = '20250201'
files = glob(os.path.join(fdir, f'*{date}*{window}*.nc'))
print(len(files))
# %%

# %%
outdir = f'plots/{date}/{window}/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

ds = xr.open_mfdataset(files)
txt_fname = os.path.join(outdir,'filelist.txt')
with open(txt_fname, 'w') as ofile:
    datads = ds.diff(dim='tstamp')
    exposures = ds.exposure.values
    del ds
    for tidx, tstamp in tqdm(enumerate(datads.tstamp.values), total=len(datads.tstamp.values), desc='Generating plots'):
        data = datads.intensity[tidx]
        exptime = exposures[tidx+1]
        fig, ax = plt.subplots(figsize=(1900/300, 1080/300), dpi=300)
        vmin = np.percentile(data.values, 1)
        vmax = np.percentile(data.values, 99.5)
        time = datetime.fromtimestamp(tstamp).astimezone(pytz.utc)
        # cmap = 'gist_ncar_r'
        cmap = 'bwr'
        plot = data.plot(vmin=vmin, vmax=vmax, cmap=cmap)
        plot.colorbar.set_label('Diff in Intensity [ADU/s]')
        ax.axvline(int(window)/10, color='k', ls='--', lw=0.5)
        ax.set_title(
            f'{time:%Y-%m-%dT%H:%M:%S}\nExposure: {exptime:.2f} s', color='k')
        outfname = os.path.join(outdir, f'{tidx}.png')
        fig.savefig(outfname, bbox_inches='tight')
        plt.close(fig)
        # plt.show()
        ofile.write(f"file '{outfname}'\n")
        # ofile.write(f'duration {exptime:.5f}\n')
        # break

# %%
ani_outdir = f'animations/{date}'
if not os.path.exists(ani_outdir):
    os.makedirs(ani_outdir)

framerate = 10 #output framerate 
# input_fname = txt_fname #input file
input_fname = '/plots/20250201/5577/filelist.txt'
pixel_format = 'yuv420p' 
ani_outfname = os.path.join(ani_outdir, f'hmsao-{window}.mp4') #output filename

command = f'ffmpeg -r {framerate} -f concat -safe 0 -i {input_fname}  -pix_fmt {pixel_format} {ani_outfname} -y'

subprocess.run(command, shell=True)

# %%
# subprocess.run(
#     'ffmpeg -f concat -i filelist.txt -vsync cfr -pix_fmt yuv420p hitmis_originlaunch.mp4 -y', shell=True)
