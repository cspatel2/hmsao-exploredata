# %%
from __future__ import annotations
import pytz
import astropy.io.fits as fits
import os
from glob import glob
import argparse
from typing import SupportsFloat as Numeric
import subprocess
from typing import Iterable
from tqdm import tqdm
import xarray as xr
import numpy as np
from datetime import datetime, timezone, timedelta
import matplotlib.pyplot as plt
from itertools import product
import sys 
from time import perf_counter_ns
from collections import Counter

#needs to be commented out for running an interactive window
import matplotlib as mpl
mpl.use('Agg') #Fixes RuntimeError: main thread is not in main loop with Matplotlib and Flask.


usetex = False
if not usetex:mpl.rcParams.update({'mathtext.fontset': 'cm'})
mpl.rc('font', **{'family': 'serif',
       'serif': ['Times' if usetex else 'Times New Roman']})
mpl.rc('text', usetex=usetex)
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)
mpl.rc('axes', labelsize=10)


# %% functions


def get_window_from_fname(fname: str) -> str:
    if isinstance(fname, list):
        return [get_window_from_fname(f) for f in fname]
    return fname.split('_')[-1].strip('.nc').split('[')[0]


def get_date_from_fname(fname: str) -> str:
    if isinstance(fname, list):
        return [get_date_from_fname(f) for f in fname]
    return fname.split('_')[-2]
# %%
parser = argparse.ArgumentParser(
    description='Creates animations of intensity and intensity difference for HMSAO L1A data.'
)
# %%
parser.add_argument(
    'rootdir',
    type=str,
    nargs='?',
    help='Root Directory containing HMSAO L1A data'
)

parser.add_argument(
    '--dest',
    metavar='dest',
    type=str,
    required=False,
    default=os.getcwd(),
    nargs='?',
    help='Root directory where animation and plots will be stored.'
)


def list_of_strings(arg: str) -> List[str]:
    a = arg.split(',')
    return [s.strip() for s in a]


parser.add_argument(
    '--windows',
    type=list_of_strings,
    required=True,
    nargs='?',
    help='Window(s) to process (i.e. --window 1235, 3456, ...).'
)
parser.add_argument(
    '--yyyymmdd',
    type=list_of_strings,
    required=False,
    default = None,
    nargs='?',
    help='dates to process (i.e. --yyyymmdd 20250101, 20250213,...).'
)

parser.add_argument(
    '--overwrite',
    type=bool,
    required=False,
    default=False,
    nargs='?',
    help='If you want to rewrite an existing file, then True. Defaults to False.'
)
# %%
args = parser.parse_args()

# check root path
rootdir = args.rootdir
if not os.path.exists(rootdir):
    raise FileNotFoundError(f'Root directory {rootdir} not found.')
elif os.path.exists(rootdir) and not os.path.isdir(rootdir):
    if not os.listdir(rootdir):
        raise ValueError(f'{rootdir} is empty.')
    else:
        raise NotADirectoryError(f'{rootdir} is not a directory.')

# get all files from rootdir
files = glob(os.path.join(rootdir, '*.nc'))
if not files:
    raise FileNotFoundError(f'No .nc files found in {rootdir}')

# check valid time and window
all_valid_windows = get_window_from_fname(files)
all_valid_windows = np.sort(list(set(all_valid_windows)))
valid_windows = [str(w) for w in args.windows if w in all_valid_windows]
print(f'Valid windows that will be processed ({len(valid_windows)}/{len(args.windows)}): {valid_windows}')

all_valid_dates = get_date_from_fname(files)
all_valid_dates = np.sort(list(set(all_valid_dates)))
if args.yyyymmdd is None:
    valid_dates = all_valid_dates
else:
    valid_dates = [str(d) for d in args.yyyymmdd if d in all_valid_dates]
print(f'Valid dates that will be processed ({len(valid_dates)}/{len(args.yyyymmdd)}): {valid_dates}')

# get all pairs of valid dates and windows
date_window_iter = list(product(valid_dates, valid_windows))

for pair in date_window_iter:
    date, window = pair

    # check/make outdir for plots
    outdir_plots = os.path.join(args.dest, 'plots', f'{date}', f'{window}')
    os.makedirs(outdir_plots, exist_ok=True, mode=777)


    # check/make outdir for animations
    outdir_ani = os.path.join(args.dest, 'animations', f'{date}')
    os.makedirs(outdir_ani, exist_ok=True, mode=777)
    ani_outfname = os.path.join(outdir_ani, f'hmsao-{window}.mp4')
    #if animation exists and not overwrite, skip (dont make new plots)
    if os.path.exists(ani_outfname) and not args.overwrite:
        print(f'{ani_outfname} exists. Skipping.')
        continue

    # filename for  .txt file inputfile to ffmpeg
    txt_fname = os.path.join(outdir_plots, f'filelist_{window}.txt')

    # get valid .nc files for date and window
    pfiles = glob(os.path.join(rootdir, f'*{date}*{window}*.nc'))
    pfiles.sort()

    with xr.open_mfdataset(pfiles) as ds:
        print('Opened dataset')
        tslen = len(ds.tstamp.values)

        #open .txt file
        with open(txt_fname, 'w') as ofile:
            for tidx in tqdm(range(1, tslen),  total=tslen - 1, desc=f'{datetime.strptime(date,'%Y%m%y'):%Y-%m-%d}/{window}'):
                time = datetime.fromtimestamp(ds.tstamp.values[tidx]).astimezone(pytz.utc)
                #create figure
                fig, ax = plt.subplots(1, 2, figsize=(
                    2002/300, 1100/300), dpi=300, gridspec_kw={'wspace': 0.5})
                fig.suptitle(f'{time:%Y-%m-%dT%H:%M:%S}', color='k')
                lw = 0.3
                color = 'k'
                ls = '--'

                #left intensity plot
                data = ds.intensity[tidx]
                vmin = np.nanpercentile(data.values, 1)
                vmax = np.nanpercentile(data.values, 99.5)
                plot = data.plot(vmin=vmin, vmax=vmax,cmap='viridis', ax=ax[0])
                ax[0].axvline(int(window)/10, color=color, ls=ls, lw=lw)
                plot.colorbar.set_label('Intensity [ADU/s]')
                ax[0].set_title('Counts')

                #right intensity difference plot
                data = ds.intensity[tidx] - ds.intensity[tidx - 1]
                vmin = np.nanpercentile(data.values, 1)
                vmax = np.nanpercentile(data.values, 99.5)
                plot = data.plot(vmin=vmin, vmax=vmax, cmap='bwr', ax=ax[1])
                ax[1].axvline(int(window)/10, color=color, ls=ls, lw=lw)
                plot.colorbar.set_label('Intensity Difference [ADU/s]')
                ax[1].set_title('Difference')
                
                #save figure
                outfname = os.path.join(outdir_plots, f'{tidx}.png')
                fig.savefig(outfname, bbox_inches='tight')
                plt.close(fig)
                # plt.show()

                #write to .txt file
                ofile.write(f"file '{outfname}'\n")
                ofile.write(f'duration {1/5}\n') # 15 fps

    # create animation
    print('Creating animation...')
    sys.stdout.flush()
    tstart = perf_counter_ns()
                
    pixel_format = 'yuv420p'
    command = f'ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i {txt_fname} -pix_fmt {pixel_format} {ani_outfname} -y'
    subprocess.run(command, shell=True)
    tend = perf_counter_ns()
    print(f'Done. [{(tend-tstart)*1e-9:.3f} s]')
    sys.stdout.flush()
