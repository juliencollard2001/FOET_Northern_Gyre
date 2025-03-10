import numpy as np
import imageio
import xarray as xr
import os
import pandas as pd
from datetime import datetime, date, time
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
from tqdm import tqdm
import os
import warnings

warnings.filterwarnings('ignore')

from data_loading import get_altimetry_data, get_2021_CTD_data
import warnings

proj = ccrs.Mercator(central_longitude=4.5, min_latitude=38.0, max_latitude=45.0, latitude_true_scale=42.0)
set = [-1, 10, 38, 45]

def plot_one_day(ds, date, path):
    sub_ds = ds.sel(time=date)

    Vgos = sub_ds['vgos'] *  np.cos( sub_ds.latitude / 180 * np.pi)
    Norm = np.sqrt(sub_ds['ugos']**2+ sub_ds['vgos']**2) / np.sqrt(sub_ds['ugos']**2 + Vgos**2)

    plt.figure(figsize=(10,10))
    ax = plt.axes(projection=proj)
    ax.set_extent(set)
    pc = ax.pcolormesh(
        sub_ds.longitude, 
        sub_ds.latitude, 
        sub_ds['gos_current_norm'], 
        transform=ccrs.PlateCarree(), 
        cmap='YlGnBu',
        alpha=0.8,
        vmax = 1.2,
        vmin= 0
    )
    #pc = ax.pcolormesh(ds_mean.longitude, ds_mean.latitude,ds_mean['ugos'], transform=ccrs.PlateCarree(), cmap='coolwarm', alpha=0.5)
    q = ax.quiver(
        sub_ds.longitude, 
        sub_ds.latitude, 
        sub_ds['ugos']*Norm, 
        Vgos*Norm, 
        transform=ccrs.PlateCarree(), 
        regrid_shape=50,
        scale=10,
        width=0.001
    )
    ax.quiverkey(q, X=0.87, Y=0.95, U=0.25, label='0.25 m/s', labelpos='E', transform=ax.transAxes)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.gridlines(draw_labels=True)
    #plt.subplots_adjust(top=0.9)
    plt.colorbar(pc, orientation='horizontal', label='speed in m/s', fraction=0.046, pad=0.04)
    plt.title(str(date), fontsize=14)
    plt.subplots_adjust(top=0.85)
    plt.savefig(path)
    plt.close()

years = [2016, 2017, 2018]


if not os.path.exists('film'):
    os.makedirs('film')
else:
    os.system('rm film/*')

for year in years:
    print(f"Processing year {year}")
    # Load altimetry data
    ds = get_altimetry_data([year],[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ds['gos_current_norm'] = np.sqrt(ds['ugos']**2 + ds['vgos']**2)

    for date in tqdm(ds.time.values):
        plot_one_day(ds, date, f"film/{date}.png")
        
film_name = f'film_{years[0]}_{years[-1]}.gif'
filenames = os.listdir('film')
filenames.sort()
images = []
for filename in tqdm(filenames):
    images.append(imageio.imread('film/'+filename))
imageio.mimsave(film_name, images, fps=15)

os.system('rm film/*')