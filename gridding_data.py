import numpy as np
import numpy as np
import scipy.io
import xarray as xr
import os
import pandas as pd
from datetime import datetime, date, time, timedelta
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm


from data_loading import *

CTD = get_CTD_all_years()

lat_min = CTD.latitude.min()
lat_max = CTD.latitude.max()
lon_min = CTD.longitude.min()
lon_max = CTD.longitude.max()
depth_min = CTD.depth.min()
depth_max = CTD.depth.max()

N_lat = 20
N_lon = 20

lat_bins = np.linspace(lat_min, lat_max, N_lat+1)
lon_bins = np.linspace(lon_min, lon_max, N_lon+1)
depth_bins = np.array([0, 40, 80, 120, 160, 200, 300, 400, 500, 600, 700, 800])

lat_centers = (lat_bins[1:] + lat_bins[:-1]) / 2
lon_centers = (lon_bins[1:] + lon_bins[:-1]) / 2
depth_centers = (depth_bins[1:] + depth_bins[:-1]) / 2

N_depth = len(depth_centers)

points_lat, points_lon = np.meshgrid(lat_centers, lon_centers)

U = np.full((N_lat, N_lon, N_depth), np.nan)
V = np.full((N_lat, N_lon, N_depth), np.nan)
T = np.full((N_lat, N_lon, N_depth), np.nan)
S = np.full((N_lat, N_lon, N_depth), np.nan)

CTD_count = np.zeros((N_lat, N_lon, N_depth))
LADCP_count = np.zeros((N_lat, N_lon, N_depth))

def gridded_data(year):
    ds_CTD = get_CTD(year)
    ds_LADCP = get_LADCP(year)

    for (i,j,k) in tqdm(np.ndindex((N_lat, N_lon, N_depth)), total=N_lat*N_lon*N_depth):
        lat_min = lat_bins[i]
        lat_max = lat_bins[i+1]
        lon_min = lon_bins[j]
        lon_max = lon_bins[j+1]
        depth_min = depth_bins[k]
        depth_max = depth_bins[k+1]
        #print(lat_min, lat_max, lon_min, lon_max, depth_min, depth_max)

        sub_CTD = ds_CTD.where(
            (ds_CTD.latitude >= lat_min) & (ds_CTD.latitude < lat_max) &
            (ds_CTD.longitude >= lon_min) & (ds_CTD.longitude < lon_max) &
            (ds_CTD.depth >= depth_min) & (ds_CTD.depth < depth_max), drop=True
        )
        sub_LADCP = ds_LADCP.where(
            (ds_LADCP.latitude >= lat_min) & (ds_LADCP.latitude < lat_max) &
            (ds_LADCP.longitude >= lon_min) & (ds_LADCP.longitude < lon_max) &
            (ds_LADCP.depth >= depth_min) & (ds_LADCP.depth < depth_max), drop=True
        )
        
        U[i,j,k] = sub_LADCP['U'].mean().values
        V[i,j,k] = sub_LADCP['V'].mean().values

        T[i,j,k] = sub_CTD['temperature'].mean().values
        S[i,j,k] = sub_CTD['salinity'].mean().values

        CTD_count[i,j,k] = len(sub_CTD['temperature'])
        LADCP_count[i,j,k] = len(sub_LADCP['U'])

    return U, V, T, S, CTD_count, LADCP_count
