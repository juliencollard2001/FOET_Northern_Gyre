import numpy as np
import scipy.io
import xarray as xr
import os
import pandas as pd
from datetime import datetime, date, time


def get_altimetry_data(years : list[int] , months : list[int]) -> xr.Dataset:
    """
    Load altimetry data from the MOOSE-Altimetry dataset for the specified years and months.
    The data is loaded from the specified folders.
    """

    base_folder = 'data/MOOSE-Altimetry'
    ds = None
    folders = [base_folder + '/' + str(year) + '/' + str(month).zfill(2) for year in years for month in months]
    for folder in folders:
        print(f'Loading data from {folder}...')
        files = os.listdir(folder)
        for file in files:
            if file.endswith('.nc'):
                try:
                    if ds is None:
                        ds = xr.open_dataset(os.path.join(folder, file))
                    else:
                        ds = xr.concat([ds, xr.open_dataset(os.path.join(folder, file))], dim='time')
                
                except Exception as e:
                    print(f'Error: {e}')
                    continue
    print('Done!')
    return ds


def get_2021_CTD_data() -> xr.Dataset:
    """
    Load the CTD data for the year 2021.
    """
    base_folder = 'data/moose-cruises'
    file = base_folder + '/CTD2021_MOOSE_GE.mat'
    data = scipy.io.loadmat(file)['all_cnv']

    columns_names = data[0,:].dtype.names
    colmmns_types = [type(data[0,i][0]) for i in range(len(columns_names))]

    new_columns = ['idx', 'file_name', 'date', 'lat', 'lon', 'depth', 'temperature', 'salinity', 'pressure']
    cols = [[] for i in range(len(new_columns))]

    idxs = np.arange(data.shape[1])
    depths = data[0,0][8][:,0]

    lats = np.empty(data.shape[1])
    lons = np.empty(data.shape[1])
    dates = np.empty(data.shape[1], dtype='datetime64[D]')
    file_names = np.empty(data.shape[1], dtype='object')
    temperatures = np.empty((data.shape[1], len(depths)))
    salinities = np.empty((data.shape[1], len(depths)))
    pressures = np.empty((data.shape[1], len(depths)))

    for i in range(data.shape[1]):
        file_names[i] = data[0,i][0][0]
        dates[i] = np.datetime64(datetime.strptime(data[0,i][3][0], '%d-%b-%Y %H:%M:%S'))
        lats[i] = float(data[0,i][6][0])
        lons[i] = float(data[0,i][5][0])
        temperatures[i] = data[0,i][9][:,0]
        salinities[i] = data[0,i][10][:,0]
        pressures[i] = data[0,i][7][:,0]

    ds = xr.Dataset(
        {
            'temperature': (['idx', 'depth'], temperatures),
            'salinity': (['idx', 'depth'], salinities),
            'pressure': (['idx', 'depth'], pressures)
        },
        coords={
            'idx': idxs,
            'depth': depths,
            'latitude': ('idx', lats),
            'longitude': ('idx', lons),
            'time': ('idx', dates),
            'file_name': ('idx', file_names)
        }
    )

    return ds
