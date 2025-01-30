
import numpy as np
import scipy.io as sio
import xarray as xr
import pandas as pd 
from datetime import datetime, date, time

def get_LADCP(year):
    file_path = f'data/moose-cruises/LADCP{year}_MOOSE_GE.mat'

    # Load the .mat file
    data = sio.loadmat(file_path)
    
    # Get LADCP data
    LADCP = data['all_LADCP'].flatten()  # Convert to 1D array
    
    # Columns to extract (in 1D)
    names_1d = ['dnum', 'lon', 'lat', 'U', 'V', 'Z', 'VELerror']
    
    # Extract and clean 1D data
    LADCP_data = {}
    for name in names_1d:
        inter = LADCP[name]  # Access the column in LADCP
        LADCP_data[name] = np.array([item[0, 0] for item in inter])  # Clean and store data
    
    # Columns to extract (in 2D, which will be flattened)
    names_2d = ['leg', 'time', 'sta']
    LADCP_data_2d = {}
    for name in names_2d:
        inter = LADCP[name]  # Access the column in LADCP
        LADCP_data_2d[name] = np.array([item[0, 0] if item.ndim == 2 else item for item in inter])  # Clean data

    # Flatten 2D columns and add to LADCP_data
    LADCP_data['leg'] = LADCP_data_2d['leg'].flatten().astype(float)
    LADCP_data['sta'] = LADCP_data_2d['sta'].flatten().astype(float)
    LADCP_data['time'] = pd.to_datetime(LADCP_data_2d['time'].flatten(), format='%d-%b-%Y %H:%M:%S')

    # Create an xarray Dataset
    ds = xr.Dataset(
        {
            'U': (['station'], LADCP_data['U']),
            'V': (['station'], LADCP_data['V']),
            'VELerror': (['station'], LADCP_data['VELerror']),
        },
        coords={
            'station': LADCP_data['sta'],
            'latitude': (['station'], LADCP_data['lat']),
            'longitude': (['station'], LADCP_data['lon']),
            'time': (['station'], LADCP_data['time']),
            'depth': (['station'], LADCP_data['Z']),
            'leg': (['station'], LADCP_data['leg']),
        },
        attrs={'year': year, 'source': 'MOOSE cruises'}
    )
    
    return ds




def get_SADCP(year):

    file_path = f'data/moose-cruises/SADCP{year}_MOOSE_GE.mat'

    # Load the .mat file
    data = sio.loadmat(file_path)
    
    # Get LADCP data
    SADCP = data['cruise_SADCP'].flatten()  # Convert to 1D array
     

    # Columns to extract (in 1D)
    names_1d = ['dnum', 'lon', 'lat', 'U', 'V', 'Z']
    
    # Create a dictionary to store cleaned data
    SADCP_data = {}
    
    # Extract and clean 1D data
    for name in names_1d:
        inter = SADCP[name]  # Access the column in LADCP
        SADCP_data[name] = np.array([item[0, 0] for item in inter])  # Clean and store data
    
    # Columns to extract (in 2D, which will be flattened)
    names_2d = ['time', 'leg']
    
    # Create a dictionary to store cleaned 2D data
    SADCP_data_2d = {}
    
    # Extract and clean 2D data
    for name in names_2d:
        inter = SADCP[name]  # Access the column in LADCP
        SADCP_data_2d[name] = np.array([item[0, 0] if item.ndim == 2 else item for item in inter])  # Clean data

    # Flatten 2D columns and add to LADCP_data
    SADCP_data['leg'] = SADCP_data_2d['leg'].flatten().astype(float)
    SADCP_data['time'] = pd.to_datetime(SADCP_data_2d['time'].flatten(), format='%d-%b-%Y %H:%M:%S')
   
     # Create an xarray Dataset
    ds = xr.Dataset(
        {
            'U': (['station'], SADCP_data['U']),
            'V': (['station'], SADCP_data['V']),
        },
        coords={
            'latitude': (['station'], SADCP_data['lat']),
            'longitude': (['station'], SADCP_data['lon']),
            'depth': (['station'], SADCP_data['Z']),
            'time': (['station'], SADCP_data['time']),
            'leg': (['station'], SADCP_data['leg']),
        },
        attrs={'year': year, 'source': 'MOOSE cruises'}
    )

    return SADCP

data = get_SADCP(2012)
print(data)