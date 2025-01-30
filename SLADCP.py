
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
    
    # get the files 
    file_path = f'data/moose-cruises/SADCP{year}_MOOSE_GE.mat'
    data = sio.loadmat(file_path)
    
    SADCP = data['cruise_SADCP'].flatten()  # Convert to 1D array

    # extract the different coordinates and variables

    SADCPZ = SADCP['Z']
    SADCP_Z = np.array(SADCPZ[0]).flatten()

    SADCPtime = SADCP['time'].flatten()
    SADCP_time1 = np.array(SADCPtime[0])
    SADCP_time = pd.to_datetime(SADCP_time1.flatten(), format='%d-%b-%Y %H:%M:%S')

    A = len(SADCP_time)     # to reshape the U and V arrays
    B = len(SADCP_Z)        # to reshape the U and V arrays

    SADCPU = SADCP['U']
    SADCP_U1 = np.array(SADCPU[0]).flatten()
    SADCP_U = np.reshape(SADCP_U1,(A,B))
   
    SADCPV = SADCP['V']
    SADCP_V1 = np.array(SADCPV[0]).flatten()
    SADCP_V = np.reshape(SADCP_V1,(A,B))
   
    SADCPlat = SADCP['lat']     
    SADCP_lat = np.array(SADCPlat[0]).flatten()

    SADCPlon = SADCP['lon']
    SADCP_lon = np.array(SADCPlon[0]).flatten()

    SADCPleg = SADCP['leg']
    SADCP_leg = np.array(SADCPleg[0]).flatten().astype(float)

    Stat = np.arange(1,len(SADCP_lon)+1)   # station number

    # creating a dataset 
    ds = xr.Dataset(
        {
            'U': (['station','depth'], SADCP_U),
            'V': (['station','depth'], SADCP_V),
        },
        coords={
            'station': Stat,
            'depth': SADCP_Z,        
            'latitude': (['station'], SADCP_lat),
            'longitude': (['station'], SADCP_lon),
            'time': (['station'], SADCP_time),
            'leg': (['station'], SADCP_leg),
        },

        attrs={'year': year, 'source': 'MOOSE cruises'}
    )
    
    return ds

data = get_SADCP(2017)
print(data)
