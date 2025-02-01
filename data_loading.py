import numpy as np
import scipy.io
import xarray as xr
import os
import pandas as pd
from datetime import datetime, date, time

import warnings
warnings.filterwarnings("ignore", message="Conversion of an array with ndim > 0 to a scalar is deprecated", category=DeprecationWarning)


def get_altimetry_data(years : list[int] , months : list[int], verbose=False) -> xr.Dataset:
    """
    Load altimetry data from the MOOSE-Altimetry dataset for the specified years and months.
    The data is loaded from the specified folders.
    """

    base_folder = 'data/MOOSE-Altimetry'
    ds = None
    folders = [base_folder + '/' + str(year) + '/' + str(month).zfill(2) for year in years for month in months]
    for folder in folders:
        if verbose:
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
                    if verbose:
                        print(f'Error: {e}')
                    continue
    if verbose:
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




def get_LADCP(year):

    # Load the data
    file_path = f'data/moose-cruises/LADCP{year}_MOOSE_GE.mat'
    data = scipy.io.loadmat(file_path)
    LADCP = data['all_LADCP'].flatten()  # Convert to 1D array

    # extract the different coordinates and variables

    # initialize the 2D arrays
    LADCPUs = LADCP['U']
    LADCPU = np.zeros((len(LADCPUs), len(LADCPUs[0].flatten())))  
    
    LADCPVs = LADCP['V']
    LADCPV = np.zeros((len(LADCPVs), len(LADCPVs[0].flatten())))    

    LADCPZs = LADCP['Z']
    LADCPZ = np.zeros((len(LADCPZs), len(LADCPZs[0].flatten())))   

    for i in range(len(LADCPUs)):
        LADCP_U = np.array(LADCPUs[i]).flatten()
        LADCPU[i, :] = LADCP_U

        LADCP_V = np.array(LADCPVs[i]).flatten()
        LADCPV[i, :] = LADCP_V

        LADCP_Z = np.array(LADCPZs[i]).flatten()
        LADCPZ[i, :] = LADCP_Z


    #LADCPVelerrors = LADCP['VELerror']
    #LADCPVelerror = np.zeros((len(LADCPVelerrors), len(LADCPVelerrors[0].flatten())))
    #for i in range(len(LADCPVelerrors)):
    #    LADCP_Velerror = np.array(LADCPVelerrors[i]).flatten()
    #    LADCPVelerror[i, :] = LADCP_Velerror

    # Initialize the 1D arrays

    LADCP_lat = LADCP['lat']
    LADCP__lat = np.array(LADCP_lat).flatten()
    LADCPlat = np.concatenate([arr.flatten() for arr in LADCP__lat])
                        
    LADCP_lon = LADCP['lon']
    LADCP__lon = np.array(LADCP_lon).flatten()
    LADCPlon = np.concatenate([arr.flatten() for arr in LADCP__lon])

 
    LADCP_leg = LADCP['leg']
    LADCPleg = np.array(LADCP_leg, dtype=float).flatten()

    LADCP_sta = LADCP['sta']
    LADCPsta = np.array(LADCP_sta, dtype=float).flatten()

    LADCP_time = LADCP['time']
    LADCP__time = np.array(LADCP_time).flatten()
    LADCP___time = np.concatenate([arr.flatten() for arr in LADCP__time])
    LADCPtime = pd.to_datetime(LADCP___time, format='%d-%b-%Y %H:%M:%S')

    # Create an xarray Dataset
    ds = xr.Dataset(
        {
            'U': (['station','depth'], LADCPU),
            'V': (['station','depth'], LADCPV),
            #'VELerror': (['station','depth'], LADCPVelerror),
        },
        coords={
            'station': LADCPsta,
            'depth': (['station','depth'],LADCPZ),

            'latitude': (['station'], LADCPlat),
            'longitude': (['station'], LADCPlon),
            'time': (['station'], LADCPtime),
            
            'leg': (['station'], LADCPleg),
        },
        attrs={'year': year, 'source': 'MOOSE cruises'}
    )
    
    return ds


def get_SADCP(year):
    
    # get the files 
    file_path = f'data/moose-cruises/SADCP{year}_MOOSE_GE.mat'
    data = scipy.io.loadmat(file_path)
    
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
