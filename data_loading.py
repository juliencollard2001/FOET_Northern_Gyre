import numpy as np
import scipy.io
import xarray as xr
import os
import pandas as pd
from datetime import datetime, date, time, timedelta

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
    dates = np.empty(data.shape[1], dtype='datetime64[ns]')
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


def get_CTD(year : int) -> xr.Dataset:
    """
    Load the CTD data for the year 2021.
    """
    base_folder = 'data/moose-cruises'
    file = base_folder + f'/CTD{year}_MOOSE_GE.mat'
    data = scipy.io.loadmat(file)['all_cnv']

    if year == 2021:

        idxs = np.arange(data.shape[1])
        depths = data[0,0][8][:,0]

        lats = np.empty(data.shape[1])
        lons = np.empty(data.shape[1])
        dates = np.empty(data.shape[1], dtype='datetime64[ns]')
        temperatures = np.empty((data.shape[1], len(depths)))
        salinities = np.empty((data.shape[1], len(depths)))
        pressures = np.empty((data.shape[1], len(depths)))

        for i in range(data.shape[1]):
            dates[i] = np.datetime64(datetime.strptime(data[0,i][3][0], '%d-%b-%Y %H:%M:%S'))
            lats[i] = float(data[0,i][6][0])
            lons[i] = float(data[0,i][5][0])
            temperatures[i] = data[0,i][9][:,0]
            salinities[i] = data[0,i][10][:,0]
            pressures[i] = data[0,i][7][:,0]

    else:

        idxs = np.arange(data.shape[1])
        depths = get_2021_CTD_data().depth.values

        station_data = scipy.io.loadmat(file)['all_stations']
        data_all = scipy.io.loadmat(file)

        lats = np.empty(data.shape[1])
        lons = np.empty(data.shape[1])
        dates = np.empty(data.shape[1], dtype='datetime64[ns]')
        temperatures = np.empty((data.shape[1], len(depths)))
        salinities = np.empty((data.shape[1], len(depths)))
        
        for i in range(data.shape[1]):

            time_like_number = data_all['sta_time'][0][i] # nb of days since 1st january 0000
            if not np.isnan(time_like_number):
                dates[i] = np.datetime64(datetime(1, 1, 1) + timedelta(days=time_like_number - 367))
            else:
                dates[i] = None
            lats[i] = data_all['sta_lat'][0][i]
            lons[i] = data_all['sta_lon'][0][i]
            temperatures[i] = data[0,i][1][:,0]
            salinities[i] = data[0,i][2][:,0]

    ds = xr.Dataset(
        {
            'temperature': (['time', 'depth'], temperatures),
            'salinity': (['time', 'depth'], salinities),
        },
        coords={
            'time': dates,
            'depth': -depths,
            'idx': ('time', idxs),
            'latitude': ('time', lats),
            'longitude': ('time', lons),
        }
    )

    

    return ds

def get_CTD_all_years() -> xr.Dataset:
    """
    Load the CTD data for all available years.
    """
    years = list(range(2010, 2020))
    ds = None
    for year in years:
        if ds is None:
            ds = get_CTD(year)

        else:
            ds = xr.concat([ds, get_CTD(year)], dim='time')
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
    LADCPsta = np.array(LADCP_sta, dtype=int).flatten()

    LADCP_time = LADCP['time']
    LADCP__time = np.array(LADCP_time).flatten()
    LADCP___time = np.concatenate([arr.flatten() for arr in LADCP__time])
    LADCPtime = pd.to_datetime(LADCP___time, format='%d-%b-%Y %H:%M:%S')

    # Create an xarray Dataset
    ds = xr.Dataset(
        {
            'U': (['time','depth'], LADCPU),
            'V': (['time','depth'], LADCPV),
            #'VELerror': (['station','depth'], LADCPVelerror),
        },
        coords={
            'time': LADCPtime,
            'depth': np.arange(500)*8.,

            'latitude': ('time', LADCPlat),
            'longitude': ('time', LADCPlon),
            'idx': ('time', LADCPsta),
            'leg': ('time', LADCPleg),
        },
        attrs={'year': year, 'source': 'MOOSE cruises'}
    )
    
    return ds

def get_LADCP_all_years() -> xr.Dataset:
    """
    Load the LADCP data for all available years.
    """
    years = [2012, 2015, 2017, 2019]
    ds = None
    for year in years:
        if ds is None:
            ds = get_LADCP(year)
        else:
            ds = xr.concat([ds, get_LADCP(year)], dim='time')
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

    SADCPU = SADCP['U'][0].T
    print(SADCPU.shape)
   
    SADCPV = SADCP['V'][0].T
   
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
            'U': (['time','depth'], SADCPU),
            'V': (['time','depth'], SADCPV),
        },
        coords={
            'time': SADCP_time,
            'depth': SADCP_Z,        
            'latitude': (['time'], SADCP_lat),
            'longitude': (['time'], SADCP_lon),
            'idx': (['time'], Stat),
            'leg': (['time'], SADCP_leg),
        },

        attrs={'year': year, 'source': 'MOOSE cruises'}
    )
    
    return ds

def get_SADCP_all_years() -> xr.Dataset:
    """
    Load the SADCP data for all available years.
    """
    years = [2012, 2015, 2017, 2019, 2021]
    ds = None
    for year in years:
        if ds is None:
            ds = get_SADCP(year)
        else:
            ds = xr.concat([ds, get_SADCP(year)], dim='time')
    return ds
