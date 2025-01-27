# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:51:32 2025

@author: georgiamcquade
"""
import scipy.io as sio # type: ignore
import numpy as np # type: ignore

def get_LADCP_data(year):

    file_path = f'data/moose-cruises/LADCP{year}_MOOSE_GE.mat'

    # Load the .mat file
    data = sio.loadmat(file_path)
    
    # Get LADCP data
    LADCP = data['all_LADCP'].flatten()  # Convert to 1D array
    
    # Columns to extract (in 1D)
    names_1d = ['dnum', 'lon', 'lat', 'U', 'V', 'Z', 'VELerror']
    
    # Create a dictionary to store cleaned data
    LADCP_data = {}
    
    # Extract and clean 1D data
    for name in names_1d:
        inter = LADCP[name]  # Access the column in LADCP
        LADCP_data[name] = np.array([item[0, 0] for item in inter])  # Clean and store data
    
    # Columns to extract (in 2D, which will be flattened)
    names_2d = ['leg', 'time', 'sta']
    
    # Create a dictionary to store cleaned 2D data
    LADCP_data_2d = {}
    
    # Extract and clean 2D data
    for name in names_2d:
        inter = LADCP[name]  # Access the column in LADCP
        LADCP_data_2d[name] = np.array([item[0, 0] if item.ndim == 2 else item for item in inter])  # Clean data

    # Flatten 2D columns and add to LADCP_data
    LADCP_data['leg'] = LADCP_data_2d['leg'].flatten()
    LADCP_data['sta'] = LADCP_data_2d['sta'].flatten()
    LADCP_data['time'] = LADCP_data_2d['time'].flatten()
    
        
    # noms : 
        #'leg','time', 'sta', 'dnum', 'lon', 'lat', 'U', 'V', 'Z', 'VELerror'
        
    #copier coller pratique :
        
    #years = [2012, 2015, 2017, 2019, 2021]
    #LADCP = {}
    #for year in years:
    #    LADCP[year] =  get_LADCP_data(year)

    return LADCP_data

#%%

def get_SADCP_data(year):

    file_path = f'data/moose-cruises/LADCP{year}_MOOSE_GE.mat'

    # Load the .mat file
    data = sio.loadmat(file_path)
    
    # Get LADCP data
    SADCP = data['all_LADCP'].flatten()  # Convert to 1D array
    
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
    SADCP_data['time'] = SADCP_data_2d['time'].flatten()
    
       
    # noms : 
        #'leg','time', 'dnum', 'lon', 'lat', 'U', 'V', 'Z', 
        
    #copier coller pratique :
        
    #years = [2012, 2015, 2017, 2019, 2021]
    #SADCP = {}
    #for year in years:
    #    SADCP[year] =  get_SADCP_data(year)

    return SADCP_data