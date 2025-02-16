import numpy as np
import gsw 
import matplotlib.pyplot as plt
import cmocean.cm as cmo



def ts_diagram(salt,temp,depth, time):
    
    ## author : Yannis Cuypers 

    '''function for TS siagram with pressure(depth) in color'''
    #generating gridded values of salinity temp
    si=np.linspace(np.nanmin(salt),np.nanmax(salt))
    ti=np.linspace(np.nanmin(temp),np.nanmax(temp))
    
    
    # Calculate how many gridcells we need in the x and y dimensions
    xdim = len(si)
    ydim = len(ti)
    
    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
    
   
    
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            #dens[j,i]=sw.dens(si[i],ti[j],0)
            dens[j,i]=gsw.rho(si[i],ti[j],0)
    
    # Substract 1000 to convert to sigma-t
    dens = dens - 1000
    
    # Plot data *********************************************
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
    sc = ax1.scatter(salt,temp,c=depth, cmap = 'jet')
    cbar = plt.colorbar(sc, pad = 0.1, label = f"profondeur (m) pour le {time}")
    cbar.ax.invert_yaxis()
    plt.title(f"Diagramme T-S le {time}")
    ax1.set_xlabel(f'Salinity')  
    ax1.set_ylabel(f'Temperature (°C)')
    plt.show()

def ts_diagram_superposees(salt,temp,depth, salt2, temp2, depth2, time, time2):

    '''function for TS siagram with pressure(depth) in color'''
    #generating gridded values of salinity temp
    si=np.linspace(np.nanmin(salt),np.nanmax(salt))
    ti=np.linspace(np.nanmin(temp),np.nanmax(temp))

    si2=np.linspace(np.nanmin(salt2),np.nanmax(salt2))
    ti2=np.linspace(np.nanmin(temp2),np.nanmax(temp2))
    
    
    # Calculate how many gridcells we need in the x and y dimensions
    xdim = len(si)
    ydim = len(ti)

    xdim2 = len(si2)
    ydim2 = len(ti2)
    
    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
    dens2 = np.zeros((ydim2,xdim2))
    
   
    
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            #dens[j,i]=sw.dens(si[i],ti[j],0)
            dens[j,i]=gsw.rho(si[i],ti[j],0)
            dens2[j,i]=gsw.rho(si2[i],ti2[j],0)
    
    # Substract 1000 to convert to sigma-t
    dens = dens - 1000
    dens2 = dens2 - 1000
    
    # Plot data *********************************************
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')

    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
    sc = ax1.scatter(salt,temp,c=depth, marker = 'o', cmap = 'viridis')
    sc2 = ax1.scatter(salt2,temp2,c=depth2, marker = 'v', cmap = 'plasma')

    cbar = plt.colorbar(sc2, pad = 0.06, label = f"profondeur (m) pour le {time2}")
    cbar.ax.invert_yaxis()
    cbar2 = plt.colorbar(sc, pad = 0.1, label = f"profondeur (m) pour le {time}")
    cbar2.ax.invert_yaxis()

    plt.title(f"Diagramme T-S le {time} et le {time2}")
    ax1.set_xlabel(f'Salinity')  
    ax1.set_ylabel(f'Temperature (°C)')
    plt.show()

def ts_diagram_video(salt,temp,depth, time):
    

    '''function for TS siagram with pressure(depth) in color'''
    #generating gridded values of salinity temp
    si=np.linspace(np.nanmin(salt),np.nanmax(salt))
    ti=np.linspace(np.nanmin(temp),np.nanmax(temp))
    
    
    # Calculate how many gridcells we need in the x and y dimensions
    xdim = len(si)
    ydim = len(ti)
    
    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
    
   
    
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            #dens[j,i]=sw.dens(si[i],ti[j],0)
            dens[j,i]=gsw.rho(si[i],ti[j],0)
    
    # Substract 1000 to convert to sigma-t
    dens = dens - 1000
    
    # Plot data *********************************************
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
    sc = ax1.scatter(salt,temp,c=depth, cmap = 'jet')
    cbar = plt.colorbar(sc, pad = 0.1, label = f"profondeur (m) pour le {time}")
    cbar.ax.invert_yaxis()
    plt.title(f"Diagramme T-S le {time}")
    ax1.set_xlabel(f'Salinity')  
    ax1.set_ylabel(f'Temperature (°C)')
    #plt.show()
