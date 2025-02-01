def ts_diagram(salt,temp,Pres):
    
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
    print(dens)
    # Plot data *********************************************
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
    #ax1.set_clim(0,10)
    #ax1.plot(salt,temp,'or',markersize=9)
    ax1.scatter(salt,temp,c=Pres,vmin=1, vmax=500)
    ax1.set_xlabel('Salinity')  
    ax1.set_ylabel('Temperature (C)')
    plt.show()

