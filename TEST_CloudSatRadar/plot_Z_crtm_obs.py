
import numpy as np
import netCDF4
import matplotlib.pyplot as pp


def read_ClousSat_Z( inputFile):

    # Time window when the Radar fly over the cyclone.
    timeRange = [17.34,17.42]

    csat = {}
    with netCDF4.Dataset( inputFile) as ds:

        t   = ds['Profile_time'][:]            #in seconds, [nray]
        t0  = ds['UTC_start'][:]               #in seconds, [nray]            
        hr  = (t+t0)/3600  #sec to hour
        win = (hr > timeRange[0]) & ( hr < timeRange[1]) 
        csat['time'] = hr[win]

        csat['Z']      = ds['Radar_Reflectivity'][win,:]      #in dBZ*100, [nray,nbin]
        csat['ht']     = ds['Height'][win,:]                  #in meters,  [nray,nbin]
        csat['binSfc'] = ds['SurfaceHeightBin'][win]          #            [nray] 
        csat['lat']    = ds['latitude'][win]                  #in degree,  [nray]            
        csat['lon']    = ds['longitude'][win]                 #in dgreee,  [nray]
    
        csat['Z']  = csat['Z'].transpose()*0.01 #[nbin,nray]
        Z3 = csat['Z']
        print(' Z ',len(Z3))
        print(' Zs ',Z3.shape)
        print(' M ',np.max(Z3), np.min(Z3))
        csat['ht'] = csat['ht'][0,:] #[nbin], All radar rays have the same vertical hight coordinates, keep one ray only.
        csat['ht'] = csat['ht']/1000 #m to km

        csat['Z']  = np.where( csat['Z']>100, np.nan, csat['Z'])
        #csat['Z']  = np.where( csat['Z']<-40, np.nan, csat['Z'])
        csat['Z']  = np.where( csat['Z']<-20, np.nan, csat['Z'])
    return csat


def read_CRTM_Z( crtmOutFile ):

    with open( crtmOutFile, 'rb') as dat:
        nchan   = np.fromfile( dat, 'i4', count=1).item()
        nprof   = np.fromfile( dat, 'i4', count=1).item()
        chList  = np.fromfile( dat, 'i4', count=nchan)
        print(' t1 ',nchan,nprof,chList)
        Zlist = []
        for ip in range(nprof):
            nlay = np.fromfile( dat, 'i4', 1).item()
            PZZ  = np.fromfile( dat, '4f8', nlay*nchan) #[nlay,nchan,3], the 3 are [press, height, Z, Zatten]
            
            PZZ[PZZ<-900] =np.nan
            Zlist.append( PZZ.reshape(nchan,nlay,4))
            
    PZZ = np.squeeze( np.array(Zlist)) #[nprof,nlay,3]

    crtm = {}
    crtm['press']  = PZZ[1,:,0]
    crtm['ht']     = PZZ[1,:,1]
    crtm['Z']      = PZZ[:,:,2].transpose() #[nlay,nprof]
    crtm['Z_attn'] = PZZ[:,:,3].transpose()

    return crtm


def plot_Z():

    vmin, vmax = -20, 30 #image color range.
    #vmin, vmax = -40, 25 #image color range.
    fig,(ax1,ax2) = pp.subplots(2, figsize=(8,6))

    # Plot CRTM result
    #
    #crtmOutFile = './cpr_cloudsat.CLOUDSAT.OCEAN.FORWARD.result.dat'
    crtmOutFile = './cpr_ref.bin'
    crtm = read_CRTM_Z( crtmOutFile)
    z1 = crtm['Z_attn']
    x1 = range(1, z1.shape[1]+1) #profile numbers 
    #y1 = crtm['press']
    y1 = crtm['ht']
    print(' height ',np.max(y1),np.min(y1))
    p1 = ax1.pcolormesh(x1,y1,z1, vmin=vmin,vmax=vmax)
    #ax1.invert_yaxis()
    ax1.set_ylim([0,20])
    ax1.grid()
    ax1.set_xlabel('Profile No.')
    ax1.set_ylabel('Pressure (hpa)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title('CRTM simulation')
    cb = fig.colorbar( p1)
    cb.ax.set_ylabel('Attenuated Reflectivity (dBZ)')

    print(' CloudSAT observation ')
    # Plot CloudSat observation
    #
    csatFile = './2009231162830_17611_CS_2B-GEOPROF_GRANULE_P1_R05_E02_F00.nc'
    csat = read_ClousSat_Z( csatFile)
    z2 = csat['Z']
    x2 = range(1, z2.shape[1]+1) #profile numbers 
    y2 = csat['ht']
    
    p2 = ax2.pcolormesh(x2,y2,z2, vmin=vmin,vmax=vmax)
    ax2.set_ylim([0,20])
    ax2.grid()
    ax2.set_xlabel('Profile No.')
    ax2.set_ylabel('Hight (km)')
    ax2.set_title('CloudSat CPR observation')
    cb = fig.colorbar( p2)
    cb.ax.set_ylabel('Attenuated Reflectivity (dBZ)')

    fig.tight_layout()
    pp.savefig( 'Z.png',dpi=200)


if __name__=='__main__':
    plot_Z()
