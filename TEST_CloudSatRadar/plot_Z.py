
import numpy as np
import matplotlib.pyplot as pp

def plot_Z():

    crtmOutFile = './cpr_cloudsat.CLOUDSAT.OCEAN.FORWARD.result.dat'

    with open( crtmOutFile, 'rb') as dat:
        nchan   = np.fromfile( dat, 'i4', count=1).item()
        nprof   = np.fromfile( dat, 'i4', count=1).item()
        chList  = np.fromfile( dat, 'i4', count=nchan)

        Zlist = []
        for ip in range(nprof):
            nlay = np.fromfile( dat, 'i4', 1).item()
            PZZ  = np.fromfile( dat, '3f8', nlay*nchan) #[nlay,nchan,3], the 3 are [press, Z, Zatten]
            
            PZZ[PZZ<-900] =np.nan
            Zlist.append( PZZ.reshape(nchan,nlay,3))

    PZZ = np.squeeze( np.array(Zlist)) #[nprof,nlay,3]
    press = PZZ[1,:,0]
    Z     = PZZ[:,:,1].transpose() #[nlay,nprof]
    Zattn = PZZ[:,:,2].transpose()


    x = range(1,nprof+1)
    y = press
    pp.figure(figsize=(10,5))
    pp.pcolormesh(x,y, Z)
    pp.gca().invert_yaxis()
    pp.grid()
    pp.xlabel('Profile No.')
    pp.ylabel('Pressure (hpa)')
    cb = pp.colorbar()
    cb.ax.set_ylabel('Reflectivity (dBZ)')
    pp.savefig( 'Z.png',dpi=200)



if __name__=='__main__':
    plot_Z()