import numpy as np
import glob
import xarray as xr
import pandas as pd
from shapely.geometry import Point, Polygon
#from storm_functions import latlon2km

# a function to check whether a given point is in the polygon or not (from storm_stats.py)
def in_or_out( lon_in, lat_in, poly_in):
        pt = Point(lon_in,lat_in)
        point_in = pt.within(poly_in)
        return point_in

for year in np.arange(1960,2022):
    print(year)
    # Load storm data
    sdir  ='/nfs/turbo/seas-hutsona/shansen/etc_tracks/'
    fname=f'storm_track_slp_{year}.npz'
    data = np.load(sdir+fname,allow_pickle=True,encoding='latin1')
    storms = data['storms']

    # Load ERA Data
    edir = '/nfs/turbo/seas-hutsona/shansen/era5_ETC_data/'
    ename= f'era5_{year-1}-{year}.nc'
    edata= xr.open_dataset(edir+ename)
    datetime = pd.to_datetime(edata.time)

    # define polygon for Great Lakes storms (from storm_stats.py)
    lats_polygon = [ 50, 50, 41, 41 ]
    lons_polygon = [ -75.5, -93.5, -93.5, -75.5 ]
    xy = zip(lons_polygon,lats_polygon)
    poly = Polygon(xy)

    count = 1 # required to initialize concatenated arrays

    print('initiated composite for' + str(year))
        
    # loop over tracked storms
    for ed in range(len(storms)):
        flagog = False # still need this "OG flag" to get through criteria below
        flags= [] # keeps track of which storm points are within great lakes region

        # loop over data points on each storm track
        for nl in range(len(storms[ed]['lon'])):
            # flag will remain false if storm point is not in GL at this time
            flag = False
            
            # judge if a point is in the polygon by calling the in_or_out function
            if ( in_or_out(storms[ed]['lon'][nl],storms[ed]['lat'][nl], poly) ):
                flag  = True
                flagog= True
            
            flags.append(flag)

        # do the following code only for storm tracks that pass through GLR
        if (flagog):
 
            # find the date/time at which storm *in* the GLR is at minimum pressure
            min_idx = np.argmin(storms[ed]['amp'][flags])
           
            # save pressure, temperature, wind, and precip at date/time of storm with minimum pressure
            yr,mth = storms[ed]['year'][flags][min_idx],storms[ed]['month'][flags][min_idx]
            dy,hr  = storms[ed]['day'][flags][min_idx],storms[ed]['hour'][flags][min_idx]
            datebool     = (datetime.year==yr)&(datetime.month==mth)&(datetime.day==dy)&(datetime.hour==hr)
          
            print(yr,mth,dy,hr)
            dtes = []
            dtes_init = datetime[datebool]  
            dtes.extend(dtes_init)
            msl = edata.msl.data[datebool]
            t2m = edata.t2m.data[datebool]
            tp  = edata.tp.data[datebool]
            d2m = edata.d2m.data[datebool]
            tcwv= edata.tcwv.data[datebool] # total column vertical integrated water content (kg m^-2)
 

            # find the index of ERA5 lat array and lon array containing center of desired storm point
            elat,elon = edata.latitude.data,edata.longitude.data
            sloc_lat  = np.argwhere(elat==storms[ed]['lat'][flags][min_idx])[0][0]
            sloc_lon  = np.argwhere(elon==storms[ed]['lon'][flags][min_idx])[0][0]

            # grab variables in 550km^2 square centered on desired storm location
            #lati = 20  #20 grid points ~= 550 km in latitude direction 
            #loni = 26  #26 grid points ~= 550 km in longitude direction
            lati = int(storm_width/0.25)
            loni = int(storm_width/0.25)

            msl_c = msl[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            t2m_c = t2m[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            tp_c  = tp[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            d2m_c = d2m[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            tcwv_c= tcwv[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
           
            west_east  = np.arange(-loni/4,(loni+1)/4,0.25)*85
            south_north = np.arange(-lati/4,(lati+1)/4,0.25)*111
            
            msl_c = np.expand_dims(msl_c,0)
            t2m_c = np.expand_dims(t2m_c,0)
            tp_c  = np.expand_dims(tp_c,0)
            d2m_c = np.expand_dims(d2m_c,0)
            tcwv_c = np.expand_dims(tcwv_c,0)
            
            t2m_rev = t2m_c[::-1,:]
            msl_rev = msl_c[::-1,:]
            tp_rev  = tp_c[::-1,:]
            d2m_rev = d2m_c[::-1,:]
            tcwv_rev= tcwv_c[::-1,:]
            
            dim_dict = {'storm': np.arange(1)}
    
            dtes_array = np.array(dtes)
            dtes_data = xr.DataArray(dtes_array, dims=('storm',), coords={'storm': np.arange(1)})
        
            avg_t2m = np.nanmean(t2m_rev, axis=(1, 2)) #{'units':'K'}
            t2m_data = xr.DataArray(avg_t2m, dims=('storm',), coords={'storm': np.arange(1)})
    
            avg_msl= np.nanmean(msl_rev,axis=(1,2))    #{'units':'Pa'}
            msl_data = xr.DataArray(avg_msl, dims=('storm',), coords={'storm': np.arange(1)})
    
            avg_tp= np.nanmean(tp_rev,axis=(1,2))      #{'units':'m'})
            tp_data = xr.DataArray(avg_tp, dims=('storm',), coords={'storm': np.arange(1)})
    
            avg_d2m= np.nanmean(d2m_rev,axis=(1,2))    #{'units':'K'}
            d2m_data = xr.DataArray(avg_d2m, dims=('storm',), coords={'storm': np.arange(1)})
    
            avg_tcwv= np.nanmean(tcwv_rev,axis=(1,2))  #{'units':'kg m**-2'}
            tcwv_data = xr.DataArray(avg_tcwv, dims=('storm',), coords={'storm': np.arange(1)})
    
                     
            data_dict = {'Date': dtes_data, 't2m': t2m_data, 'msl': msl_data, 'tp': tp_data, 
                 'd2m': d2m_data, 'tcwv': tcwv_data} 
            #data_dict2 = dict(Date = (dtes_data), t2m= (t2m_data,{'units':'K'}), msl= (msl_data,{'units':'Pa'}),
                              #tp= (tp_data, {'units':'m'}), d2m= (d2m_data, {'units':'K'}),
                              #tcwv= (tcwv_data, {'units':'kg m**-2'})
            
            ds = xr.Dataset(data_vars=data_dict, coords = dim_dict)
            outname = f'storm_{storm_width}deg_{year-1}-{year}_{str(count)}.nc'
            turbodir = '/nfs/turbo/seas-hutsona/shansen/sh_outputs/1d_storms/'
            ds.to_netcdf(turbodir + outname)

            count = count + 1
