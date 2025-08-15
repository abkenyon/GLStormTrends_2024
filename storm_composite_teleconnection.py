import numpy as np
import glob
import xarray as xr
import pandas as pd
from shapely.geometry import Point, Polygon
import sys
import scipy.stats as stats

# a function to check whether a given point is in the polygon or not (from storm_stats.py)
def in_or_out( lon_in, lat_in, poly_in):
        pt = Point(lon_in,lat_in)
        point_in = pt.within(poly_in)
        return point_in

def concat_data(edata,datebool):
    print(yr,mth,dy,hr)
    msl = edata.msl.data[datebool]
    t2m = edata.t2m.data[datebool]
    u10,v10 = edata.u10.data[datebool],edata.v10.data[datebool]
    tp  = edata.tp.data[datebool]
    d2m = edata.d2m.data[datebool]
    tcwv= edata.tcwv.data[datebool] # total column vertical integrated water content (kg m^-2)


    # find the index of ERA5 lat array and lon array containing center of desired storm point
    elat,elon = edata.latitude.data,edata.longitude.data
    sloc_lat  = np.argwhere(elat==storms[ed]['lat'][flags][min_idx])[0][0]
    sloc_lon  = np.argwhere(elon==storms[ed]['lon'][flags][min_idx])[0][0]

    # grab variables in 550km^2 square centered on desired storm location
    lati = 20  #20 grid points ~= 550 km in latitude direction
    loni = 26  #26 grid points ~= 550 km in longitude direction
    msl_c = msl[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    t2m_c = t2m[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    u10_c = u10[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    v10_c = v10[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    tp_c  = tp[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    d2m_c = d2m[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
    tcwv_c= tcwv[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]

    # create new dimension to concatenate each storm's data onto one array
    msl_c = np.expand_dims(msl_c,0)
    t2m_c = np.expand_dims(t2m_c,0)
    u10_c = np.expand_dims(u10_c,0)
    v10_c = np.expand_dims(v10_c,0)
    tp_c  = np.expand_dims(tp_c,0)
    d2m_c = np.expand_dims(d2m_c,0)
    tcwv_c = np.expand_dims(tcwv_c,0)

    return msl_c,t2m_c,u10_c,v10_c,tp_c,d2m_c,tcwv_c,lati,loni

#################################
# Choose teleconnection
telecon = 'ao'

tcdict = {'ao':'AO','nao':'NAO','npgo':'NPGO','enso':'ENSO','pna':'PNA','pdo':'PDO'}
####################################

turbodir = '/nfs/turbo/seas-hutsona/shansen/sh_outputs/'

if telecon=='ao':
    enso = pd.read_csv(turbodir+f'{telecon}.csv',delim_whitespace=True)
else:
    enso = pd.read_csv(turbodir + f'{telecon}_cold.csv')
enso_std = enso.std()
enso_avg = enso.mean()
enso_pos_thres = enso_avg + enso_std
enso_neg_thres = enso_avg - enso_std
ensopos = enso_pos_thres['Index']
ensoneg = enso_neg_thres['Index']
poscutoff = ensopos
negcutoff = ensoneg
enso = enso.set_index('Year')

# Loop over all cyclogenesis locations, including 'False', which will create a composite for all cyclones
for cyclogen in ['False','greatlakes','alberta','colorado']:
    
    stable = pd.read_csv(turbodir+'stormtable.csv')
    stable = stable[(stable['year']>=1970)]
    locs = stable['loc'].to_numpy()
    phase= stable[f'{tcdict[telecon]} Phase'].to_numpy()

    neg_flag = 0
    pos_flag = 0
    storm_idx = -1
    pcount,ncount = 0,0
    for year in np.arange(1970,2022):# I CHANGED THIS FOR AO TESTING!!!!!!! ORIGINAL: 1960,2022
        print(year)
        # Load storm data
        etcdir  ='/nfs/turbo/seas-hutsona/shansen/etc_tracks/'
        fname=f'storm_track_slp_{year}.npz'
        data = np.load(etcdir+fname,allow_pickle=True,encoding='latin1')
        storms = data['storms']

        # define polygon for Great Lakes storms (from storm_stats.py)
        lats_polygon = [ 50, 50, 41, 41 ]
        lons_polygon = [ -75.5, -93.5, -93.5, -75.5 ]
        xy = zip(lons_polygon,lats_polygon)
        poly = Polygon(xy)

        # loop over tracked storms
        for ed in range(len(storms)):

            flagog = False # still need this "OG flag" to get through criteria below
            flags= [] # keeps track of which storm points are within great lakes region

            # Use year/month from start of storm track to define teleconnection phase
            #yr_tele,mn_tele = storms[ed]['year'][0],storms[ed]['month'][0]

            #tvalue = enso[enso['Month'] == mn_tele].loc[yr_tele]['Index']

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
                storm_idx = storm_idx+1

                # if teleconnection value doesn't meet cutoff value,
                # it's not considered strong enough to be "in phase", so we don't count the storm
                # absolute value!!!
#                 if float(tvalue) >= negcutoff and float(tvalue) <= poscutoff:
                if (phase[storm_idx]=='Neutral'):
                    continue
            
                if cyclogen is not False:
                    print(phase[storm_idx],locs[storm_idx],cyclogen,locs[storm_idx] == cyclogen)
                    if locs[storm_idx] == cyclogen:
                        print(storm_idx,locs[storm_idx])
                    else:
                        continue

                # Load ERA Data
                edir = '/nfs/turbo/seas-hutsona/shansen/era5_ETC_data/'
                ename= f'era5_{year-1}-{year}.nc'
                edata= xr.open_dataset(edir+ename)
                datetime = pd.to_datetime(edata.time)

                # find the date/time at which storm *in* the GLR is at minimum pressure
                min_idx = np.argmin(storms[ed]['amp'][flags])

                # save pressure, temperature, wind, and precip at date/time of storm with minimum pressure
                yr,mth = storms[ed]['year'][flags][min_idx],storms[ed]['month'][flags][min_idx]
                dy,hr  = storms[ed]['day'][flags][min_idx],storms[ed]['hour'][flags][min_idx]
                datebool     = (datetime.year==yr)&(datetime.month==mth)&(datetime.day==dy)&(datetime.hour==hr)

                msl_c,t2m_c,u10_c,v10_c,tp_c,d2m_c,tcwv_c,lati,loni = concat_data(edata,datebool)
                
                # concatenate new storm on previous storms' arrays
                if (pos_flag==1)&(phase[storm_idx] == 'Positive'):
                    msl_p = np.concatenate((msl_p,msl_c),axis=0)
                    t2m_p = np.concatenate((t2m_p,t2m_c),axis=0)
                    u10_p = np.concatenate((u10_p,u10_c),axis=0)
                    v10_p = np.concatenate((v10_p,v10_c),axis=0)
                    tp_p  = np.concatenate((tp_p,tp_c),axis=0)
                    d2m_p = np.concatenate((d2m_p,d2m_c),axis=0)
                    tcwv_p= np.concatenate((tcwv_p,tcwv_c),axis=0)
                    pcount = pcount+1

                if (neg_flag==1)&(phase[storm_idx] == 'Negative'):
                    msl_n = np.concatenate((msl_n,msl_c),axis=0)
                    t2m_n = np.concatenate((t2m_n,t2m_c),axis=0)
                    u10_n = np.concatenate((u10_n,u10_c),axis=0)
                    v10_n = np.concatenate((v10_n,v10_c),axis=0)
                    tp_n  = np.concatenate((tp_n,tp_c),axis=0)
                    d2m_n = np.concatenate((d2m_n,d2m_c),axis=0)
                    tcwv_n= np.concatenate((tcwv_n,tcwv_c),axis=0)
                    ncount = ncount+1

                # Initialize arrays for first storm
                if (pos_flag==0):
                    if phase[storm_idx] == 'Positive':
                        msl_p = msl_c
                        t2m_p = t2m_c
                        u10_p = u10_c
                        v10_p = v10_c
                        tp_p  = tp_c
                        d2m_p = d2m_c
                        tcwv_p= tcwv_c
                        pcount = 1
                        pos_flag=1 #indicates positive arrays have initialized

                if (neg_flag == 0):
                    if phase[storm_idx] == 'Negative':
                        msl_n = msl_c
                        t2m_n = t2m_c
                        u10_n = u10_c
                        v10_n = v10_c
                        tp_n  = tp_c
                        d2m_n = d2m_c
                        tcwv_n= tcwv_c
                        ncount = 1
                        neg_flag=1 #indicates negative arrays have initialized
                     
                print(phase[storm_idx])
                print('Positive Count',pcount,'Negative Count',ncount)
                    

        # Save averaged storm composite data for whole season
    # Calculate p-value from two-sample t-test
    for i in np.arange(msl_p.shape[1]):
        for j in np.arange(msl_p.shape[2]):
            pv_msl = stats.ttest_ind(msl_p,msl_n).pvalue
            pv_t2m = stats.ttest_ind(t2m_p,t2m_n).pvalue
            pv_tp  = stats.mannwhitneyu(tp_p,tp_n).pvalue
            pv_d2m = stats.ttest_ind(d2m_p,d2m_n).pvalue
            pv_wv  = stats.ttest_ind(tcwv_p,tcwv_n).pvalue

    # coordinates centered on middle of storm (in meters; negative is west/south, positive is north/east)
    west_east  = np.arange(-loni/4,(loni+1)/4,0.25)*85
    south_north = np.arange(-lati/4,(lati+1)/4,0.25)*111
    dim_dict_p  = dict(south_north=south_north,west_east=west_east,pcount=pcount)
    data_dict_pos = dict(t2m=(["south_north", "west_east"], np.nanmean(t2m_p,axis=0)[::-1,:],{'units':'K'}),
                     msl=(["south_north", "west_east"], np.nanmean(msl_p,axis=0)[::-1,:],{'units':'Pa'}),
                     u10=(["south_north", "west_east"], np.nanmean(u10_p,axis=0)[::-1,:],{'units':'m s**-1'}),
                     v10=(["south_north", "west_east"], np.nanmean(v10_p,axis=0)[::-1,:],{'units':'m s**-1'}),
                     tp =(["south_north", "west_east"], np.nanmean(tp_p,axis=0)[::-1,:],{'units':'m'}),
                     d2m=(["south_north", "west_east"], np.nanmean(d2m_p,axis=0)[::-1,:],{'units':'K'}),
                     tcwv=(["south_north", "west_east"], np.nanmean(tcwv_p,axis=0)[::-1,:],{'units':'kg m**-2'}),
                     t2m_std=(["south_north", "west_east"], np.nanstd(t2m_p,axis=0)[::-1,:],{'units':'K'}),
                     msl_std=(["south_north", "west_east"], np.nanstd(msl_p,axis=0)[::-1,:],{'units':'Pa'}),
                     u10_std=(["south_north", "west_east"], np.nanstd(u10_p,axis=0)[::-1,:],{'units':'m s**-1'}),
                     v10_std=(["south_north", "west_east"], np.nanstd(v10_p,axis=0)[::-1,:],{'units':'m s**-1'}),
                     tp_std =(["south_north", "west_east"], np.nanstd(tp_p,axis=0)[::-1,:],{'units':'m'}),
                     d2m_std=(["south_north", "west_east"], np.nanstd(d2m_p,axis=0)[::-1,:],{'units':'K'}),
                     tcwv_std=(["south_north", "west_east"], np.nanstd(tcwv_p,axis=0)[::-1,:],{'units':'kg m**-2'}),
                        )
    ds = xr.Dataset(data_vars=data_dict_pos, coords=dim_dict_p)
    outname=f'composite_{telecon}_pos_stddev_1970.nc'
    if cyclogen is not False:
        outname = outname.split('.')[0]+f'_{cyclogen}'+'_1970.nc'
    ds.to_netcdf(turbodir+outname)
    
    
    dim_dict_n  = dict(south_north=south_north,west_east=west_east,ncount=ncount)
    data_dict_neg = dict(t2m=(["south_north", "west_east"], np.nanmean(t2m_n,axis=0)[::-1,:],{'units':'K'}),
                     msl=(["south_north", "west_east"], np.nanmean(msl_n,axis=0)[::-1,:],{'units':'Pa'}),
                     u10=(["south_north", "west_east"], np.nanmean(u10_n,axis=0)[::-1,:],{'units':'m s**-1'}),
                     v10=(["south_north", "west_east"], np.nanmean(v10_n,axis=0)[::-1,:],{'units':'m s**-1'}),
                     tp =(["south_north", "west_east"], np.nanmean(tp_n,axis=0)[::-1,:],{'units':'m'}),
                     d2m=(["south_north", "west_east"], np.nanmean(d2m_n,axis=0)[::-1,:],{'units':'K'}),
                     tcwv=(["south_north", "west_east"], np.nanmean(tcwv_n,axis=0)[::-1,:],{'units':'kg m**-2'}),
                     t2m_std=(["south_north", "west_east"], np.nanstd(t2m_n,axis=0)[::-1,:],{'units':'K'}),
                     msl_std=(["south_north", "west_east"], np.nanstd(msl_n,axis=0)[::-1,:],{'units':'Pa'}),
                     u10_std=(["south_north", "west_east"], np.nanstd(u10_n,axis=0)[::-1,:],{'units':'m s**-1'}),
                     v10_std=(["south_north", "west_east"], np.nanstd(v10_n,axis=0)[::-1,:],{'units':'m s**-1'}),
                     tp_std =(["south_north", "west_east"], np.nanstd(tp_n,axis=0)[::-1,:],{'units':'m'}),
                     d2m_std=(["south_north", "west_east"], np.nanstd(d2m_n,axis=0)[::-1,:],{'units':'K'}),
                     tcwv_std=(["south_north", "west_east"], np.nanstd(tcwv_n,axis=0)[::-1,:],{'units':'kg m**-2'}),
                        )
    ds = xr.Dataset(data_vars=data_dict_neg, coords=dim_dict_n)
    #outname=f'composite_{telecon}_neg_stddev.nc'
    outname=f'composite_{telecon}_neg_stddev_1970.nc'
    if cyclogen is not False:
        outname = outname.split('.')[0]+f'_{cyclogen}'+'_1970.nc'
    ds.to_netcdf(turbodir+outname)
    
    
    dim_dict = dict(south_north=south_north,west_east=west_east)
    data_dict_pval = dict(t2m=(["south_north", "west_east"], pv_t2m[::-1,:],{'units':'K'}),
                     msl=(["south_north", "west_east"], pv_msl[::-1,:],{'units':'Pa'}),
                     tp =(["south_north", "west_east"], pv_tp[::-1,:],{'units':'m'}),
                     d2m=(["south_north", "west_east"], pv_d2m[::-1,:],{'units':'K'}),
                     tcwv=(["south_north", "west_east"], pv_wv[::-1,:],{'units':'kg m**-2'}),
                    )
    ds = xr.Dataset(data_vars=data_dict_pval, coords=dim_dict)
    outname=f'composite_{telecon}_pval_1970.nc'
    if cyclogen is not False:
        outname = outname.split('.')[0]+f'_{cyclogen}'+'_1970.nc'
    ds.to_netcdf(turbodir+outname)
