##########################################################
# storm_composite.py is a script that creates the storm-centered 
# composites used in Hutson et al. (2024) (AGU GRL publication)
#
# Requires: ERA5 (or other gridded data from which to see meteorological data)
#           storm tracks resulting from the stormTracking algorithm (in a .npz structured file)
#
# Last edited by Abby Hutson (2024-06-28)
#
###########################################################

import numpy as np
import glob
import xarray as xr
import pandas as pd
from shapely.geometry import Point, Polygon
from storm_functions import latlon2km

# a function to check whether a given point is in the polygon or not (from storm_stats.py)
def in_or_out( lon_in, lat_in, poly_in):
        pt = Point(lon_in,lat_in)
        point_in = pt.within(poly_in)
        return point_in

#################################################################
# "True" means that we will create composites based on track length
track_length = False
storm_width  = 10
precip_intensity=True
#################################################################
    
for year in np.arange(1960,2022):
    # Load storm data
    dir  ='/nfs/ayumif/projects/hutsona/stormTracking_GL/bigdomain_150_24/'
    fname=f'storm_track_slp_{year}.npz'
    data = np.load(dir+fname,allow_pickle=True,encoding='latin1')
    storms = data['storms']

    # Load ERA Data
    edir = '/home/hutsona/data_big_domain/'
    ename= f'era5_{year-1}-{year}.nc'
    edata= xr.open_dataset(edir+ename)
    datetime = pd.to_datetime(edata.time)

    # define polygon for Great Lakes storms (from storm_stats.py)
    lats_polygon = [ 50, 50, 41, 41 ]
    lons_polygon = [ -75.5, -93.5, -93.5, -75.5 ]
    xy = zip(lons_polygon,lats_polygon)
    poly = Polygon(xy)

    count=0 # required to initialize concatenated arrays

    if track_length is True:
        long_tracks = []
        short_tracks= []
        short_flag = False
        long_flag  = False
        
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
            
            if track_length is True:
                # Calculate total track length
                d_tot = 0
                for k in range(len(storms[ed]['lon'])-1):
                    d_tot += latlon2km(storms[ed]['lon'][k], storms[ed]['lat'][k], storms[ed]['lon'][k+1], storms[ed]['lat'][k+1])
                if d_tot <= 2000:
                    short_tracks.append(d_tot)
                elif d_tot > 2000:
                    long_tracks.append(d_tot)
                    
            # find the date/time at which storm *in* the GLR is at minimum pressure
            min_idx = np.argmin(storms[ed]['amp'][flags])
           
            # save pressure, temperature, wind, and precip at date/time of storm with minimum pressure
            yr,mth = storms[ed]['year'][flags][min_idx],storms[ed]['month'][flags][min_idx]
            dy,hr  = storms[ed]['day'][flags][min_idx],storms[ed]['hour'][flags][min_idx]
            datebool     = (datetime.year==yr)&(datetime.month==mth)&(datetime.day==dy)&(datetime.hour==hr)
            if track_length is True:
                print(yr,mth,dy,hr,d_tot)
            else:
                print(yr,mth,dy,hr)
            msl = edata.msl.data[datebool]
            t2m = edata.t2m.data[datebool]
            u10,v10 = edata.u10.data[datebool],edata.v10.data[datebool]
            tp  = edata.tp.data[datebool]
            d2m = edata.d2m.data[datebool]
            e   = edata.e.data[datebool] # evaporation, in water equivalent meters
            ptype=edata.ptype.data[datebool] #precip type, units: code table (4.201)
            sf  = edata.sf.data[datebool] # snowfall, meters of water equivalent
            tcwv= edata.tcwv.data[datebool] # total column vertical integrated water content (kg m^-2)
            
            # averaging ptype doesn't make sense. So I will create mean rain and mean snow fields.
            pt   = np.rint(ptype)
            snow,rain,fzrn,mix = np.copy(tp),np.copy(tp),np.copy(tp),np.copy(tp)
            snow[(pt==8)|(pt==7)|(pt==3)|(pt==1)|(pt==0)] = 0.0
            rain[(pt==8)|(pt==7)|(pt==6)|(pt==5)|(pt==3)|(pt==0)] = 0.0 # np.nan creates unrealistic means...
            fzrn[(pt==8)|(pt==7)|(pt==6)|(pt==5)|(pt==1)|(pt==0)] = 0.0 # ...especially when comparing to total precip
            mix[(pt==8)|(pt==1)|(pt==6)|(pt==5)|(pt==3)|(pt==0)]  = 0.0

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
            u10_c = u10[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            v10_c = v10[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            tp_c  = tp[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            d2m_c = d2m[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            e_c   = e[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            sf_c  = sf[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            tcwv_c= tcwv[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            
            snow_c= snow[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            rain_c= rain[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            fzrn_c= fzrn[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]
            mix_c = mix[0,sloc_lat-lati:sloc_lat+(lati+1),sloc_lon-loni:sloc_lon+(loni+1)]

            # create new dimension to concatenate each storm's data onto one array
            msl_c = np.expand_dims(msl_c,0)
            t2m_c = np.expand_dims(t2m_c,0)
            u10_c = np.expand_dims(u10_c,0)
            v10_c = np.expand_dims(v10_c,0)
            tp_c  = np.expand_dims(tp_c,0)
            d2m_c = np.expand_dims(d2m_c,0)
            e_c   = np.expand_dims(e_c,0)
            sf_c  = np.expand_dims(sf_c,0)
            tcwv_c = np.expand_dims(tcwv_c,0)
            
            snow_c = np.expand_dims(snow_c,0)
            rain_c = np.expand_dims(rain_c,0)
            fzrn_c = np.expand_dims(fzrn_c,0)
            mix_c  = np.expand_dims(mix_c,0)
 
            # Initialize arrays for first storm
            if track_length is True:
                if (short_flag is False)&(d_tot <=2000):
                    msl_short = msl_c
                    t2m_short = t2m_c
                    u10_short = u10_c
                    v10_short = v10_c
                    tp_short  = tp_c
                    d2m_short = d2m_c
                    e_short   = e_c
                    sf_short  = sf_c
                    tcwv_short= tcwv_c
                    snow_short= snow_c
                    rain_short= rain_c
                    fzrn_short= fzrn_c
                    mix_short = mix_c
                    short_flag=True
                    
                if (long_flag is False)&(d_tot >2000):
                    msl_long = msl_c
                    t2m_long = t2m_c
                    u10_long = u10_c
                    v10_long = v10_c
                    tp_long  = tp_c
                    d2m_long = d2m_c
                    e_long   = e_c
                    sf_long  = sf_c
                    tcwv_long= tcwv_c
                    snow_long= snow_c
                    rain_long= rain_c
                    fzrn_long= fzrn_c
                    mix_long = mix_c
                    long_flag=True
            else:
                if count==0:             
                    msl_all = msl_c
                    t2m_all = t2m_c
                    u10_all = u10_c
                    v10_all = v10_c
                    tp_all  = tp_c
                    d2m_all = d2m_c
                    e_all   = e_c
                    sf_all  = sf_c
                    tcwv_all= tcwv_c
                    snow_all= snow_c
                    rain_all= rain_c
                    fzrn_all= fzrn_c
                    mix_all = mix_c
                
                newflag=1 #indicates new arrays have initialized

            # concatenate new storm on previous storms' arrays
            if ((track_length is False)&(count>0))|((track_length is True)):
                if track_length is True:
                    if (d_tot <= 2000)&(short_flag is True):
                        msl_short = np.concatenate((msl_short,msl_c),axis=0)
                        t2m_short = np.concatenate((t2m_short,t2m_c),axis=0)
                        u10_short = np.concatenate((u10_short,u10_c),axis=0)
                        v10_short = np.concatenate((v10_short,v10_c),axis=0)
                        tp_short  = np.concatenate((tp_short,tp_c),axis=0)
                        d2m_short = np.concatenate((d2m_short,d2m_c),axis=0)
                        e_short   = np.concatenate((e_short,e_c),axis=0)
                        sf_short  = np.concatenate((sf_short,sf_c),axis=0)
                        tcwv_short= np.concatenate((tcwv_short,tcwv_c),axis=0)
                        snow_short= np.concatenate((snow_short,snow_c),axis=0)
                        rain_short= np.concatenate((rain_short,rain_c),axis=0)
                        fzrn_short= np.concatenate((fzrn_short,fzrn_c),axis=0)
                        mix_short= np.concatenate((mix_short,mix_c),axis=0)
                    if (d_tot > 2000)&(long_flag is True):
                        msl_long = np.concatenate((msl_long,msl_c),axis=0)
                        t2m_long = np.concatenate((t2m_long,t2m_c),axis=0)
                        u10_long = np.concatenate((u10_long,u10_c),axis=0)
                        v10_long = np.concatenate((v10_long,v10_c),axis=0)
                        tp_long  = np.concatenate((tp_long,tp_c),axis=0)
                        d2m_long = np.concatenate((d2m_long,d2m_c),axis=0)
                        e_long   = np.concatenate((e_long,e_c),axis=0)
                        sf_long  = np.concatenate((sf_long,sf_c),axis=0)
                        tcwv_long= np.concatenate((tcwv_long,tcwv_c),axis=0)
                        snow_long= np.concatenate((snow_long,snow_c),axis=0)
                        rain_long= np.concatenate((rain_long,rain_c),axis=0)
                        fzrn_long= np.concatenate((fzrn_long,fzrn_c),axis=0)
                        mix_long = np.concatenate((mix_long,mix_c),axis=0)
                else:
                    msl_all = np.concatenate((msl_all,msl_c),axis=0)
                    t2m_all = np.concatenate((t2m_all,t2m_c),axis=0)
                    u10_all = np.concatenate((u10_all,u10_c),axis=0)
                    v10_all = np.concatenate((v10_all,v10_c),axis=0)
                    tp_all  = np.concatenate((tp_all,tp_c),axis=0)
                    d2m_all = np.concatenate((d2m_all,d2m_c),axis=0)
                    e_all   = np.concatenate((e_all,e_c),axis=0)
                    sf_all  = np.concatenate((sf_all,sf_c),axis=0)
                    tcwv_all= np.concatenate((tcwv_all,tcwv_c),axis=0)
                    snow_all= np.concatenate((snow_all,snow_c),axis=0)
                    rain_all= np.concatenate((rain_all,rain_c),axis=0)
                    fzrn_all= np.concatenate((fzrn_all,fzrn_c),axis=0)
                    mix_all= np.concatenate((mix_all,mix_c),axis=0)
 
            count = count+1

    
# I MOVED THIS CODE/PROCESS TO PRECIP_HISTOGRAMS.PY

#     if precip_intensity is True:
#         if track_length is False:
#             bins = np.asarray([.001,.002,.003,.004,.005,.006,.007,.008,.009,.010,.015,.020,.030,.1])
#             tp_intensity  = np.histogram(tp_all,bins=bins)[0]
#             snow_intensity= np.histogram(snow_all,bins=bins)[0]
#             rain_intensity= np.histogram(rain_all,bins=bins)[0]
            
#             if year==1960:
#                 pfile = open(dir+'precip_hist.txt','w')
#                 sfile = open(dir+'snow_hist.txt','w')
#                 rfile = open(dir+'rain_hist.txt','w')
#                 for tfile in [pfile,sfile,rfile]:
#                     tfile.write(f'{bins*1000} \n')
#             else:
#                 pfile = open(dir+'precip_hist.txt','a')
#                 sfile = open(dir+'snow_hist.txt','a')
#                 rfile = open(dir+'rain_hist.txt','a')
                
#             for i,ini in enumerate(tp_intensity):
#                 pfile.write(f'{tp_intensity[i]}, ')
#                 sfile.write(f'{snow_intensity[i]}, ')
#                 rfile.write(f'{rain_intensity[i]}, ')
#             pfile.write('\n')
#             sfile.write('\n')
#             rfile.write('\n')
#             pfile.close()
#             sfile.close()
#             rfile.close()
              
            
        if track_length is True:
            tp_intensity_l  = np.histogram(tp_long,bins=np.arange(1/1000,65/1000,5/1000))
            snow_intensity_l= np.histogram(snow_long,bins=np.arange(1/1000,65/1000,5/1000))
            rain_intensity_l= np.histogram(rain_long,bins=np.arange(1/1000,65/1000,5/1000))
            
            tp_intensity_s  = np.histogram(tp_short,bins=np.arange(1/1000,65/1000,5/1000))
            snow_intensity_s= np.histogram(snow_short,bins=np.arange(1/1000,65/1000,5/1000))
            rain_intensity_s= np.histogram(rain_short,bins=np.arange(1/1000,65/1000,5/1000))
            
        
    if precip_intensity is False:
        # Save averaged storm composite data for whole season

        # coordinates centered on middle of storm (in meters; negative is west/south, positive is north/east)
        west_east  = np.arange(-loni/4,(loni+1)/4,0.25)*85
        south_north = np.arange(-lati/4,(lati+1)/4,0.25)*111
        if track_length is False:
            dim_dict  = dict(south_north=south_north,west_east=west_east,count=count)
            data_dict = dict(t2m=(["south_north", "west_east"], np.nanmean(t2m_all,axis=0)[::-1,:],{'units':'K'}),
                             msl=(["south_north", "west_east"], np.nanmean(msl_all,axis=0)[::-1,:],{'units':'Pa'}),
                             u10=(["south_north", "west_east"], np.nanmean(u10_all,axis=0)[::-1,:],{'units':'m s**-1'}),
                             v10=(["south_north", "west_east"], np.nanmean(v10_all,axis=0)[::-1,:],{'units':'m s**-1'}),
                             tp =(["south_north", "west_east"], np.nanmean(tp_all,axis=0)[::-1,:],{'units':'m'}),
                             d2m=(["south_north", "west_east"], np.nanmean(d2m_all,axis=0)[::-1,:],{'units':'K'}),
                             e  =(["south_north", "west_east"], np.nanmean(e_all,axis=0)[::-1,:],{'units':'m'}),
                             sf  =(["south_north", "west_east"], np.nanmean(sf_all,axis=0)[::-1,:],{'units':'m'}),
                             tcwv=(["south_north", "west_east"], np.nanmean(tcwv_all,axis=0)[::-1,:],{'units':'kg m**-2'}),
                             snow=(["south_north", "west_east"], np.nanmean(snow_all,axis=0)[::-1,:],{'units':'m'}),
                             rain=(["south_north", "west_east"], np.nanmean(rain_all,axis=0)[::-1,:],{'units':'m'}),
                             fzrn=(["south_north", "west_east"], np.nanmean(fzrn_all,axis=0)[::-1,:],{'units':'m'}),
                             mix =(["south_north", "west_east"], np.nanmean(mix_all,axis=0)[::-1,:],{'units':'m'}),
                            )
            ds = xr.Dataset(data_vars=data_dict, coords=dim_dict)

        if track_length is True:
            dim_dict_s  = dict(south_north=south_north,west_east=west_east,count=len(short_tracks))
            data_dict_s = dict(t2m=(["south_north", "west_east"], np.nanmean(t2m_short,axis=0)[::-1,:],{'units':'K'}),
                             msl=(["south_north", "west_east"], np.nanmean(msl_short,axis=0)[::-1,:],{'units':'Pa'}),
                             u10=(["south_north", "west_east"], np.nanmean(u10_short,axis=0)[::-1,:],{'units':'m s**-1'}),
                             v10=(["south_north", "west_east"], np.nanmean(v10_short,axis=0)[::-1,:],{'units':'m s**-1'}),
                             tp =(["south_north", "west_east"], np.nanmean(tp_short,axis=0)[::-1,:],{'units':'m'}),
                             d2m=(["south_north", "west_east"], np.nanmean(d2m_short,axis=0)[::-1,:],{'units':'K'}),
                             e  =(["south_north", "west_east"], np.nanmean(e_short,axis=0)[::-1,:],{'units':'m'}),
                             sf  =(["south_north", "west_east"], np.nanmean(sf_short,axis=0)[::-1,:],{'units':'m'}),
                             tcwv=(["south_north", "west_east"], np.nanmean(tcwv_short,axis=0)[::-1,:],{'units':'kg m**-2'}),
                             snow=(["south_north", "west_east"], np.nanmean(snow_short,axis=0)[::-1,:],{'units':'m'}),
                             rain=(["south_north", "west_east"], np.nanmean(rain_short,axis=0)[::-1,:],{'units':'m'}),
                             fzrn=(["south_north", "west_east"], np.nanmean(fzrn_short,axis=0)[::-1,:],{'units':'m'}),
                             mix =(["south_north", "west_east"], np.nanmean(mix_short,axis=0)[::-1,:],{'units':'m'}),
                            )
            ds_s = xr.Dataset(data_vars=data_dict_s, coords=dim_dict_s)

            dim_dict_l  = dict(south_north=south_north,west_east=west_east,count=len(long_tracks))
            data_dict_l = dict(t2m=(["south_north", "west_east"], np.nanmean(t2m_long,axis=0)[::-1,:],{'units':'K'}),
                             msl=(["south_north", "west_east"], np.nanmean(msl_long,axis=0)[::-1,:],{'units':'Pa'}),
                             u10=(["south_north", "west_east"], np.nanmean(u10_long,axis=0)[::-1,:],{'units':'m s**-1'}),
                             v10=(["south_north", "west_east"], np.nanmean(v10_long,axis=0)[::-1,:],{'units':'m s**-1'}),
                             tp =(["south_north", "west_east"], np.nanmean(tp_long,axis=0)[::-1,:],{'units':'m'}),
                             d2m=(["south_north", "west_east"], np.nanmean(d2m_long,axis=0)[::-1,:],{'units':'K'}),
                             e  =(["south_north", "west_east"], np.nanmean(e_long,axis=0)[::-1,:],{'units':'m'}),
                             sf  =(["south_north", "west_east"], np.nanmean(sf_long,axis=0)[::-1,:],{'units':'m'}),
                             tcwv=(["south_north", "west_east"], np.nanmean(tcwv_long,axis=0)[::-1,:],{'units':'kg m**-2'}),
                             snow=(["south_north", "west_east"], np.nanmean(snow_long,axis=0)[::-1,:],{'units':'m'}),
                             rain=(["south_north", "west_east"], np.nanmean(rain_long,axis=0)[::-1,:],{'units':'m'}),
                             fzrn=(["south_north", "west_east"], np.nanmean(fzrn_long,axis=0)[::-1,:],{'units':'m'}),
                             mix =(["south_north", "west_east"], np.nanmean(mix_long,axis=0)[::-1,:],{'units':'m'}),
                            )
            ds_l = xr.Dataset(data_vars=data_dict_l, coords=dim_dict_l)

        if track_length is False:
            outname=f'storm_composite_{storm_width}deg_{year-1}-{year}.nc'
            ds.to_netcdf(dir+outname)
        if track_length is True:
            outname_l=f'storm_composite_{storm_width}deg_{year-1}-{year}_long.nc'
            outname_s=f'storm_composite_{storm_width}deg_{year-1}-{year}_short.nc'
            ds_l.to_netcdf(dir+outname_l)
            ds_s.to_netcdf(dir+outname_s)

