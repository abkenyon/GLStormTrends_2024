'''

    Calculate storm stats (# of storms passing over the GL region)
    And plot storm tracks passing over the GL region
    
    Originally created by Eric J. Oliver for the stormTracking algorithm
    (https://github.com/ecjoliver/stormTracking)
  
    Modified by Abby Hutson and Ayumi Fujisaki-Manome for use with ERA5 
    and for storms in the Great Lakes Region.
    
    NOTE: the plotting code is very basic. It is used for temporary visualization of seasonal tracks. 

'''

# Last edited by Abby Hutson (2024-06-28)

# Load required modules

import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfea
from shapely.geometry import Point, Polygon
import sys
import cmocean
import matplotlib
from storm_functions import cartopy_gl_map
import pandas as pd

# a fuction to check whether a given point is in the polygon or not.
def in_or_out( lon_in, lat_in, poly_in):
        pt = Point(lon_in,lat_in)
        point_in = pt.within(poly_in)
        return point_in


# decide projection
#projection=ccrs.PlateCarree()

fig = plt.figure()
ax,datacrs = cartopy_gl_map(fig,1,1,0,leftlbls=[0],botlbls=[0])

# year input
count_all = 0
for yearid in np.arange(1960,2022):
    # yearid=sys.argv[1]

    # Load storm data
    track_type = 'bigdomain'
    #dir=f'/nfs/ayumif/projects/hutsona/stormTracking_GL/storm_tracks_{track_type}/'
    dir=f'/nfs/ayumif/projects/hutsona/gl_storm_tracking/bigdomain_150_24/'
    fname=dir+'storm_track_slp_'+str(yearid)+'.npz'
    data = np.load(fname,allow_pickle=True,encoding='latin1')
    storms = data['storms']
    #print storms

    # Plot storm tracks
 #   # Regional zoom in on the NW Atlantic
 #   fig=plt.figure()
 #   ax=fig.add_subplot(1,1,1,projection=projection)
 #   ax.set_extent([-120.0,-65.0,35.0,60.0],projection)
 #   ax.add_feature(cfea.COASTLINE,lw=.5)
 #   ax.add_feature(cfea.LAKES,alpha=0.3,lw=.5)

    # create polygon for mapping
    # this can be updated to use more sophisticated polygon (e.g., state shapefile)
    lats_polygon = [ 50, 50, 41, 41 ]
    lons_polygon = [ -75.5, -93.5, -93.5, -75.5 ]
    xy = zip(lons_polygon,lats_polygon)
    poly = Polygon(xy)

    # plot a polygon on a map
#    ax.add_geometries([poly], crs=datacrs, alpha=0.9)

    # initialize a counter for a number of storm tracks
    count = 0

    # print header
    print("\nperiod, count, min MSLP [Pa]")

    # loop over the tracked storms
    for ed in range(len(storms)):

        # reset a flag to judge if a track passes over the Polygon
        flag   = False
        glflags= []

        # select for: cyclonic storms 
        if (storms[ed]['type'] == 'cyclonic'):

            # look over data points on each storm track
            for nl in range(len(storms[ed]['lon'])):
                glflag = False
                # judge if a point is in the polygon by calling the in_or_out function
                if ( in_or_out(storms[ed]['lon'][nl],storms[ed]['lat'][nl], poly) ):
                    flag = True # the track passes over the polygon
                    glflag= True
                    
                glflags.append(glflag)

            # if the track passes over the Polygon, add to the count and plot the track on a map
            if ( flag ):
                count = count + 1
                count_all = count_all+1

                # find first recorded lat/lon, before points are removed at the boundaries
                start_lon = storms[ed]['lon'][0]
                start_lat = storms[ed]['lat'][0]
                start_month= storms[ed]['month'][0]
                start_day = storms[ed]['day'][0]
                start_hour= storms[ed]['hour'][0]
                
                # find time at which storm reaches minimum surface pressure while in GL polygon
                min_idx = np.argmin(storms[ed]['amp'][glflags])
                yr,mth  = storms[ed]['year'][glflags][min_idx],storms[ed]['month'][glflags][min_idx]
                dy,hr   = storms[ed]['day'][glflags][min_idx],storms[ed]['hour'][glflags][min_idx]
                min_msl = storms[ed]['amp'][glflags][min_idx]
                
                # find each storm's duration in the GL region
                durstart_yr = storms[ed]['year'][glflags][0]
                durstart_mth= storms[ed]['month'][glflags][0]
                durstart_dy = storms[ed]['day'][glflags][0]
                durstart_hr = storms[ed]['hour'][glflags][0]
                start_str   = pd.to_datetime(f'{durstart_yr}-{durstart_mth:02d}-{durstart_dy:02d} {durstart_hr:02d}:00:00')
                
                durend_yr = storms[ed]['year'][glflags][-1]
                durend_mth= storms[ed]['month'][glflags][-1]
                durend_dy = storms[ed]['day'][glflags][-1]
                durend_hr = storms[ed]['hour'][glflags][-1]
                end_str   = pd.to_datetime(f'{durend_yr}-{durend_mth:02d}-{durend_dy:02d} {durend_hr:02d}:00:00')
                
                duration = (end_str-start_str).total_seconds()/3600
                
                # assign a category based on start lat/lon
                alberta   = np.asarray([43.1,60,-120,-94.1]) #lat_s, lat_n, lon_w, lon_e
                colorado  = np.asarray([31,43,-110,-94.1])
                gulf      = np.asarray([25,36,-97,-85])
                greatlakes= np.asarray([36.5,51,-94,-75])
                atlantic  = np.asarray([25,40,-85,-72])
                origin_dict = ['alberta','colorado','gulf','greatlakes','atlantic']
                start_loc = 'N/A'
                for o,og in enumerate([alberta,colorado,gulf,greatlakes,atlantic]):
                    if ((start_lat<=og[1])&(start_lat>=og[0])) & ((start_lon<=og[-1])&(start_lon>=og[-2])):
                        start_loc = origin_dict[o]
                
                # drop when data points are close to the boundaries
                leftedge = 120 #120 or 103.5
                topedge  = 60 #50 or 56
                lons = np.ma.masked_where((storms[ed]['lon'][:]<-leftedge) | (storms[ed]['lat'][:] > topedge), storms[ed]['lon'][:])
                lats = np.ma.masked_where((storms[ed]['lon'][:]<-leftedge) | (storms[ed]['lat'][:] > topedge), storms[ed]['lat'][:])
                amps = np.ma.masked_where((storms[ed]['lon'][:]<-leftedge) | (storms[ed]['lat'][:] > topedge), storms[ed]['amp'][:])
                
                # print minimum mean sea level pressure (MSLP) for each storm that passes over the Polygon
                print(str(int(yearid)-1)+'Oct-'+str(yearid)+"Mar, ",str(count)+",",np.nanmean(amps))
                # save the above print statement in a text file
                if count_all==1:
                    tfile = open(dir+'storm_count_data.txt','w')
                else:
                    tfile = open(dir+'storm_count_data.txt','a')
                tfile.write(f'{yearid},{count},{np.nanmean(amps)},{np.nanmin(amps)} \n')
                tfile.close()
                
                if count_all==1:
                    cfile = open(dir+'storm_cyclogen_locs.txt','w')
                else:
                    cfile = open(dir+'storm_cyclogen_locs.txt','a')
                cfile.write(f'{yearid},{start_month},{start_day},{start_hour},{start_lon},{start_lat},{start_loc}\n')
                cfile.close()
                
                if count_all==1:
                    pfile = open(dir+'date_of_strongest.txt','w')
                else:
                    pfile = open(dir+'date_of_strongest.txt','a')
                pfile.write(f'{yr},{mth},{dy},{hr},{min_msl}\n')
                pfile.close()
                
                if count_all==1:
                    durfile = open(dir+'duration.txt','w')
                else:
                    durfile = open(dir+'duration.txt','a')
                durfile.write(f'{duration},{durstart_yr},{durstart_mth},{durstart_dy},{durstart_hr},{durend_yr},{durend_mth},{durend_dy},{durend_hr}\n')
                durfile.close()
                
                cmap = matplotlib.cm.get_cmap(cmocean.cm.matter)
                rgba = cmap(count_all/1000)
                ax.plot(lons,lats,'-', color=rgba,linewidth=1, alpha=0.6,transform=datacrs)

    #plt.title('Storm Tracks Crossing the Great Lakes Region ('+str(int(yearid)-1)+' Oct-'+str(yearid)+' Mar)')
#     plt.savefig(dir+f'figures/storm_tracks_GreatLakes_bigdomain', bbox_inches='tight', pad_inches=0.05, dpi=300)
    #plt.savefig(dir+f'figures/storm_tracks_GreatLakes_cleaned_{track_type}_'+str(yearid), bbox_inches='tight', pad_inches=0.05, dpi=300)

