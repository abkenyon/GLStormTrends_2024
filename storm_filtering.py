'''

  Software used to filter anomalous points resulting from the storm_detection.py software
  
  Purpose: Remove any points erroneously identified as ETC storm centers within 
  a 3 degree lon/lat box surrounding the strongest cyclone point (lowest MSLP). 

'''

# Last edited by Abby Hutson (2024-06-28)

import numpy as np
import xarray as xr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
from itertools import repeat
import sys


# This script is used to filter out the storm-detection points created by
# storm_detection.py. It will find the points with the lowest MSLP associated 
# with different storms in the Great Lakes Region.  

# afm 20220616 
yearid=sys.argv[1]

# Open ERA5 Reference file for lat/lon grid
reffile = '/nfs/ayumif/projects/hutsona/era_grid_ref.nc'
era_fh  = xr.open_dataset(reffile)
elon,elat = era_fh.longitude.data,era_fh.latitude.data

########################################################
# Read in .npz file containing cyclonic and anticyclonic points
detfile   = './storm_det_slp_'+yearid+'.npz'
########################################################
storm_det = np.load(detfile,allow_pickle=True,encoding='latin1')
storms    = storm_det['storms']
day       = storm_det['day']
month     = storm_det['month']
hour      = storm_det['hour']
year      = storm_det['year']

# Maximum distance (in degrees) covered by one cyclone
max_lon_dist=7 # 7*0.25*2 = 3.5 degrees longitude
max_lat_dist=6 # 6*0.25*2 = 3.0 degrees latitude

# Set empty arrays to fill as we loop through storm points
years,months,days,hours = [],[],[],[]
storms_file = []

# Loop through storm points
for s in np.arange(0,len(storms)):
    storms_lat,storms_lon, amp = [],[],[]
    storm   = storms[s]
    
    # create 2D arrays of ERA Lat and ERA Lon points
    Lon,Lat = np.meshgrid(elon,elat)
    storm_bool = np.zeros_like(Lon)
    
    # create grid and assign storm points/amplitude to specific grid cells
    for i,stype in enumerate(storm['type']):
        
        # only filter and save cyclonic storms
        if stype == 'cyclonic':
            
            # find which lat/lon gridpoint is closest to storm
            sidx_lat = np.argmin(np.absolute(elat-storm['lat'][i]))
            sidx_lon = np.argmin(np.absolute(elon-storm['lon'][i]))
            slat,slon = elat[sidx_lat],elon[sidx_lon]
            
            # fill meshgrid cell with storm point amplitude (i.e., its pressure)
            mesh_idx = np.argwhere((Lat==slat)&(Lon==slon))[0] 
            storm_bool[mesh_idx[0],mesh_idx[1]] = storm['amp'][i]
    
    # assign NaN anwhere there isn't a storm point
    storm_bool[storm_bool==0.0] = np.nan
    
    # This while loop does all the work. 
    # First, it finds the storm point with the lowest MSLP.
    # Then, it masks out a distance around said point to remove all points associated with same storm
    # It repeats the process again until all points are masked out and there are no more minimums
    alldone = False 
    counter=0
    while alldone is False:
        # Find point with lowest MSLP
        minstorm,minp = np.argwhere(storm_bool==np.nanmin(storm_bool)),np.nanmin(storm_bool)
        
        # If "argwhere" finds a point, assign it to variable
        if len(minstorm)>0:
            minstorm=minstorm[0]
            print('Found a minimum',int(day[s]),int(hour[s]),minp)
        # if "argwhere" does not find a minimum point, then the while loop is broken
        else:
            alldone = True
            break
        
        # Save lat/lon of strongest storm point
        latmin,lonmin = elat[minstorm[0]],elon[minstorm[1]]

        # Append lat, lon, and pressure value of storm point to universal arrays
        storms_lat.append(latmin)
        storms_lon.append(lonmin)
        amp.append(minp)

        # make array of distances centered on strongest storm point
        dist_lon = np.abs(lonmin-elon) # this finds the distance (in degrees) of each gridpoint away from strong storm
        dist_lat = np.abs(latmin-elat)
        dlon,dlat = np.meshgrid(dist_lon,dist_lat) # needs to be 2D
        
        # make a new array, removing any point associated with main storm
        storm_bool[(dlat<max_lat_dist)&(dlon<max_lon_dist)]=np.nan
        counter=counter+1
        
    
    # There are times in which no points are detected...
    # ...this piece of code tells the loop to continue without doing anything below it
    if len(storms_lon)==0:
        continue
        
    # save date associated with filtered storm points
    years.append(int(year[s]))
    months.append(int(month[s]))
    days.append(int(day[s]))
    hours.append(int(hour[s]))
    
    # Save in same format as other .npz files
    storm_tmp = {}
    storm_tmp['lon'] = np.asarray(storms_lon)
    storm_tmp['lat'] = np.asarray(storms_lat)
    storm_tmp['amp'] = np.asarray(amp)
    storm_tmp['type'] = list(repeat('cyclonic',len(storms_lon)))
    storm_tmp['N'] = len(storms_lon)
    storms_file.append(storm_tmp)
    
print('saving')
np.savez('./storm_det_filtered_'+yearid+'.npz', 
         storms=storms_file,year=years,month=months,day=days,hour=hours)
