#########################################
# cut_and_stitch.py is designed to create seasonal ERA5 netCDF files
# by "cutting and stitching" together two different ERA5 netCDF files
# from two different years. The result is a data file with dates 
# October - March
#
# Last edited by Abby Hutson (2024-06-28)
############################################

#import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
import os
import sys

# Ryan's ERA-5 data directory
dir="/nfs/turbo/seas-ayumif/hutsona/data_big_domain/"


# year to process
year=int(sys.argv[1])

# time interval of ERA-5 data (6 hours)
dt = timedelta(hours=6)

# load the NCO library
os.system("module load nco")


print("\nProcessing data for "+str(year-1)+"-"+str(year)+"....\n")

# Get time indices where data files need to be split. 
# Open the fiest file for "year" (Jan-Mar), read the time info.
fn=dir+"era5_"+str(year)+".nc"
ncfile = Dataset(fn, 'r')
time=ncfile.variables['time'][:]  # days since 1970-01-01 00:00
timeunit=ncfile.variables['time'].units 
ncfile.close()
# convert to datetime object
date=num2date(time,units=timeunit,calendar='standard')
# loop in time and find an index where a time gap (>6hr) happens
for nn in range(1,len(time)):
    # calculate time interval 
    time_difference = date[nn]-date[nn-1]
    if time_difference > dt: 
        winter_end_index = nn-1
        print("In "+dir+"era5_"+str(year)+".nc, winter_end_index is",winter_end_index,"and the end time is",str(date[winter_end_index]),"\n")

# Use a NCO command to extract data for Jan-Mar in "year"
os.system("ncks -O -h -d time,0,"+str(winter_end_index)+" "+fn+" ./tmp_jan-mar.nc")
os.system("ncks -O -h --mk_rec_dmn time ./tmp_jan-mar.nc ./tmp_jan-mar.nc")

# Open the second file for "year"-1 (Oct-Dec)
fn=dir+"era5_"+str(year-1)+"_snowdepth.nc"
ncfile = Dataset(fn, 'r')
time=ncfile.variables['time'][:]  # days since 1970-01-01 00:00
timeunit=ncfile.variables['time'].units
ncfile.close()
# convert to datetime object
date=num2date(time,units=timeunit,calendar='standard')
# loop in time abd find an index where a time gap >6hr happends
for nn in range(1,len(time)):
    # calculate time interval 
    time_difference = date[nn]-date[nn-1]
    if time_difference > dt:
        winter_begin_index = nn
        print("In "+dir+"era5_"+str(year-1)+".nc, winter_begin_index is",winter_begin_index,"and the begin time is",str(date[winter_begin_index]),"\n")

# Use a NCO command to extract data for Oct-Dec in "year"-1
os.system("ncks -O -h -d time,"+str(winter_begin_index)+",-1  "+fn+" ./tmp_oct-dec.nc")
os.system("ncks -O -h --mk_rec_dmn time ./tmp_oct-dec.nc ./tmp_oct-dec.nc")



# Use a NCO command to stitch the two files
os.system("ncrcat -O -h tmp_oct-dec.nc tmp_jan-mar.nc era5_"+str(year-1)+"-"+str(year)+".nc")

# Use a NCO command to remove temperary files
os.system("rm tmp_oct-dec.nc tmp_jan-mar.nc")

# Check time info
ncfile = Dataset("./"+"era5_"+str(year-1)+"-"+str(year)+".nc", 'r')
time=ncfile.variables['time'][:]  # days since 1970-01-01 00:00
timeunit=ncfile.variables['time'].units
# convert to datetime object
date=num2date(time,units=timeunit,calendar='standard')
print("The merged file era5_"+str(year-1)+"-"+str(year)+".nc has ",len(time),
" time records that range from "+str(date[0])+" to "+str(date[len(time)-1]),"\n")
print("Completed.\n")

# Use below if you want to check time in the merged file
# loop in time
#for nn in range(1,len(time)):
#    # calculate time interval 
#    time_difference = date[nn]-date[nn-1]
#    print(nn,str(date[nn]),time_difference)
#    if time_difference > dt:
#       print("invalid time interval")
#       exit






