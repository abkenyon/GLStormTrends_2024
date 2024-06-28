'''

    Software for the tracking of storms
    based on detected and filtered storm position data.
  
    Originally created by Eric J. Oliver for the stormTracking algorithm
    (https://github.com/ecjoliver/stormTracking)
  
    Modified by Abby Hutson and Ayumi Fujisaki-Manome for use with ERA5 
    and for storms in the Great Lakes Region.

'''

# Load required modules

import numpy as np
import storm_functions as storm
import sys

#
# Automated storm tracking
#

yearid=sys.argv[1]

# Load in detected positions and date/hour information
filename = './storm_det_filtered_'+yearid # NOTE: Be sure to read in filtered data (output from storm_filtering.py)
data = data=np.load(filename+'.npz',allow_pickle=True,encoding='latin1')
det_storms = data['storms']
year = data['year']
month = data['month']
day = data['day']
hour = data['hour']

# Initialize storms discovered at first time step

storms = storm.storms_init(det_storms, year, month, day, hour)

# Stitch storm tracks together at future time steps

T = len(det_storms) # number of time steps
for tt in range(1, T-1):
    # Track storms from time step tt-1 to tt and update corresponding tracks and/or create new storms
    storms = storm.track_storms(storms, det_storms, tt, year, month, day, hour, dt=6)
    #print det_storms

# Add keys for storm age and flag if storm was still in existence at end of run
for ed in range(len(storms)):
    storms[ed]['age'] = len(storms[ed]['lon'])

# Strip storms based on track lengths
storms = storm.strip_storms(storms, dt=6, d_tot_min=300, d_ratio=0.6, dur_min=36)


# Save tracked storm data
np.savez('storm_track_slp_'+yearid, storms=storms)

