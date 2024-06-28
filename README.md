# GLStormTrends_2024

This repository contains the code necessary to create all data used in AGU GRL publication Hutson et al. 2024: _Historical Trends in Cold-Season Mid-Latitude Cyclones in the Great Lakes Region_. 

**cut_and_stitch.py**: Use this script to take two yearly ERA5 netCDF files and combine them into one cold-season ERA5 file. ;

Once all data is organized into seasonal files, the workflow to detect, filter, and track extratropical cyclones will be:
storm_detection.py &rarr storm_filtering.py &rarr storm_tracking.py
