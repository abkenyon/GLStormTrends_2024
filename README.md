# GLStormTrends_2024

This repository contains the code necessary to create all data used in AGU GRL publication Hutson et al. 2024: _Historical Trends in Cold-Season Mid-Latitude Cyclones in the Great Lakes Region_. 

 - **cut_and_stitch.py**: Use this script to take two yearly ERA5 netCDF files and combine them into one cold-season ERA5 file. ;

### Once all data is organized into seasonal files, the workflow to detect, filter, and track extratropical cyclones will be:

_storm_detection.py_ $\rightarrow$ _storm_filtering.py_ $\rightarrow$ _storm_tracking.py_

 - **storm_detection.py**: Detects potential cyclonic and anticyclonic storm centers using the mean sea level pressure (MSLP) field. 

 - **storm_filtering.py**: Filters points that have been erroneously identified as storm centers by removing any points within a 3 degree square area surrounding the strongest (i.e., lowest MSLP) point.

 - **storm_tracking.py**: Connects the cyclone center points into coherent tracks. The tracks are filtered here using the requirements described in Hutson et al. (2024). 

 - **storm_stats.py**: Can be used to loop through all storm tracks and identify/save useful statistics (e.g., track length, track duration, lowest MSLP, etc.) 

 - **storm_functions.py**: A library of functions used throughout the storm tracking workflow. Also contains a useful plotting function that creates a figure map of the Great Lakes Region, along with an axis handle and projection class allowing additional data to be plotted on the map. 

 - **storm_composite.py** Creates seasonal storm-centered composites of extratropical cyclones (ETCs) that pass through the GLR.


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
