# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

# Define where my python is located
#!/usr/bin/python

# Load the appropriate components
import os, datetime, pylab
import numpy as np
import matplotlib as mpl

# <codecell>

# Set user options
datenow = '2012082200'
fcst_hr = '000'
archive_inlocation = 1 # 1 for real-time, 2 for near real-time
datadir = '/home/xxxx/data/'

# <codecell>

# parse the  user input
yyyy = datenow[0:4]
mm = datenow[4:6]
dd = datenow[6:8]
hh = datenow[8:10]

print yyyy, mm, dd, hh

# <codecell>

# Determine the GFS file name
if archive_inlocation == 1:
     base_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.' # real-time location  
     fcst_hr_num = int(fcst_hr)
    
     if fcst_hr_num < 100:
         fcst_hr_str = '0' + str(fcst_hr_num)
     else:
         fcst_hr_str = fcst_hr
     
     fnamein = 'gfs.t' + hh + 'z.pgrbf' + fcst_hr_str + '.grib2'
     fpathin = base_url + datenow + '/' + fnamein           
else:
     base_url = 'ftp://nomads.ncdc.noaa.gov/GFS/Grid4/' # near real-time, longer archive 
     fnamein = 'gfs_4_' + yyyy + mm + dd + '_00' + hh + '_' + fcst_hr + '.grb2' 
     fpathin = base_url + yyyy + mm + '/' + yyyy + mm + dd + '/' + fnamein

print fnamein
print "File source:"
print fpathin

# <codecell>

# Get the file (only will work on SoM machines unless you have 'wget' on your computer)
# If you do not have wget, cut, copy and paste the statement last printed above this block
#   into the url box of a web browser.  This will download the file to a default location
#   on your computer, probably the 'Downloads' directory.
cmd1 = 'wget -N ' 

os.chdir(datadir) # Tell the operaring system to change directories to the data directory

os.system(cmd1+fpathin) 

# <codecell>


