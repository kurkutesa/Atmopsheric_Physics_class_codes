# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/python

# imports
import netCDF4

import os, datetime
import numpy as np
import matplotlib as mpl
mpl.use('Agg')              # Switch to handle no X11 display.
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Add a couple of user defined functions
# from weather_modules import *
# from utilities_modules import *
import weather_modules as wm
import utilities_modules as um

# <codecell>

# Set user options
date_string = '2012082312' # yyyymmddhh
level_option = 50000 # Pascals; set to -1 for sea level pressure

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/pmarsh/metr4424/lab01/data'




# <codecell>

# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, 'gfs_4_%s_%s00_000.nc' % (
                     dt.strftime('%Y%m%d'), dt.strftime('%H')))
print fpath
f = netCDF4.Dataset(fpath,'r')
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]
if ( level_option == -1 ) :
   plotvar  = f.variables['PRMSL_P0_L101_GLL0'][::-1,:]/100
else:
    levelindex = np.ravel(levs==level_option)
    plotvar = f.variables['HGT_P0_L100_GLL0'][levelindex,::-1,:].squeeze() # Reverse latitude dimension
f.close

# <codecell>

lonin = lons
plotvar, lons = um.addcyclic(plotvar, lonin)

# Refresh the dimensions
X, Y = np.meshgrid(lons, lats)

levelh = level_option / 100 # Convert level to hPa for title

if (level_option == -1 ):
   titletext = 'Sea level pressure valid %s at %s UTC' % (
                dt.strftime('%d %b %Y'), dt.strftime('%H00'))
else:
   titletext = '%s hPa geopotential heights valid %s at %s UTC' % (
                levelh, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
print titletext

# <codecell>

# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

# <codecell>
# Setting contour interval done here
if (level_option == -1):
    cint = 4
    cbar_min = 1012-20*cint
    cbar_max = 1012+20*cint
else:
    base_cntr = 5580 # a contour close to standard atmospheric value
    cint = 60 # contour interval
    cbar_min = base_cntr-20*cint
    cbar_max = base_cntr+20*cint

cflevs = np.arange(cbar_min, cbar_max+1, cint)

fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])

m = Basemap(projection='ortho', lat_0 = 50, lon_0 = 260,
               resolution = 'l', area_thresh = 1000.,ax=ax1)
# draw countries, states, and differentiate land from water areas.
m.drawcountries()
m.drawstates()
m.drawlsmask(land_color='0.7', ocean_color='white', lakes=True)

# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 30))


x, y = m(X, Y)
CS = m.contour(x, y, plotvar, cflevs, colors='k', linewidths=1.5)
plt.clabel(CS, inline=1, fontsize=10, fmt='%i')
ax1.set_title(titletext)

figname = "example"
save_name = figname + ".png"
plt.savefig(save_name, bbox_inches='tight')
# plt.show()

