#!/usr/bin/python 


# imports
import netCDF4

import os, datetime
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um


import warnings
warnings.filterwarnings("ignore")
###################################
# Set user options
###################################
date_string = '1900010100'
map_projection = 'lcc' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
plot_barbs = 'true' # 'true' to plot windbarbs, 'false' otherwise
compute_geostrophic_wind = 'true' # If true, converts true wind to geostrophic wind
record_num = 0 # time index to plot.  0 means the first time, 1 the second time...etc.
figname = "gfs_trth_analysis" # will save an image named this with a .png extension


# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/xxxx/data'
fname = 'files1234.nc';
###################################
# END user options
###################################

# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, fname)
print fpath



f = netCDF4.Dataset(fpath,'r')
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]
trtempin = f.variables['TMP_P0_L109_GLL0'][record_num,1,::-1,:].squeeze()
trpresin = f.variables['PRES_P0_L109_GLL0'][record_num,1,::-1,:].squeeze()
#slp = f.variables['PRMSL_P0_L101_GLL0'][record_num,::-1,:].squeeze()/10**2
f.close


trth = wm.temp_to_theta(trtempin, trpresin)



cbar_min_trth = 280
cbar_max_trth = 360
cint_trth = 2
cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)

base_cntr_slp = 1012
cint_slp = 2
cbar_min_slp = base_cntr_slp-20*cint_slp
cbar_max_slp = base_cntr_slp+20*cint_slp
cbar_max_slp = cbar_max_slp + (cint_slp/2)   
cflevs_slp =  np.arange(cbar_min_slp, cbar_max_slp, cint_slp)
cflevs_slp_ticks = np.arange(cbar_min_slp,cbar_max_slp,4*cint_slp)


lonin = lons
trth, lons = um.addcyclic(trth, lonin)
#slp, lons = um.addcyclic(slp, lonin)

X, Y = np.meshgrid(lons, lats)

titletext1 = 'Tropopause potential temperature valid %s at %s UTC' % (
        dt.strftime('%d %b %Y'), dt.strftime('%H00'))  

# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

# Figure 1
fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
if map_projection == 'ortho':
   m = Basemap(projection='ortho', lat_0 = 50, lon_0 = 260,
               resolution = 'l', area_thresh = 1000.,ax=ax1)
elif map_projection == 'lcc':
#   m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
    m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax1)	       
# draw countries, states, and differentiate land from water areas.
m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
m.drawcountries(linewidth=1, color='#444444', zorder=5)
m.drawstates(linewidth=0.66, color='#444444', zorder=4)
m.drawmapboundary

# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 30))

x, y = m(X, Y)

# Contour tropopause potential temperature
CS1 = m.contourf(x,y,trth,cmap=plt.cm.jet,levels=cflevs_trth, extend='both',zorder=1)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
clabs = ['%i K' % f for f in cflevs_trth_ticks]
cbar.ax.set_yticklabels(clabs, size=10) 
     
CS2 = m.contour(x, y, trth, cflevs_trth, colors='k', linewidths=0.5)
#CS3 = m.contour(x, y, slp, cflevs_slp, colors='k', linewidths=1.5)
#plt.clabel(CS3, inline=1, fontsize=10, fmt='%i')   


ax1.set_title(titletext1)
save_name = figname + "_01.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()
