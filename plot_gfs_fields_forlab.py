# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/python

# imports
import netCDF4

import os, datetime
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')              # Switch to handle no X11 display.

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um


# <codecell>

###################################
# Set user options
###################################
date_string = 'yyyymmddhh' # yyyymmddhh
level_option = 50000 # Pascals; set to -1 for sea level pressure
map_projection = 'lcc' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
plot_barbs = 'true' # 'true' to plot windbarbs, 'false' otherwise
plot_contours = 'true' # True to plot contours, false otherwise
figname = "gfs_500hPa_analysis" # will save an image named this with a .png extension

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/Users/xxxx/data'

###################################
# END user options
###################################

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
        
   temperature_plot = f.variables['TMP_P0_L100_GLL0'][levelindex,::-1,:].squeeze()
   u_plot =  f.variables['UGRD_P0_L100_GLL0'][levelindex,::-1,:].squeeze()
   v_plot =  f.variables['VGRD_P0_L100_GLL0'][levelindex,::-1,:].squeeze() 

f.close

# Convert temperature to degrees Celsius
temperature_plot = temperature_plot - 273.15

# Convert wind to knots
u_plot = u_plot * 1.94384
v_plot = v_plot * 1.94384

lonin = lons
plotvar, lons = um.addcyclic(plotvar, lonin)
temperature_plot, lons = um.addcyclic(temperature_plot, lonin)
u_plot, lons = um.addcyclic(u_plot, lonin)
v_plot, lons = um.addcyclic(v_plot, lonin)

# Refresh the dimensions
X, Y = np.meshgrid(lons, lats)

levelh = level_option / 100 # Convert level to hPa for title

# <codecell>

# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

# <codecell>
# Setting contour interval done here
if (level_option == -1):
    base_cntr = 1012
    cint = 4    
elif (level_option == 50000):
    base_cntr = 5580 # a contour close to standard atmospheric value
    cint = 60 # contour interval
elif (level_option == 30000):
    base_cntr = 9180 # a contour close to standard atmospheric value
    cint = 60 # contour interval
    
    cflevs_wind = np.arange(0, 101, 10)
    cflevs_shade = np.arange(70,101,10)
else:
    base_cntr = 5580 # Base contour on surfaces other than 500 hPa or 300 hPa.  Change this for your 850 hPa plot    
    cint = 60 # Contour interval on surfaces other than 500 hPa or 300 hPa.  

cbar_min = base_cntr-20*cint
cbar_max = base_cntr+20*cint
cflevs = np.arange(cbar_min, cbar_max+1, cint)

# Temperature contour interval from -40 degrees Celsius to positive 40 degrees Celsius, with a contour interval of 2 degrees Celsius
cflevs_temp = np.arange(-40,41,2)


if (level_option == -1 ):
   titletext = 'Sea level pressure valid %s at %s UTC' % (
                dt.strftime('%d %b %Y'), dt.strftime('%H00'))
else:
   titletext = '%s hPa geopotential heights valid %s at %s UTC' % (
                levelh, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
print titletext

fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])

if map_projection == 'ortho':
   m = Basemap(projection='ortho', lat_0 = 50, lon_0 = 260,
               resolution = 'l', area_thresh = 1000.,ax=ax1)
elif map_projection == 'lcc':
   m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax1)	       
# draw countries, states, and differentiate land from water areas.
m.drawcountries()
m.drawstates()
m.drawlsmask(land_color='0.7', ocean_color='white', lakes=True)

# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 30))


x, y = m(X, Y)

# Set up to plot wind vectors on projection grid.
ugrid,newlons = shiftgrid(180.,u_plot,lons,start=False) 
vgrid,newlons = shiftgrid(180.,v_plot,lons,start=False)
# Rotate and interpolate wind from lat/lon grid to map projection grid.
uproj,vproj, xx, yy = m.transform_vector(ugrid,vgrid,newlons,lats,41,41,returnxy=True)


if ( plot_contours == 'true' ):
   CS = m.contour(x, y, plotvar, cflevs, colors='k', linewidths=1.5)
   plt.clabel(CS, inline=1, fontsize=10, fmt='%i')
   

if plot_barbs == "true":
   barbs = m.barbs(xx,yy,uproj,vproj,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)

if (level_option == 30000):
   
   windmag = np.sqrt(u_plot**2 + v_plot**2)
   
   if ( plot_contours == 'true' ):
      CS2 = m.contour(x, y, windmag, cflevs_wind, colors='b', linewidths=1.0)
      plt.clabel(CS2, inline=1, fontsize=10, fmt='%i')      
   
      CS3 = m.contourf(x,y,windmag,cflevs_shade,cmap=plt.cm.cool)
elif (level_option == 50000):
    if ( plot_contours == 'true' ):  
        CS2 = m.contour(x, y, temperature_plot, cflevs_temp, colors='r', linestyles='dashed',linewidths=1.5)
        plt.clabel(CS2, inline=1, fontsize=10, fmt='%i')      


ax1.set_title(titletext)
save_name = figname + ".png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()

