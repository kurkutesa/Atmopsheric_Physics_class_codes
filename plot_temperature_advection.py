# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/python
# v0.2 9 March 2012

import netCDF4

sys.path.append("/Users/scavallo/scripts/python_scripts")
from weather_modules import *
from utilities_modules import *

import os, datetime, pylab
import numpy as np
import matplotlib as plt
from pylab import *
from mpl_toolkits.basemap import Basemap, shiftgrid

# <codecell>

fpath = "/Users/scavallo/data/gfs.t00z.pgrbf00_2012080900_f000.nc"
level_option = 100000 # Pascals 

# <codecell>

# Open the Netcdf file
f = netCDF4.Dataset(fpath,'r')

# <codecell>

# Read in select netcdf variables and close the file
lons = f.variables['lon_0'][:] 
lats = f.variables['lat_0'][::-1] # Must read latitudes in reverse direction
levs = f.variables['lv_ISBL0'][:] 
temperature = f.variables['TMP_P0_L100_GLL0'][:]
u = f.variables['UGRD_P0_L100_GLL0'][:]
v = f.variables['VGRD_P0_L100_GLL0'][:]
f.close  

# <codecell>

# Compute the gradients on a sphere
dfdp,dfdlat,dfdlon = gradient_sphere(temperature, levs, lats, lons)        

# Compute temperature advection
temperature_advection = -1*u*dfdlon + -1*v*dfdlat

# Pull out the data at the user specified level, and reverse the latitude dimension
levelindex = pylab.find(levs==level_option)
temperature_plot = temperature[levelindex,::-1,:].squeeze()
temperature_advection_plot = temperature_advection[levelindex,::-1,:].squeeze()
u_plot = u[levelindex,::-1,:].squeeze()
v_plot = v[levelindex,::-1,:].squeeze()

# Wrap the longtitudes around so they are cyclic
lonin = lons
temperature_advection_plot, lons = addcyclic(temperature_advection_plot, lonin)
temperature_plot, lons = addcyclic(temperature_plot, lonin)
u_plot, lons = addcyclic(u_plot, lonin)
v_plot, lons = addcyclic(v_plot, lonin)

# Refresh the dimensions
[X,Y] = np.meshgrid(lons,lats)
[ny,nx] = np.shape(X)

levelh = level_option / 100 # Convert level to hPa for title

# <codecell>

# Set global figure properties
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden),dpi=128)
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top = 0.93, wspace=0.2, hspace=0.2)

# <codecell>

# Plot temperature and winds
clevs = np.arange(-40,41,5) # Contours from -40 to 40 with a contour interval of 5

temperature_plot = temperature_plot - 273.15 # Plot in degrees Celsius

fig = plt.figure(**figprops)   # New figure   
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
    
map = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax1)	       
map.drawcoastlines(linewidth=2, color='#444444', zorder=6)
map.drawcountries(linewidth=1, color='#444444', zorder=5)
map.drawstates(linewidth=0.66, color='#444444', zorder=4)
map.drawmapboundary

x,y = map(X,Y)
cf = map.contour(x,y,temperature_plot,levels=clevs, extend='both',zorder=1,linewidths=2)
plt.clabel(cf, inline=1, fontsize=10)

# plot wind vectors on projection grid.
# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
ugrid,newlons = shiftgrid(180.,u_plot,lons,start=False)
vgrid,newlons = shiftgrid(180.,v_plot,lons,start=False)
# transform vectors to projection grid and plot windbarbs.
uproj,vproj, xx, yy = map.transform_vector(ugrid,vgrid,newlons,lats,41,41,returnxy=True)
barbs = map.barbs(xx,yy,uproj,vproj,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)
# Title
ax1.set_title(str(levelh)+' hPa temperature and wind barbs')

# <codecell>

# Temperature advection plot
cbar_min = -0.0001
cbar_max = 0.0001   
cfint = 0.00001    
cflevs = np.arange(cbar_min,cbar_max+cfint-(cfint/2),cfint)
ncolors = np.size(cflevs)-1
cmap = cm.get_cmap('jet')
djet = cmap_whitezero(cmap,ncolors,1,-1)    


plt.clf()
fig = plt.figure(**figprops)
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
map = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax = ax1)	       
map.drawcoastlines(linewidth=2, color='#444444', zorder=6)
map.drawcountries(linewidth=1, color='#444444', zorder=5)
map.drawstates(linewidth=0.66, color='#444444', zorder=4)
map.drawmapboundary

x,y = map(X,Y)
cf = map.contourf(x,y,temperature_advection_plot,cmap=djet,levels=cflevs, extend='both',zorder=1)
cbar = plt.colorbar(cf, shrink=0.95, orientation='horizontal',extend='both')



plt.show()

# <codecell>


