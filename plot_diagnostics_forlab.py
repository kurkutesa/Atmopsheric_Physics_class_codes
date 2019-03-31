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


###################################
# Set user options
###################################
date_string = '2012091812' # yyyymmddhh
map_projection = 'lcc' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
plot_barbs = 'true' # 'true' to plot windbarbs, 'false' otherwise
compute_geostrophic_wind = 'true' # If true, converts true wind to geostrophic wind
levels_x1 = [10000,20000] # Set the level(s) you want to plot vorticity at
levels_x2 = [15000] 
figname = "gfs_qg_diagnostics" # will save an image named this with a .png extension



# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/xxx/data'

###################################
# END user options
###################################

# <codecell>

# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, 'gfs_4_%s_%s00_000.nc' % (
                     dt.strftime('%Y%m%d'), dt.strftime('%H')))
print fpath

# Open netcdf file 
f = netCDF4.Dataset(fpath,'r')
# Read in latitude, longitude, and pressure vectors
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]

# Find the indices where the pressure vector equals the desired pressure levels for your diagnostics
levelindex1 = np.ravel(levs==levels_x1[0])
levelindex2 = np.ravel(levs==levels_x1[1])
levelindex3 = np.ravel(levs==levels_x2[0])

# u-component of wind corresponding to first entry of levels_x1
uin1 =  f.variables['UGRD_P0_L100_GLL0'][levelindex1,::-1,:].squeeze()
# u-component of wind corresponding to second entry of levels_x1
uin2 =  f.variables['UGRD_P0_L100_GLL0'][levelindex2,::-1,:].squeeze()

# Same as for u, but for v component of wind
vin1 =  f.variables['VGRD_P0_L100_GLL0'][levelindex1,::-1,:].squeeze()
vin2 =  f.variables['VGRD_P0_L100_GLL0'][levelindex2,::-1,:].squeeze()

# Same as for wind, but for geopotential heights
ghgtin1 = f.variables['HGT_P0_L100_GLL0'][levelindex1,::-1,:].squeeze()
ghgtin2 = f.variables['HGT_P0_L100_GLL0'][levelindex2,::-1,:].squeeze()

# Omega
omegain = f.variables['VVEL_P0_L100_GLL0'][levelindex3,::-1,:].squeeze()

# Close netcdf file
f.close

# Convert wind to geostrophic wind
if compute_geostrophic_wind == 'true' :
   uin1,vin1 = wm.geostrophic_latlon(uin1, vin1, ghgtin1, lats, lons)
   uin2,vin2 = wm.geostrophic_latlon(uin2, vin2, ghgtin2, lats, lons)

# Compute vertical voriticity
zeta_a1 = wm.vertical_vorticity_latlon(uin1, vin1, lats, lons, 1)
zeta_a2 = wm.vertical_vorticity_latlon(uin2, vin2, lats, lons, 1)

# Make all fields to be plotted cyclic
lonin = lons
u_plot1, lons = um.addcyclic(uin1, lonin)
v_plot1, lons = um.addcyclic(vin1, lonin)
u_plot2, lons = um.addcyclic(uin2, lonin)
v_plot2, lons = um.addcyclic(vin2, lonin)
zeta_a1, lons = um.addcyclic(zeta_a1, lonin)
zeta_a2, lons = um.addcyclic(zeta_a2, lonin)
omegain, lons = um.addcyclic(omegain, lonin)

# For the title
levelh1 = levels_x1[0] / 100 # Convert level to hPa for title
levelh2 = levels_x1[1] / 100 # Convert level to hPa for title

# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

# Colorbar settings
# Setting contour interval done here
base_cntr = 0
cint = 1*10**-5   
cbar_min = 8*cint
cbar_max = base_cntr+26*cint
cbar_max = cbar_max + (cint/2)
cflevs = np.arange(cbar_min, cbar_max, cint)
cflevs_ticks = np.arange(cbar_min,cbar_max,4*cint)


base_cntr2 = 0
cint2 = 2*10**-10   
cbar_min2 = base_cntr2-14*cint2
cbar_max2 = base_cntr2+14*cint2
cbar_max2 = cbar_max2 + (cint2/2)
cflevs2 = np.arange(cbar_min2, cbar_max2, cint2)
cflevs2_ticks = np.arange(cbar_min,cbar_max,4*cint)



# Set titletexts here
titletext1 = '%s hPa geostrophic wind and vorticity valid %s at %s UTC' % (
             levelh1, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
titletext2 = '%s hPa geostrophic wind and vorticity valid %s at %s UTC' % (
             levelh2, dt.strftime('%d %b %Y'), dt.strftime('%H00'))


# Refresh the dimensions
X, Y = np.meshgrid(lons, lats)	     
	     
# Set up to plot wind vectors on projection grid.
ugrid1,newlons = shiftgrid(180.,u_plot1,lons,start=False) 
vgrid1,newlons = shiftgrid(180.,v_plot1,lons,start=False)
ugrid2,newlons = shiftgrid(180.,u_plot2,lons,start=False) 
vgrid2,newlons = shiftgrid(180.,v_plot2,lons,start=False)


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
CS = m.contourf(x,y,zeta_a1,cmap=plt.cm.hot_r,levels=cflevs, extend='both',zorder=1)
cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
CS2 = m.contour(x, y, zeta_a1, cflevs, colors='k', linewidths=0.5)

if plot_barbs == "true":
    # Rotate and interpolate wind from lat/lon grid to map projection grid.
   uproj1,vproj1, xx1, yy1 = m.transform_vector(ugrid1,vgrid1,newlons,lats,41,41,returnxy=True)    
   barbs = m.barbs(xx1,yy1,uproj1,vproj1,length=5,barbcolor='k',flagcolor='r',linewidth=1.0)


cbar.set_ticks(cflevs_ticks)
clabs = ['%i K' % f for f in cflevs_ticks]
cbar.ax.set_yticklabels(clabs, size=10)
ax1.set_title(titletext1)
save_name = figname + "_01.png"
plt.savefig(save_name, bbox_inches='tight')


# Figure 2
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
CS = m.contourf(x,y,zeta_a2,cmap=plt.cm.hot_r,levels=cflevs, extend='both',zorder=1)
cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
CS2 = m.contour(x, y, zeta_a2, cflevs, colors='k', linewidths=0.5)

if plot_barbs == "true":
    # Rotate and interpolate wind from lat/lon grid to map projection grid.
   uproj2,vproj2, xx2, yy2 = m.transform_vector(ugrid2,vgrid2,newlons,lats,41,41,returnxy=True)    
   barbs = m.barbs(xx2,yy2,uproj2,vproj2,length=5,barbcolor='k',flagcolor='r',linewidth=1.0)
cbar.set_ticks(cflevs2_ticks)
clabs = ['%i K' % f for f in cflevs2_ticks]
cbar.ax.set_yticklabels(clabs, size=10)
ax1.set_title(titletext2)
save_name = figname + "_02.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()

	     
	     

