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
date_string = '2012091812'
fname2 = 'gfs_4_20120918_1200_ghgtonly.nc'
map_projection = 'lcc' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
plot_barbs = 'true' # 'true' to plot windbarbs, 'false' otherwise
compute_geostrophic_wind = 'true' # If true, converts true wind to geostrophic wind
diagnostic = 1 # 1 for height tendency equation, 2 for traditional omega
levels_upperlower = [70000,30000] # Set the level(s) you want to plot vorticity at
levels_mid = [50000] # Set the level(s) you want to plot temperature advection at
figname = "gfs_500hPa_analysis" # will save an image named this with a .png extension



# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/xxxx/data'

###################################
# END user options
###################################
if diagnostic == 1:
   fpath2 = os.path.join(fpath,fname2)
   print fpath2   


# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, 'gfs_4_%s_%s00_000.nc' % (
                     dt.strftime('%Y%m%d'), dt.strftime('%H')))
print fpath



f = netCDF4.Dataset(fpath,'r')
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]


levelindex_lower = np.ravel(levs==levels_upperlower[0])
levelindex_upper = np.ravel(levs==levels_upperlower[1])
levelindex_mid = np.ravel(levs==levels_mid[0])


uin_lower =  f.variables['UGRD_P0_L100_GLL0'][levelindex_lower,::-1,:].squeeze()
uin_upper =  f.variables['UGRD_P0_L100_GLL0'][levelindex_upper,::-1,:].squeeze()

vin_lower =  f.variables['VGRD_P0_L100_GLL0'][levelindex_lower,::-1,:].squeeze()
vin_upper =  f.variables['VGRD_P0_L100_GLL0'][levelindex_upper,::-1,:].squeeze()

uin_mid = f.variables['UGRD_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()
vin_mid = f.variables['UGRD_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()
if diagnostic == 1:
   # Read in temperature at the lower level and call it temp_lower
   # Also read in temperature at the upper level and call it temp_upper
   putstuffhere = 1

elif diagnostic == 2:
   temp_mid = f.variables['TMP_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()
else:
   temp_mid = f.variables['TMP_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()  
   
ghgtin_lower = f.variables['HGT_P0_L100_GLL0'][levelindex_lower,::-1,:].squeeze()
ghgtin_upper = f.variables['HGT_P0_L100_GLL0'][levelindex_upper,::-1,:].squeeze()
ghgtin_mid = f.variables['HGT_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()

omegain = f.variables['VVEL_P0_L100_GLL0'][levelindex_mid,::-1,:].squeeze()

f.close

if diagnostic == 1:
   # Reading geopotential height from the second file.  Index 0 is the analysis field, and index 1 is the 3 hour forecast.  We already have the analysis field from above, so no need to read that in again.   
   f = netCDF4.Dataset(fpath2,'r')
   ghgtin_mid2 = f.variables['HGT_P0_L100_GLL0'][1,levelindex_mid,::-1,:].squeeze()
   f.close
   geoph_tend = ghgtin_mid2 - ghgtin_mid



# Compute the geostrophic wind
if compute_geostrophic_wind == 'true' :
   uin_lower,vin_lower = wm.geostrophic_latlon(uin_lower, vin_lower, ghgtin_lower, lats, lons)
   uin_upper,vin_upper = wm.geostrophic_latlon(uin_upper, vin_upper, ghgtin_upper, lats, lons)
   uin_mid,vin_mid = wm.geostrophic_latlon(uin_mid, vin_mid, ghgtin_mid, lats, lons)



lonin = lons
if diagnostic == 1: # height tendency
   # Follow the example for diagnostic==2 below.  The key is to make sure your variable names end up as:
   # var_mid, var_lower, var_upper, var_adv_lower, var_adv_upper, and var_adv_mid
   putstuffhere = 1

   # compute absolute geostrophic vorticity at the middle level (call it var_mid)
   
   # var_lower is temp_lower and var_upper is temp_upper
   var_lower = temp_lower
   var_upper = temp_upper
   
   # Calculate temperature advection at lower and upper level for differential thermal advection term later

   
   # Calculate absolute vorticity advection at middle level     

   # Need to wrap the tendency array
   geoph_tend, lons = um.addcyclic(geoph_tend, lonin)
   
   # Convert temperature to degrees Celsius
   var_lower = var_lower - 273.15
   var_upper = var_upper - 273.15
elif diagnostic == 2: # traditional omega
   # Temperature at 700 hPa
   var_mid = temp_mid   

   # Vertical relative and absolute voriticity
   var_lower = wm.vertical_vorticity_latlon(uin_lower, vin_lower, lats, lons, 1)
   var_upper = wm.vertical_vorticity_latlon(uin_upper, vin_upper, lats, lons, 1)

   # for differential vorticity advection 
   var_adv_lower = wm.hadvection_latlon(uin_lower, vin_lower, var_lower, lats, lons)
   var_adv_upper = wm.hadvection_latlon(uin_upper, vin_upper, var_upper, lats, lons)
   
   # temperature advection level middle level
   var_adv_mid = wm.hadvection_latlon(uin_mid, vin_mid, temp_mid, lats, lons)
   
   # convert temperature to celsius
   var_mid = var_mid - 273.15
  
else:
   var_lower = wm.vertical_vorticity_latlon(uin_lower, vin_lower, lats, lons, 1)   


var_lower, lons = um.addcyclic(var_lower, lonin)
var_upper, lons = um.addcyclic(var_upper, lonin) 
var_mid, lons = um.addcyclic(var_mid, lonin) 

u_plot_lower, lons = um.addcyclic(uin_lower, lonin)
u_plot_upper, lons = um.addcyclic(uin_upper, lonin)
u_plot_mid, lons = um.addcyclic(uin_mid, lonin)
v_plot_lower, lons = um.addcyclic(vin_lower, lonin)
v_plot_upper, lons = um.addcyclic(vin_upper, lonin)
v_plot_mid, lons = um.addcyclic(vin_mid, lonin)

ghgt_plot_lower, lons = um.addcyclic(ghgtin_lower, lonin)
ghgt_plot_upper, lons = um.addcyclic(ghgtin_upper, lonin)
ghgt_plot_mid, lons = um.addcyclic(ghgtin_mid, lonin)


var_adv_lower, lons = um.addcyclic(var_adv_lower, lonin)
var_adv_upper, lons = um.addcyclic(var_adv_upper, lonin)
var_adv_mid, lons = um.addcyclic(var_adv_mid, lonin)
omegain, lons = um.addcyclic(omegain, lonin)


diff_varadv = var_adv_upper - var_adv_lower


# Convert wind to knots
uin_lower = uin_lower * 1.94384
vin_lower = vin_lower * 1.94384
uvin_lower = np.sqrt(uin_lower**2 + vin_lower**2)

uin_upper = uin_upper * 1.94384
vin_upper = vin_upper * 1.94384
uvin_upper = np.sqrt(uin_upper**2 + vin_upper**2)

uin_mid = uin_mid * 1.94384
vin_mid = vin_mid * 1.94384
uvin_mid = np.sqrt(uin_mid**2 + vin_mid**2)


# Refresh the dimensions
X, Y = np.meshgrid(lons, lats)

level_hlower = levels_upperlower[0] / 100 # Convert level to hPa for title
level_hupper = levels_upperlower[1] / 100
level_hmid = levels_mid[0] / 100


# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

# Setting contour interval done here
#base_cntr_upper = 0
base_cntr_diffadv = 0
base_cntr_midadv = 0
base_cntr_omega = 0


if diagnostic == 1:
   base_cntr_mid = 0
   base_cntr_geophtend = 0
   base_cntr_upper = -40
   base_cntr_lower = 0

   cint_upper = 2          
   cint_lower = 2 
   cint_mid = 8*10**-6     
   cint_diffadv = 4*10**-6   
   cint_midadv = 5*10**-11   
   cint_geophtend = 2
     

   cbar_min_upper = base_cntr_upper-40*cint_upper
   cbar_max_upper = base_cntr_upper+40*cint_upper
   cbar_max_upper = cbar_max_upper + (cint_upper/2)
   
   cbar_min_lower = base_cntr_lower-40*cint_lower
   cbar_max_lower = base_cntr_lower+40*cint_lower
   cbar_max_lower = cbar_max_lower + (cint_lower/2)
   
   cbar_min_mid = 8*cint_mid
   cbar_max_mid = cbar_min_mid+14*cint_mid
   cbar_max_mid = cbar_max_mid + (cint_mid/2)

   cbar_min_diffadv = base_cntr_diffadv-14*cint_diffadv
   cbar_max_diffadv = base_cntr_diffadv+14*cint_diffadv
   cbar_max_diffadv = cbar_max_diffadv + (cint_diffadv/2)   

   cbar_min_midadv = base_cntr_midadv-14*cint_midadv
   cbar_max_midadv = base_cntr_midadv+14*cint_midadv
   cbar_max_midadv = cbar_max_midadv + (cint_midadv/2)   
   
   cbar_min_geophtend = base_cntr_geophtend-15*cint_geophtend
   cbar_max_geophtend = base_cntr_geophtend+15*cint_geophtend
   cbar_max_geophtend = cbar_max_geophtend + (cint_geophtend/2)      

   # Set titletexts here

   titletext1 = '%s hPa temperature and wind valid %s at %s UTC' % (
        	level_hupper, dt.strftime('%d %b %Y'), dt.strftime('%H00'))   
   titletext2 = '%s hPa abs. geostrophic vorticity and wind %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))

   titletext3 = '%s hPa - %s hPa diff. temperature adv. valid %s at %s UTC' % (
        	level_hupper, level_hlower, dt.strftime('%d %b %Y'), dt.strftime('%H00'))

   titletext4 = '%s hPa abs. geostrophic vort. adv. valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
		
   titletext5 = '%s hPa height tendency valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))   		

elif diagnostic == 2:
   base_cntr_mid = -10
   base_cntr_upper = 0
   base_cntr_lower = 0

   cint_upper = 8*10**-6      
   cint_lower = 8*10**-6         
   cint_mid = 4
   cint_diffadv = 6*10**-11   
   cint_midadv = 4*10**-6   
   cint_omega = 0.2*10**-1
   
   
   cbar_min_upper = 8*cint_upper
   cbar_max_upper = base_cntr_upper+26*cint_upper
   cbar_max_upper = cbar_max_upper + (cint_upper/2)

   cbar_min_lower = 8*cint_lower
   cbar_max_lower = base_cntr_lower+26*cint_lower
   cbar_max_lower = cbar_max_lower + (cint_lower/2)   

   cbar_min_mid = base_cntr_mid-14*cint_mid
   cbar_max_mid = base_cntr_mid+14*cint_mid
   cbar_max_mid = cbar_max_mid + (cint_mid/2)   

   cbar_min_diffadv = base_cntr_diffadv-14*cint_diffadv
   cbar_max_diffadv = base_cntr_diffadv+14*cint_diffadv
   cbar_max_diffadv = cbar_max_diffadv + (cint_diffadv/2)   

   cbar_min_midadv = base_cntr_midadv-14*cint_midadv
   cbar_max_midadv = base_cntr_midadv+14*cint_midadv
   cbar_max_midadv = cbar_max_midadv + (cint_midadv/2)   

   
   cbar_min_omega = base_cntr_omega-14*cint_omega
   cbar_max_omega = base_cntr_omega+14*cint_omega
   cbar_max_omega = cbar_max_omega + (cint_omega/2)

   # Set titletexts here
   titletext1 = '%s hPa geostrophic wind and vorticity valid %s at %s UTC' % (
        	level_hupper, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
   titletext2 = '%s hPa temperature and wind valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
   titletext3 = '%s hPa - %s hPa diff. geo. vort. adv. valid %s at %s UTC' % (
        	level_hupper, level_hlower, dt.strftime('%d %b %Y'), dt.strftime('%H00'))

   titletext4 = '%s hPa temperature advection valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))


   titletext5 = '%s hPa omega valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))   
else:
   cint_upper = 8*10**-6      
   cint_diffadv = 6*10**-11   
   cint_midadv = 4*10**-6   
   cint_omega = 0.2*10**-1
   
   cbar_min_upper = 8*cint_upper
   cbar_max_upper = base_cntr_upper+26*cint_upper
   cbar_max_upper = cbar_max_upper + (cint_upper/2)

   cbar_min_diffadv = base_cntr_diffadv-14*cint_diffadv
   cbar_max_diffadv = base_cntr_diffadv+14*cint_diffadv
   cbar_max_diffadv = cbar_max_diffadv + (cint_diffadv/2)   

   cbar_min_midadv = base_cntr_midadv-14*cint_midadv
   cbar_max_midadv = base_cntr_midadv+14*cint_midadv
   cbar_max_midadv = cbar_max_midadv + (cint_midadv/2)   
   
   cbar_min_omega = base_cntr_omega-14*cint_omega
   cbar_max_omega = base_cntr_omega+14*cint_omega
   cbar_max_omega = cbar_max_omega + (cint_omega/2)

   # Set titletexts here
   titletext1 = '%s hPa geostrophic wind and vorticity valid %s at %s UTC' % (
        	level_hupper, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
   titletext2 = '%s hPa - %s hPa diff. geo. vort. adv. valid %s at %s UTC' % (
        	level_hupper, level_hlower, dt.strftime('%d %b %Y'), dt.strftime('%H00'))

   titletext3 = '%s hPa temperature advection valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))
   titletext4 = '%s hPa omega valid %s at %s UTC' % (
        	level_hmid, dt.strftime('%d %b %Y'), dt.strftime('%H00'))



cflevs_upper = np.arange(cbar_min_upper, cbar_max_upper, cint_upper)
cflevs_lower = np.arange(cbar_min_lower, cbar_max_lower, cint_lower)
cflevs_mid = np.arange(cbar_min_mid, cbar_max_mid, cint_mid)
cflevs_diffadv = np.arange(cbar_min_diffadv, cbar_max_diffadv, cint_diffadv)
cflevs_midadv = np.arange(cbar_min_midadv, cbar_max_midadv, cint_midadv)



cflevs_upper_ticks = np.arange(cbar_min_upper,cbar_max_upper,4*cint_upper)
cflevs_mid_ticks = np.arange(cbar_min_mid,cbar_max_mid,4*cint_mid)
cflevs_diffadv_ticks = np.arange(cbar_min_diffadv,cbar_max_diffadv,4*cint_diffadv)
cflevs_midadv_ticks = np.arange(cbar_min_midadv,cbar_max_midadv,4*cint_midadv)

if diagnostic == 1:
   cflevs_geophtend = np.arange(cbar_min_geophtend, cbar_max_geophtend, cint_geophtend)
   cflevs_geophtend_ticks = np.arange(cbar_min_geophtend,cbar_max_geophtend,4*cint_geophtend)

if diagnostic == 2:
   cflevs_omega = np.arange(cbar_min_omega, cbar_max_omega, cint_omega)
   cflevs_omega_ticks = np.arange(cbar_min_omega,cbar_max_omega,4*cint_omega)

base_cntr = 9180 # a contour close to standard atmospheric value
cint = 60 # contour interval
cbar_min = base_cntr-20*cint
cbar_max = base_cntr+20*cint
cflevs_cntr = np.arange(cbar_min, cbar_max+1, cint)

base_cntr = 5400 # a contour close to standard atmospheric value
cint = 60 # contour interval
cbar_min = base_cntr-20*cint
cbar_max = base_cntr+20*cint
cflevs_ghgts = np.arange(cbar_min, cbar_max+1, cint)




# Set up to plot wind vectors on projection grid.
ugrid_lower,newlons = shiftgrid(180.,u_plot_lower,lons,start=False) 
vgrid_lower,newlons = shiftgrid(180.,v_plot_lower,lons,start=False)
ugrid_upper,newlons = shiftgrid(180.,u_plot_upper,lons,start=False) 
vgrid_upper,newlons = shiftgrid(180.,v_plot_upper,lons,start=False)
ugrid_mid,newlons = shiftgrid(180.,u_plot_mid,lons,start=False) 
vgrid_mid,newlons = shiftgrid(180.,v_plot_mid,lons,start=False)




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
if diagnostic == 1:
   # Plot temperature and wind at either the upper or lower level
   putstuffhere = 1
   
else:
   # Plot geostrophic vorticity and wind at either the upper or lower level
   CS1 = m.contourf(x,y,var_upper,cmap=plt.cm.hot_r,levels=cflevs_upper, extend='both',zorder=1)
   cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
   CS2 = m.contour(x, y, var_upper, cflevs_upper, colors='k', linewidths=0.5)
   cbar.set_ticks(cflevs_upper_ticks)
   clabs = ['%i K' % f for f in cflevs_upper_ticks]
   cbar.ax.set_yticklabels(clabs, size=10)      
if plot_barbs == "true":
   # Rotate and interpolate wind from lat/lon grid to map projection grid.
   uproj_lower,vproj_lower, xx_lower, yy_lower = m.transform_vector(ugrid_lower,vgrid_lower,newlons,lats,41,41,returnxy=True)
   uproj_upper,vproj_upper, xx_upper, yy_upper = m.transform_vector(ugrid_upper,vgrid_upper,newlons,lats,41,41,returnxy=True)  
   uproj_mid,vproj_mid, xx_mid, yy_mid = m.transform_vector(ugrid_mid,vgrid_mid,newlons,lats,41,41,returnxy=True)     
   barbs = m.barbs(xx_upper,yy_upper,uproj_upper,vproj_upper,length=5,barbcolor='k',flagcolor='r',linewidth=1.0)
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

if diagnostic == 1:
   # Plot geostrophic vorticity and wind at the middle level 
   
   CS1 = m.contourf(x,y,var_mid,cmap=plt.cm.hot_r,levels=cflevs_mid, extend='both',zorder=1)
   cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
   CS2 = m.contour(x, y, var_mid, cflevs_mid, colors='k', linewidths=0.5)
   cbar.set_ticks(cflevs_mid_ticks)
   clabs = ['%i K' % f for f in cflevs_mid_ticks]      
else:
   # Plot temperature and wind at the middle level
   CS1 = m.contour(x,y,var_mid,levels=cflevs_mid, colors='k', linestyles='solid',linewidths=2.0)
   plt.clabel(CS1, inline=1, fontsize=10, fmt='%i')   

if plot_barbs == "true":
   barbs = m.barbs(xx_mid,yy_mid,uproj_mid,vproj_mid,length=5,barbcolor='k',flagcolor='r',linewidth=1.0)

ax1.set_title(titletext2)
save_name = figname + "_02.png"
plt.savefig(save_name, bbox_inches='tight')


# Figure 3
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

CS = m.contourf(x,y,diff_varadv,cmap=plt.cm.bwr,levels=cflevs_diffadv, extend='both',zorder=1)
cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
cbar.set_ticks(cflevs_diffadv_ticks)
clabs = ['%i K' % f for f in cflevs_diffadv_ticks]
cbar.ax.set_yticklabels(clabs, size=10)   
ax1.set_title(titletext3)
save_name = figname + "_03.png"
plt.savefig(save_name, bbox_inches='tight')


# Figure 4
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

CS = m.contourf(x,y,var_adv_mid,cmap=plt.cm.bwr,levels=cflevs_midadv, extend='both',zorder=1)
cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
if plot_barbs == "true":
    # Rotate and interpolate wind from lat/lon grid to map projection grid.
   barbs = m.barbs(xx_mid,yy_mid,uproj_mid,vproj_mid,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)
cbar.set_ticks(cflevs_midadv_ticks)
clabs = ['%i K' % f for f in cflevs_midadv_ticks]
cbar.ax.set_yticklabels(clabs, size=10)   
ax1.set_title(titletext4)
save_name = figname + "_04.png"
plt.savefig(save_name, bbox_inches='tight')


if diagnostic == 1:
   # Figure 5
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

   CS = m.contourf(x,y,geoph_tend,cmap=plt.cm.bwr,levels=cflevs_geophtend, extend='both',zorder=1)
   cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
   cbar.set_ticks(cflevs_geophtend_ticks)
   clabs = ['%i K' % f for f in cflevs_geophtend_ticks]
   cbar.ax.set_yticklabels(clabs, size=10)
   ax1.set_title(titletext5)
   save_name = figname + "_05.png"
   plt.savefig(save_name, bbox_inches='tight')


if diagnostic == 2:
   # Figure 5
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

   CS = m.contourf(x,y,omegain,cmap=plt.cm.bwr,levels=cflevs_omega, extend='both',zorder=1)
   cbar = plt.colorbar(CS, shrink=0.95, orientation='horizontal',extend='both')
   cbar.set_ticks(cflevs_omega_ticks)
   clabs = ['%i K' % f for f in cflevs_omega_ticks]
   cbar.ax.set_yticklabels(clabs, size=10)
   ax1.set_title(titletext5)
   save_name = figname + "_05.png"
   plt.savefig(save_name, bbox_inches='tight')


plt.show()
