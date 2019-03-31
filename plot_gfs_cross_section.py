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

# from mstats import *
# <codecell>

###################################
# Set user options
###################################
#date_string = '2012081112' # yyyymmddhh
date_string = '2012090800' # yyyymmddhh
#varname2fill = 'TMP_P0_L100_GLL0' # Pascals; set to -1 for sea level pressure
varname2fill ='HGT_P0_L100_GLL0'
varname2cntr ='TMP_P0_L100_GLL0'
plot_theta = 'true'
plot_temperature = 'true'
plot_anomaly = 'true'
plot_winds = 'false' # 'true' to plot windbarbs, 'false' otherwise
plot_contours = 'true' # True to plot contours, false otherwise
lat_range = [37,37]
lon_range = [260,280]
pres_range = [200,1000]
figname = "gfs_cross_section_analysis_ew" # will save an image named this with a .png extension

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/xxxx/data'

###################################
# END user options
###################################

# <codecell>

# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, 'gfs_4_%s_%s00_000.nc' % (
                     dt.strftime('%Y%m%d'), dt.strftime('%H')))
print fpath


if (lat_range[0] == lat_range[1]):
   cross_orient = 'east-west'
else:
   cross_orient = 'north-south'

# If plotting geopotential height, automatically plot anomaly instead since the full field varies strongly with height
if varname2fill == 'HGT_P0_L100_GLL0':
   plot_anomaly = 'true'

f = netCDF4.Dataset(fpath,'r')
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]
var2fill = f.variables[varname2fill][:,::-1,:]
var2cntr = f.variables[varname2cntr][:,::-1,:]
u_plot =  f.variables['UGRD_P0_L100_GLL0'][:,::-1,:]
v_plot =  f.variables['VGRD_P0_L100_GLL0'][:,::-1,:]
if ( plot_theta == 'true' or plot_temperature == 'true' ):
   temperature = f.variables['TMP_P0_L100_GLL0'][:,::-1,:]

f.close

nx = np.size(lons)
ny = np.size(lats)


#if plot_winds == 'true':
uv_plot = np.sqrt(u_plot**2 + v_plot**2)


# Convert temperature to potential temperature
if plot_theta == 'true' :
   pres = np.zeros_like(temperature).astype('f')
   for kk in range(0,len(levs)):
      pres[kk,:,:] = levs[kk]

   theta = wm.temp_to_theta(temperature, pres)


# Convert wind to knots
u_plot = u_plot * 1.94384
v_plot = v_plot * 1.94384

if cross_orient == 'east-west':
   latind =  np.ravel(lats==lat_range[0])
   loninds = np.ravel((lons<=lon_range[1])&(lons>=lon_range[0]))
   x_cross = lons[loninds]
   var_cross = var2fill[:,latind,loninds].squeeze()
   var_cross_cntr = var2cntr[:,latind,loninds].squeeze()
   u_cross = u_plot[:,latind,loninds].squeeze()
   v_cross = v_plot[:,latind,loninds].squeeze()
   uv_cross = uv_plot[:,latind,loninds].squeeze()
   if plot_theta == 'true':
      theta_cross = theta[:,latind,loninds].squeeze()
   if plot_temperature == 'true':
      temperature_cross = temperature[:,latind,loninds].squeeze()

   if plot_anomaly == 'true':
      varcross_mean = np.average(var_cross,1)
      tmp = np.zeros_like(var_cross).astype('f')
      for kk in range(0,len(levs)):
         tmp[kk,:] = var_cross[kk,:] - varcross_mean[kk]
      del var_cross
      var_cross = tmp
      del tmp

else:
   lonind =  np.ravel(lons==lon_range[0])
   latinds = np.ravel((lats<=lat_range[1])&(lats>=lat_range[0]))
   x_cross = lats[latinds]
   var_cross = var2fill[:,latinds,lonind].squeeze()
   var_cross_cntr = var2cntr[:,latinds,lonind].squeeze()
   u_cross = u_plot[:,latinds,lonind].squeeze()
   v_cross = v_plot[:,latinds,lonind].squeeze()
   uv_cross = uv_plot[:,latinds,lonind].squeeze()
   if plot_theta == 'true':
      theta_cross = theta[:,latinds,lonind].squeeze()
   if plot_temperature == 'true':
      temperature_cross = temperature[:,latinds,lonind].squeeze()
   if plot_anomaly == 'true':
      varcross_mean = np.average(var_cross,1)
      tmp = np.zeros_like(var_cross).astype('f')
      for kk in range(0,len(levs)):
         tmp[kk,:] = var_cross[kk,:] - varcross_mean[kk]
      del var_cross
      var_cross = tmp
      del tmp



# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)

if varname2fill == 'TMP_P0_L100_GLL0':
   base_cntr = 275 # Base contour
   cint = 5 # Contour interval
   nconts = 15 # number of contours
if varname2fill == 'HGT_P0_L100_GLL0':
   cint = 10 # Contour interval
if plot_anomaly == 'true':
   if cross_orient == 'east-west':
      cint = 5
      base_cntr = 0
      nconts = 10
   else:
      cint = 10
      base_cntr = 0
      nconts = 10

cbar_min = base_cntr-nconts*cint
cbar_max = base_cntr+nconts*cint
cflevs = np.arange(cbar_min, cbar_max+1, cint)
cnlevs = np.arange(cbar_min,cbar_max,4*cint)

# If the color interval is outside the range of data, limit the data displayed.
#   Otherwise, it will shade white which is confusing because it also means zero.
ncy,ncx = np.shape(var_cross)
for jj in range(0,ncy):
   for ii in range(0,ncx):
      if var_cross[jj,ii]>cbar_max:
         var_cross[jj,ii] = cbar_max
      if var_cross[jj,ii]<cbar_min:
         var_cross[jj,ii] = cbar_min



cflevs_cntr = np.arange(200,501,4)

cflevs_winds = np.arange(20,150,5)
cflevs_temp = np.arange(200,501,3)
cflevs_theta = np.arange(250,501,3)

ncolors = np.size(cflevs)-1
djet1 = plt.get_cmap('jet')
djet2 = um.cmap_whitezero(djet1,ncolors,1,-1)
if plot_anomaly == 'false':
   djet = djet1
else:
   djet = djet2

label_fontsize = 14


levs_hPa = levs/100

fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
cf = plt.contourf(x_cross,levs_hPa,var_cross,cmap=djet,levels=cflevs)
#cf2 = plt.contour(x_cross,levs_hPa,var_cross_cntr,levels=cflevs_cntr,colors='k', linestyles='solid',linewidths=1.5)
if plot_winds == 'true':
   cf3 = plt.contour(x_cross,levs_hPa,uv_cross,levels=cflevs_winds,colors='k', linestyles='dashed',linewidths=1.0)
   plt.clabel(cf3, inline=1, fontsize=10, fmt='%i')
if plot_theta == 'true':
   cf4 = plt.contour(x_cross,levs_hPa,theta_cross,levels=cflevs_theta,colors='k', linestyles='solid',linewidths=1.5)
   plt.clabel(cf4, inline=1, fontsize=10, fmt='%i')
if plot_temperature == 'true':
   cf5 = plt.contour(x_cross,levs_hPa,temperature_cross,levels=cflevs_temp,colors='r', linestyles='dashed',linewidths=1.5)
   plt.clabel(cf5, inline=1, fontsize=10, fmt='%i')

plt.ylim((pres_range[0],pres_range[1]))
ax1.set_ylim(ax1.get_ylim()[::-1])
#ax1.set_yscale('log')
if cross_orient == 'east-west':
   plt.xlabel('Longitude (Degrees)',fontsize=label_fontsize)
else:
   plt.xlabel('Latitude (Degrees North)',fontsize=label_fontsize)
plt.ylabel('Pressure (hPa)',fontsize=label_fontsize)

cbar = plt.colorbar(cf, shrink=0.95, orientation='horizontal',extend='both')
cbar.set_ticks(cnlevs)
clabs = ['%i K' % f for f in cnlevs]
cbar.ax.set_yticklabels(clabs, size=10)
save_name = figname + ".png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()
