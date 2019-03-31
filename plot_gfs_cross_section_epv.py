# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

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
date_string = '1900010100' # yyyymmddhh
plot_anomaly = 'true'
record_num = 2
lat_range = [31.0,31.0]
lon_range = [260,285]
#lat_range = [25,40]
#lon_range = [272.5,272.5]
pres_range = [100,1000]
label_fontsize = 14
figname = "gfs_cross_section_analysis_ew" # will save an image named this with a .png extension

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fpath = '/home/xxxx/data'
fname = 'files1234.nc'
###################################
# END user options
###################################


# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, fname)
print fpath


if (lat_range[0] == lat_range[1]):
   cross_orient = 'east-west'
else:
   cross_orient = 'north-south'   



f = netCDF4.Dataset(fpath,'r')
lons = f.variables['lon_0'][:]
lats = f.variables['lat_0'][::-1] # Read in reverse direction
levs = f.variables['lv_ISBL0'][:]
levs_omega = f.variables['lv_ISBL6'][:]
u_plot =  f.variables['UGRD_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
v_plot =  f.variables['VGRD_P0_L100_GLL0'][record_num,:,::-1,:].squeeze() 
temperature = f.variables['TMP_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
omegain = f.variables['VVEL_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
f.close

nz, ny, nx = temperature.shape

pres = np.zeros_like(temperature).astype('f')   
for kk in range(0,len(levs)):      
   pres[kk,:,:] = levs[kk]

theta = wm.temp_to_theta(temperature, pres)

epv = wm.epv_sphere(theta,pres,u_plot,v_plot,lats,lons)


# Convert wind to knots
u_plot = u_plot * 1.94384
v_plot = v_plot * 1.94384

uv_plot = np.sqrt(u_plot**2 + v_plot**2)

if cross_orient == 'east-west':
   latind =  np.ravel(lats==lat_range[0])
   loninds = np.ravel((lons<=lon_range[1])&(lons>=lon_range[0]))
   x_cross = lons[loninds]
   epv_cross = epv[:,latind,loninds].squeeze()
   theta_cross = theta[:,latind,loninds].squeeze()
   omega_cross = omegain[:,latind,loninds].squeeze()
   u_cross = u_plot[:,latind,loninds].squeeze()
   v_cross = v_plot[:,latind,loninds].squeeze()
   uv_cross = uv_plot[:,latind,loninds].squeeze()
   
   if plot_anomaly == 'true':
      var_cross = epv_cross
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
   epv_cross = epv[:,latinds,lonind].squeeze()
   theta_cross = theta[:,latinds,lonind].squeeze()
   omega_cross = omegain[:,latinds,lonind].squeeze()
   u_cross = u_plot[:,latinds,lonind].squeeze()
   v_cross = v_plot[:,latinds,lonind].squeeze()
   uv_cross = uv_plot[:,latinds,lonind].squeeze()

   if plot_anomaly == 'true':
      var_cross = epv_cross
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


if plot_anomaly == 'true':
   cbar_min_epv = -4
   cbar_max_epv = 4.2
   cint_epv = 0.4
else:
   cbar_min_epv = 2
   cbar_max_epv = 10
   cint_epv = 0.5


cflevs_epv =  np.arange(cbar_min_epv, cbar_max_epv, cint_epv)
cflevs_epv_ticks = np.arange(cbar_min_epv,cbar_max_epv,4*cint_epv)

cflevs_trop = [2.0,2.0] 

dflevs_winds = np.arange(20,150,5)
cflevs_temp = np.arange(200,501,3)
cflevs_theta = np.arange(250,501,3)
#cflevs_omega = np.arange(-0.6,0.61,0.06)
cflevs_omega = np.arange(-1.4,-0.06,0.05)


levs_hPa = levs/100
levs_hPa_omega = levs_omega/100

fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
if plot_anomaly == 'false':
   CS1 = plt.contourf(x_cross,levs_hPa,epv_cross,cmap=plt.cm.winter,levels=cflevs_epv)
else:
   CS1 = plt.contourf(x_cross,levs_hPa,var_cross,cmap=plt.cm.bwr,levels=cflevs_epv)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
clabs = ['%i K' % f for f in cflevs_epv_ticks]
cbar.ax.set_yticklabels(clabs, size=10) 

CS2 = plt.contour(x_cross,levs_hPa,epv_cross,levels=cflevs_trop,colors='k', linestyles='solid',linewidths=2.5)  
plt.clabel(CS2, inline=1, fontsize=10, fmt='%i')   

CS3 = plt.contour(x_cross,levs_hPa,theta_cross,levels=cflevs_theta,colors='k', linestyles='solid',linewidths=1.5)  
plt.clabel(CS3, inline=1, fontsize=10, fmt='%i')   

CS4 = plt.contour(x_cross,levs_hPa_omega,omega_cross,levels=cflevs_omega,colors='r', linewidths=1.5)  
#plt.clabel(CS4, inline=1, fontsize=10, fmt='%i')   

plt.ylim((pres_range[0],pres_range[1]))
ax1.set_ylim(ax1.get_ylim()[::-1])
#ax1.set_yscale('log')
if cross_orient == 'east-west':
   plt.xlabel('Longitude (Degrees)',fontsize=label_fontsize)
else:
   plt.xlabel('Latitude (Degrees North)',fontsize=label_fontsize)
plt.ylabel('Pressure (hPa)',fontsize=label_fontsize)

save_name = figname + ".png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()
