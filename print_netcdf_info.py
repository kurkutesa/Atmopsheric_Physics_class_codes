# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/python

# imports
import netCDF4

import os, datetime


###################################
# Set user options
###################################
fpath = '/home/scavallo/data'
date_string = '2011120700' # yyyymmddhh
###################################
# END user options
###################################

# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, 'gfs_4_%s_%s00_000.nc' % (
                     dt.strftime('%Y%m%d'), dt.strftime('%H')))
print fpath
f = netCDF4.Dataset(fpath,'r')
plotvar = f.variables['HGT_P0_L100_GLL0'][:]

# Print the netcdf file global attributes
for name in f.ncattrs():
   print 'Global attr', name, '=', getattr(f,name)

# Print the variable name, description, and dimensions
for dim in f.variables:
   print dim, f.variables[dim].long_name, f.variables[dim].shape


# Print elements of variable name `lv_ISBL0'
print f.variables['lv_ISBL0'][:]

# print the shape of the field we named `plotvar'
print 'shape (dimensions) of array plotvar is ',plotvar.shape
