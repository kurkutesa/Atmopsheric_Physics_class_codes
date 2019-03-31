#!/usr/bin/python

import numpy as np
from scipy import ndimage

Cp = 1004.5;
Cv = 717.5;
Rd = 287.04;
Rv = 461.6;
RvRd = Rv / Rd;
g = 9.81;
L = 2.50e6;
Talt = 288.1500;
Tfrez = 273.1500;
To = 300;
Po = 101325;
Pr = 1000.;
lapsesta = 6.5 / 1000;
kappa = Rd / Cp;
epsil = Rd/Rv;
pi = 3.14159265;
pid = pi/180;
R_earth = 6371200;
omeg_e = (2*pi) / (24*3600);
eo = 6.11;
missval = -9999;
eps = 2.2204e-16

def temp_to_theta(temp, pres):
    ''' Compute potential temperature '''
    ''' '''
    ''' theta: Input potential temperature (K) '''
    ''' pres:  Input pressure (Pa)'''
    ''' temp:  Output temperature (K)'''
    return temp * (100000. / pres) ** 0.286

def theta_to_temp(theta, pres):
    ''' Compute temperature '''
    ''' '''
    ''' temp:  Input temperature (K)'''
    ''' pres:  Input pressure (Pa)'''
    ''' theta: Output potential temperature (K)'''
    return theta * (pres / 100000.) ** 0.286

def td_to_mixrat(tdew, pres):
    ''' Convert from dewpoint temperature to water vapor mixing ratio '''
    ''' '''
    ''' tdew:   Input dewpoint temperature (K)'''
    ''' pres:   Input pressure (Pa)'''
    ''' mixrat: Output water vapor mixing ratio (kg/kg)'''
    pres = pres/100
    mixrat = eo / (pres * RvRd)  * np.exp( (L/Rv)*((1/Tfrez) - (1 / tdew) ) )
    return mixrat

def mixrat_to_td(qvap, pres):
    ''' Convert from water vapor mixing ratio to dewpoint temperature '''
    ''' '''
    ''' qvap: Input water vapor mixing ratio (kg/kg)'''
    ''' pres: Input pressure (Pa)'''
    ''' tdew: Output dewpoint temperature (K)'''
    pres = pres/100
    evap = qvap * pres * RvRd;
    tdew = 1/((1/Tfrez) - (Rv/L)*np.log(evap/eo))
    return tdew

def spechum_to_td(spechum, pres):
    ''' Convert from water vapor mixing ratio to dewpoint temperature '''
    ''' '''
    ''' spechum: Input specific humidity in (kg/kg)'''
    ''' pres: Input pressure (Pa)'''
    ''' tdew: Output dewpoint temperature (K)'''
    qvap = (spechum/(1-spechum))
    
    pres = pres/100
    evap = qvap * pres * RvRd;
    tdew = 1/((1/Tfrez) - (Rv/L)*np.log(evap/eo))
    return tdew
    
def claus_clap(temp):
    ''' Compute saturation vapor pressure '''
    ''' '''
    ''' temp: Input temperature (K)  '''
    ''' esat: Output satuation vapor pressure (Pa)'''
    esat = (eo * np.exp( (L / Rv) * ( 1/Tfrez - 1/temp) ) ) * 100
    return esat

def thetae(thta, temp, esat):
    ''' Compute equivalent potential temperature '''
    ''' '''
    ''' thta:   Input potential temperature (K) '''
    ''' temp:   Input temperature (K) '''
    ''' esat:   Input saturation vapor pressure (Pa)'''
    ''' thetae: Output equivalent potential temperature (K)'''
    thout = thta * np.exp( (L * esat) / (Cp * temp) )
    return thout

def calc_gradient(fldin, dx, dy, dz):
    '''
    Computes the horizontal gradient of any scalar given a constant
    grid spacing in the x, y, and z directions.

    fldin: Input scalar
    dx: Input x grid spacing (must be single value)
    dy: Input y grid spacing (must be single value)
    dz: Input z grid spacing (must be single value)
    '''
    dfdx, dfdy, dfdz = np.gradient(fldin, dx, dy, dz)
    return dfdx, dfdy, dfdz

def latlon_to_dlatdlon(lats,lons):
    """
    Return arrays with the spacing between latitude and longitudes

    The gradients are computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    lats : vector of latitudes in degrees
    lons : vector of longitudes in degrees

    Returns
    -------
    dlat : array with differences in latitudes between grid points with size (lats,lons)
    dlon : array with differences in longitudes between grid points with size (lats,lons)

    Examples
    --------
    dlat,dlon = latlon_to_dlatdlon(lats,lons)

    """

    nlat = len(lats)
    nlon = len(lons)
    latarr = np.zeros((nlat,nlon))
    lonarr = np.zeros((nlat,nlon))
    dlatarr = np.zeros((nlat,nlon))
    dlonarr = np.zeros((nlat,nlon))


    for jj in range(0,nlat):
       for ii in range(0,nlon):
          latarr[jj,ii] = lats[jj]
          lonarr[jj,ii] = lons[ii]


    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    otype = latarr.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    dlats = -1*dlats

    for jj in range(0,nlat):
       for ii in range(0,nlon):
          dlonarr[jj,ii] = dlons[ii]
	  dlatarr[jj,ii] = dlats[jj]


    return dlatarr, dlonarr

def gradient_cartesian(f, *varargs):
    """
    Return the gradient of an N-dimensional array on an evenly spaced grid.

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : N-dimensional array containing samples of a scalar function.
        If 2-D, must be ordered as f(y,x)
	If 3-D, must be ordered as f(z,y,x) or f(p,y,x)
    `*varargs` : scalars
          0, 1, or N scalars specifying the sample distances in each direction,
          that is: `dz`, `dy`, `dx`, ... The default distance is 1.

	  If a vector is specified as the first argument of three, then the difference
	     of this vector will be taken here.


    Returns
    -------
    g : ndarray
      N arrays of the same shape as `f` giving the derivative of `f` with
      respect to each dimension.

    Examples
    --------
    temperature = temperature(pressure,y,x)
    levs = pressure vector
    dy = scalar or array of grid spacing in y direction
    dx = scalar or vector of grid spacing in x direction
    >>> dfdz, dfdy = gradient_evenspaced(fldin, dy, dx)
    >>> dfdz, dfdy, dfdx = gradient_evenspaced(fldin, dz, dy, dx)
    >>> dfdp, dfdy, dfdx = gradient_evenspaced(fldin, levs, dy, dx)
    """
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining arguments")

    df = np.gradient(f)

    if n == 1:
        dy = argsin[0]

	dfdy = df[0]
    elif n == 2:
        dy = argsin[0]
	dx = argsin[1]

	dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
	dy = argsin[1]
        dx = argsin[2]

	dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'


    try:
       M = len(dx.shape)
    except:
       M = 1
    dyarr = np.zeros_like(f).astype(otype)
    dxarr = np.zeros_like(f).astype(otype)
    if M == 1:
       dyarr[:] = dy
       dxarr[:] = dx
       if N == 1:
          ny = np.shape(f)
       elif N == 2:
          ny, nx = np.shape(f)
       else:
          nz, ny, nx = np.shape(f)

    else:
       if N == 1:
          ny = np.shape(f)
	  for jj in range(0,ny):
	     dyarr[jj,ii] = dy[jj]
       elif N == 2:
	  ny, nx = np.shape(f)
	  for jj in range(0,ny):
             for ii in range(0,nx):
        	dyarr[jj,ii] = dy[jj]
		dxarr[jj,ii] = dx[ii]
       else:
	  nz, ny, nx = np.shape(f)
	  for kk in range(0,nz):
             for jj in range(0,ny):
        	for ii in range(0,nx):
                   dyarr[kk,jj,ii] = dy[jj]
		   dxarr[kk,jj,ii] = dx[ii]

    if n==1:
       dfdy = dfdy/dx

       return dfdy
    elif n==2:
       dfdy = dfdy/dy
       dfdx = dfdx/dx

       return dfdy,dfdx
    elif n==3:
       dfdy = dfdy/dy
       dfdx = dfdx/dx

       nzz = np.shape(levs)
       
       if not nzz:
          nzz=0

       if nzz>1:
	    zin = levs
	    dz = np.zeros_like(zin).astype(otype)
            dz[1:-1] = (zin[2:] - zin[:-2])/2
            dz[0] = (zin[1] - zin[0])
            dz[-1] = (zin[-1] - zin[-2])	    	    	    
	    if zin[1] < zin[0]:
	       dz = dz*-1 # assume the model top is the first index and the lowest model is the last index
	    
	    dx3 = np.ones_like(f).astype(otype)
	    for kk in range(0,nz):
	       dx3[kk,:,:] = dz[kk]
       else:
            dx3 = np.ones_like(f).astype(otype)
            dx3[:] = dx[0]

       dfdz = dfdz/dx3
       return dfdz,dfdy,dfdx

def gradient_sphere(f, *varargs):
    """
    Return the gradient of a 2-dimensional array on a sphere given a latitude
    and longitude vector.

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : A 2-dimensional array containing samples of a scalar function.
    latvec: latitude vector
    lonvec: longitude vector

    Returns
    -------
    g : dfdx and dfdy arrays of the same shape as `f` giving the derivative of `f` with
        respect to each dimension.

    Examples
    --------
    temperature = temperature(pressure,latitude,longitude)
    levs = pressure vector
    lats = latitude vector
    lons = longitude vector
    >>> tempin = temperature[5,:,:]
    >>> dfdlat, dfdlon = gradient_sphere(tempin, lats, lons)

    >>> dfdp, dfdlat, dfdlon = gradient_sphere(temperature, levs, lats, lons)

    based on gradient function from /usr/lib64/python2.6/site-packages/numpy/lib/function_base.py
    """

    R_earth = 6371200.
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining argumens")

    df = np.gradient(f)

    if n == 1:
        lats = argsin[0]

	dfdy = df[0]
    elif n == 2:
        lats = argsin[0]
	lons = argsin[1]

	dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
	lats = argsin[1]
        lons = argsin[2]

	dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    latarr = np.zeros_like(f).astype(otype)
    lonarr = np.zeros_like(f).astype(otype)
    if N == 1:
       nlat = np.shape(f)
       for jj in range(0,nlat):
          latarr[jj,ii] = lats[jj]
       lonarr = latarr
       lons = lats
    elif N == 2:
       nlat, nlon = np.shape(f)
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             latarr[jj,ii] = lats[jj]
	     lonarr[jj,ii] = lons[ii]
    else:
       nz, nlat, nlon = np.shape(f)
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                latarr[kk,jj,ii] = lats[jj]
		lonarr[kk,jj,ii] = lons[ii]

    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    outvals = []

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    #dlats = -1*dlats

    dlatarr = np.tile(dlats,[nlon,1])
    dlatarr = np.reshape(dlatarr,[nlat,nlon])

    dlonarr = np.zeros_like(f).astype(otype)
    if N==2:
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             dlonarr[jj,ii] = dlons[ii]
    elif N==3:
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
	        dlonarr[kk,jj,ii] = dlons[ii]

    dlatsrad = dlatarr*(pi/180)
    dlonsrad = dlonarr*(pi/180)
    latrad = latarr*(pi/180)

    if n==1:
       dx1 = R_earth * dlatsrad
       dfdy = dfdy/dx1

       return dfdy
    elif n==2:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       return dfdy,dfdx
    elif n==3:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       nzz = np.shape(levs)
       if not nzz:
          nzz=0

       if nzz>1:
	    zin = levs
	    dz = np.zeros_like(zin).astype(otype)
            dz[1:-1] = (zin[2:] - zin[:-2])/2
            dz[0] = (zin[1] - zin[0])
            dz[-1] = (zin[-1] - zin[-2])
	    if zin[0,1,1] > zin[1,1,1]:
	       dz = dz*-1 # assume the model top is the first index and the lowest model is the last index

	    dx3 = np.ones_like(f).astype(otype)
	    for kk in range(0,nz):
	       dx3[kk,:,:] = dz[kk]
       else:
            dx3 = np.ones_like(f).astype(otype)
            dx3[:] = dx[0]

       dfdz = dfdz/dx3
       return dfdz,dfdy,dfdx





def _get_gradients(u, v, dx, dy):
    #Helper function for getting convergence and vorticity from 2D arrays
    dudx, dudy = np.gradient(u, dx, dy)
    dvdx, dvdy = np.gradient(v, dx, dy)
    return dudx, dudy, dvdx, dvdy

def vertical_vorticity(u, v, dx, dy, grid_opt):
    '''
Calculate the vertical vorticity of the horizontal wind. The grid
must have a constant spacing in each direction.

u, v : 2 dimensional arrays
Arrays with the x and y components of the wind, respectively.
X must be the first dimension and y the second.

dx : scalar or array
The grid spacing in the x-direction

dy : scalar or array
The grid spacing in the y-direction

grid_opt: 1 for cartesian grid
          2 for lat/lon grid

Returns : 2 dimensional array
The vertical vorticity
'''

    if grid_opt == 1:
       dudy,dudx = gradient_cartesian(u, dy, dx)
       dvdy,dvdx = gradient_cartesian(v, dy, dx)
    else:
       dudy,dudx = gradient_sphere(u, dy, dx)
       dvdy,dvdx = gradient_sphere(v, dy, dx)

    return dvdx - dudy

def h_convergence(u, v, dx, dy):
    '''
Calculate the horizontal convergence of the horizontal wind. The grid
must have a constant spacing in each direction.

u, v : 2 dimensional arrays
Arrays with the x and y components of the wind, respectively.
X must be the first dimension and y the second.

dx : scalar
The grid spacing in the x-direction

dy : scalar
The grid spacing in the y-direction

Returns : 2 dimensional array
The horizontal convergence
'''
    dudx, dudy, dvdx, dvdy = _get_gradients(u, v, dx, dy)
    return dudx + dvdy

def geostrophic_latlon(u, v, ghgt, lats, lons):
    '''
Calculate horizontal advection on a latitude/longitude grid

Input:
   
   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                             dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
			     components of the wind, respectively.
			     
   ghgt(lats,lons): 2 dimensional array of geopotential height
   
   lats(lats) : latitude vector
   lons(lons) : longitude array

Output:

   ug(lats,lons), vg(lats,lons): 2 dimensional geostrophic u and v 
                             wind arrays dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
			     components of the geostrophic wind, 
			     respectively.
'''
    # 2D latitude array
    glats = np.zeros_like(u).astype('f')      
    for jj in range(0,len(lats)):
	for ii in range(0,len(lons)):    
	   glats[jj,ii] = lats[jj]

    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    

    ghgt = ndimage.gaussian_filter(ghgt,0.75)
    
    geop = ghgt * 9.81
    dphidy,dphidx = gradient_sphere(geop, lats, lons)
    
    ug = -dphidy/f
    vg = dphidx/f
    
    return ug, vg

def hadvection_latlon(u, v, datain, lats, lons):
    '''
Calculate horizontal advection on a latitude/longitude grid

Input:

   datain(lats,lons): 2 dimensional array to compute advection of

   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                             dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
			     components of the wind, respectively.

   lats(lats) : latitude vector
   lons(lons) : longitude array

Output:

   dataout(lats,lons): Two dimensional array with the
                        horizontal advection of data
'''
    datady,datadx = gradient_sphere(datain, lats, lons)
    dataout = -u*datadx -v*datady
    
    return dataout
def vertical_vorticity_latlon(u, v, lats, lons, abs_opt):
    '''
Calculate the vertical vorticity on a latitude/longitude grid

Input:
   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                                dimensioned by (lats,lons).
                                Arrays correspond to the x and y 
		   	        components of the wind, respectively.

   lats(lats) : latitude vector
   lons(lons) : longitude array

   abs_opt: 1 to compute absolute vorticity
            0 for relative vorticity only

Output: 

   vert_vort(lats,lons): Two dimensional array of vertical voriticity
'''
    dudy,dudx = gradient_sphere(u, lats, lons)
    dvdy,dvdx = gradient_sphere(v, lats, lons)
    
    if abs_opt == 1 :
       # 2D latitude array
       glats = np.zeros_like(u).astype('f')      
       for jj in range(0,len(lats)):
	   for ii in range(0,len(lons)):    
	      glats[jj,ii] = lats[jj]

       # Coriolis parameter
       f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    
       
    else:   
       f = 0.
    
    zeta = dvdx - dudy
    vert_vort = zeta + f  

    return vert_vort

def vertical_vorticity_cartesian(u, v, lats, deltax, deltay, abs_opt):
    '''
Calculate the vertical vorticity on a Cartesian grid

Input:
   u(y,x), v(y,x) : 2 dimensional u and v wind arrays                                 
                    Arrays correspond to the x and y 
		    components of the wind, respectively.

   lats(lats) : 2D latitude array
   deltax     : horizontal x grid spacing in meters
   deltay     : horiztional y grid spacing in meters

   abs_opt: 1 to compute absolute vorticity
            0 for relative vorticity only

Output: 

   vert_vort(y,x): Two dimensional array of vertical voriticity
'''
    dudy,dudx = gradient_cartesian(u, deltax, deltay)
    dvdy,dvdx = gradient_cartesian(v, deltax, deltay)
        
    iy, ix = u.shape
    
    if abs_opt == 1 :
       # Coriolis parameter
       f = 2*(7.292e-05)*np.sin(np.deg2rad(lats))    
       
    else:   
       f = 0.
    
    zeta = dvdx - dudy
    vert_vort = zeta + f  

    return vert_vort

def thermal_wind_sphere(thickness_in, lats, lons):
    '''
Calculate the thermal wind on a latitude/longitude grid

Input:
   thickness_in(lats,lons) : 2 dimensional geopotential height thickness array 
                             dimensioned by (lats,lons).


   lats(lats) : latitude vector
   lons(lons) : longitude array


Output: 

   thermal_wind_u(lats,lons), thermal_wind_v(lats,lons): Two dimensional arrays of 
                                                         u- and v- components of 
							 thermal wind vector
'''
    
    # Smooth the thicknesses
    thickness_in = ndimage.gaussian_filter(thickness_in,0.75)
    
    dthickdy,dthickdx = gradient_sphere(thickness_in, lats, lons)

       
    # 2D latitude array
    glats = np.zeros_like(thickness_in).astype('f')      
    for jj in range(0,len(lats)):
	for ii in range(0,len(lons)):    
	   glats[jj,ii] = lats[jj]

    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    
    
           
    thermal_wind_u = -1*(9.81/f)*dthickdy
    thermal_wind_v = (9.81/f)*dthickdx

    return thermal_wind_u, thermal_wind_v

def epv_sphere(theta,pres,u,v,lats,lons):
    """

   Computes the Ertel Potential Vorticity (PV) on a latitude/longitude grid

   Input:    

       theta:       3D potential temperature array on isobaric levels
       pres:        3D pressure array
       u,v:         3D u and v components of the horizontal wind on isobaric levels
       lats,lons:   1D latitude and longitude vectors

   Output:
      
      epv: Ertel PV in potential vorticity units (PVU)
   

   Steven Cavallo
   October 2012
   University of Oklahoma
    
    """
    iz, iy, ix = theta.shape
    

    dthdp, dthdy, dthdx = gradient_sphere(theta, pres, lats, lons)
    dudp, dudy, dudx = gradient_sphere(u, pres, lats, lons)
    dvdp, dvdy, dvdx = gradient_sphere(v, pres, lats, lons)    

    avort = np.zeros_like(theta).astype('f')   
    for kk in range(0,iz):       
       avort[kk,:,:] = vertical_vorticity_latlon(u[kk,:,:].squeeze(), v[kk,:,:].squeeze(), lats, lons, 1)

    epv = (-9.81*(-dvdp*dthdx - dudp*dthdy + avort*dthdp))*10**6


    return epv
    
def epv_cartesian(theta,pres,u,v,lats,deltax,deltay):
    """

   Computes the Ertel Potential Vorticity (PV) on a Cartesian grid

   Input:    

       theta:       3D potential temperature array on isobaric levels
       pres:        1D pressure vector
       u,v:         3D u and v components of the horizontal wind on isobaric levels
       lats:        2D latitude array
       deltax, deltay: x and y horizontal grid spacing in meters

   Output:
      
      epv: Ertel PV in potential vorticity units (PVU)
   

   Steven Cavallo
   October 2012
   University of Oklahoma
    
    """
    iz, iy, ix = theta.shape
    
    dthdp, dthdy, dthdx = gradient_cartesian(theta, pres, deltax, deltay)
    dudp, dudy, dudx = gradient_cartesian(u, pres, deltax, deltay)
    dvdp, dvdy, dvdx = gradient_cartesian(v, pres, deltax, deltay)

    avort = np.zeros_like(theta).astype('f')   
    for kk in range(0,iz):       
       avort[kk,:,:] = vertical_vorticity_cartesian(u[kk,:,:].squeeze(), v[kk,:,:].squeeze(), lats, deltax, deltay, 1)

    epv = (-9.81*(-dvdp*dthdx - dudp*dthdy + avort*dthdp))*10**6


    return epv


    
def interp2pv(pv, fval, pv_surf):
    """

   Linearly interpolates a field to a PV surface

   Input:    

       pv - 3D array that contains the PV at common levels (P or theta)
       fval - 3D field to interpolate
       pv_surf - potential vorticity surface to interpolate onto

    Steven Cavallo
    October 2012
    University of Oklahoma
    
    """
    iz, iy, ix = pv.shape
    
    # Scan  from the top of the model downward.
    # The zeroth index is assumed to correspond to the top of the model.
    trop = fval[0,:,:].squeeze() 

    for jj in range(iy):
       for ii in range(ix):  

           aa = np.ravel(pv[:,jj,ii]>pv_surf)
	   pvcol = pv[:,jj,ii].squeeze()
	   minpv = np.min(pvcol)

           if ( minpv >= pv_surf ):
	      # If there are no PV values in the column less than what is desired to interpolate onto, then use value closest to the surface
              trop[jj,ii] = fval[-1,jj,ii]
           elif ( pv[0,jj,ii] <= pv_surf ):
	      # If PV at the model top is less than what is desired to interpolate onto, then use the value at the top of the model
              trop[jj,ii] = fval[0,jj,ii]
           else:               
	       for kk in range(1,iz+1):      
	          # linearly interpolate between the closest levels
	          if pv[kk,jj,ii] < pv_surf:
                      m = (fval[kk-1,jj,ii] - fval[kk,jj,ii]) / (pv[kk-1,jj,ii] - pv[kk,jj,ii])
                      trop[jj,ii] = m * (pv_surf - pv[kk,jj,ii]) + fval[kk,jj,ii]
                      break


    return trop
