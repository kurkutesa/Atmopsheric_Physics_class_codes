#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib as plt
from mpl_toolkits.basemap import Basemap, addcyclic 

def nan2zero(data):
    ''' Convert NaNs to zero '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    temp = np.reshape(data,np.prod(np.size(data)), 1)       
    
    # Find indices with NaNs
    inds = np.argwhere(np.isnan(temp))    
    
    # Replace NaNs with zero
    temp[inds] = 0.                 
    
    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
 
    return data

def zero2nan(data):
    ''' Convert zeros to Nans '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    temp = np.reshape(data,np.prod(np.size(data)), 1)       
    
    # Find indices with NaNs
    inds = np.argwhere(temp==0)    
    
    # Replace zeros with NaNs
    temp[inds] = float('NaN')                 
    
    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
 
    return data

def filter_numeric_nans(data,thresh,repl_val,high_or_low) :
    ''' Filter numerical nans above or below a specified value'''
    ''' '''
    ''' data:        (Input) array to filter '''
    ''' thresh:      (Input) threshold value to filter above or below '''
    ''' repl_val:    (Input) replacement value'''
    ''' high_or_low: (Input)''' 


    dimens = np.shape(data);
    temp = np.reshape(data,np.prod(np.size(data)), 1)
    if high_or_low=='high':        	
	inds = np.argwhere(temp>thresh) 	
	temp[inds] = repl_val	  
    elif high_or_low=='low':    
        inds = np.argwhere(temp<thresh) 
	temp[inds] = repl_val	  
    elif high_or_low =='both':
       	inds = np.argwhere(temp>thresh) 	
	temp[inds] = repl_val
	del inds
       	inds = np.argwhere(temp<-thresh) 	
	temp[inds] = -repl_val	                 
    else:
        inds = np.argwhere(temp>thresh) 
	temp[inds] = repl_val	  

    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
 
    return data    
    
def bold_labels(ax,fontsize=None):
    if fontsize is None:
        fontsize = 14
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

def draw_map_background(m, ax=mpl.gca()):
    ''' Setup the map background '''
    m.drawcoastlines(ax=ax, linewidth=2, color='#444444', zorder=6)
    m.drawcountries(ax=ax, linewidth=1, color='#444444', zorder=5)
    m.drawstates(ax=ax, linewidth=0.66, color='#444444', zorder=4)
    m.drawmapboundary
    
def lonswap(d,subtract=0.):
        sh = np.shape(d)	
	midl = sh[1]/2	
	midl = np.round(midl)
	h=d[:,midl:].copy()-subtract
	d[:,midl:]=d[:,:midl].copy()
	d[:,:midl]=h
	return d
	    
def periodic(d,add=0.):
	return np.append( d, (d[:,0].copy()+add).reshape(-1,1) , 1)    
	
def advance_time(timestrin,timeinc):
    ''' Advances or reverses a time by timeinc'''
    ''' '''
    ''' timestrin: (Input) time string in yyyymmddhh format'''
    ''' timeinc:   (Input) number of hours to increment or decrement'''
    '''             Use a negative sign to decrement '''
    
    import datetime
        
    yyyy = timestrin[0:4]
    mm = timestrin[4:6]
    dd = timestrin[6:8]
    hh = timestrin[8:10]	
	
    date=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
    date += datetime.timedelta(hours=timeinc)
    tt = date.timetuple()
    
    yyyy = str(tt[0])
    mm = str(tt[1])
    dd = str(tt[2])
    hh = str(tt[3])
    
    if tt[0]<1000: yy = '0'+mm 
    if tt[1]<10: mm = '0'+mm 
    if tt[2]<10: dd = '0'+dd
    if tt[3]<10: hh = '0'+hh
    
    timestrin = yyyy+mm+dd+hh        
    
    return timestrin
def get_cmap_cust():
    ''' Setup a custom colortable. '''
    cdict = {'red': ((0.00, 240/255., 220/255.),
                         (0.25, 40/255., 20/255.),
                         (0.50, 225/255., 255/255.),
                         (0.75, 150/255., 150/255.),
                         (1.00, 255/255., 255/255.)),

             'green': ((0.00, 240/255., 220/255.),
                         (0.25, 0/255., 50/255.),
                         (0.50, 255/255., 255/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.)),

             'blue': ((0.00, 255/255., 255/255.),
                         (0.25, 160/255., 150/255.),
                         (0.50, 255/255., 170/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.))}
    return plt.colors.LinearSegmentedColormap('cool2warm', cdict, 256)
def cmap_discretize(cmap, N):
   """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)
   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_whitezero(cmap,N,Nwhites,pos):
   """Whites out middle index of a colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
	Nwhites: Number of whites
	pos: Position for white bar; if -1, will place in middle
	                             if  0, will place at bottom

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   if ( pos == -1 ):
      mid = np.round(N/2)
      mid = int(mid)
   else:
      mid = pos
   
   nadd = Nwhites - 1
   
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)

   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       colors[mid] = 1.
       isodd = 0
       if ( np.mod(N,2) == 0 ) :           
	  colors[mid-1] = 1.
	  isodd = 1
       kk=mid-nadd
       kend=mid+nadd 
       if ( kk < kend ) :      
          while (kk <= kend) :
             colors[kk] = 1.
             kk += 1
       if (isodd == 1 ): colors[kk] = 1.   	  
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors       
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)
