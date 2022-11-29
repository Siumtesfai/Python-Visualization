#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 17:48:07 2020

@author: 
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import xarray as xr
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as feature
import numpy as np
from wrf import (to_np, getvar,interplevel, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords,ALL_TIMES)

# Open the NetCDF file

#ncfile = xr.open_dataset('/joko/s0/gsj1/build2/WRF-3.7/test/em_real/wrfout_d01_2021-11-11_12:00:00')#open the file here
#ncfile = Dataset('/joko/s0/gsj1/build2/WRF-3.7/test/em_real/wrfout_d01_2022-04-28_12:00:00')
ncfile = Dataset('E:\Greg-Sium-WRF-CHEM\wrf02')
time =getvar(ncfile,'times',timeidx=ALL_TIMES)
ntimes = np.size(time)
timename = time.astype(str)
hr =0
for i in range(0,ntimes,1):
# Get sthe sea level pressure
    pm25 = getvar(ncfile, 'PM2_5_DRY',timeidx=i)
    pm10 = getvar(ncfile,'PM10',timeidx=i)
    slp = getvar(ncfile, "slp",timeidx=i)
    z = getvar(ncfile, "z", units="dm",timeidx=i)
    ua = getvar(ncfile, "ua", units="kt",timeidx=i)
    va = getvar(ncfile, "va", units="kt",timeidx=i)
    lats, lons = latlon_coords(pm25)
    pm25a = pm25[0,:,:]
    pm10a =pm10[0,:,:]
    ua = ua[0,:,:]
    va = va[0,:,:]
   
# Get the latitude and longitude points
    

# Smooth the sea level pressure since it tends to be noisy near the
# mountains

# Get the latitude and longitude points
    lats, lons = latlon_coords(pm25)


# Get the cartopy mapping object
    cart_proj = get_cartopy(pm25)

# Create a figure
    fig = plt.figure(figsize=(9,5.5))
# Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="white",
                             name="admin_1_states_provinces_shp")
#states = NaturalEarthFeature(scale="50m",
#                             facecolor="black")
                             
    ax.coastlines('50m', linewidth=0.8)
    ax.add_feature(feature.BORDERS)
    colorlist = ['white','yellow','orange','red','brown','purple']
    slevels = [1008,1012,1016,1020,1024,1028,1032]
    clevels = [0, 35, 55, 150]
    #,250,350,450,1000]
# Make the contour outlines and filled contours for the smoothed sea level
# pressure.
# Add the 500 hPa geopotential height contours
#    levels = np.arange(900., 980., 10.)
    test = plt.contour(to_np(lons), to_np(lats), to_np(slp),
                       levels = slevels, colors="black",
                       transform=crs.PlateCarree())
    plt.clabel(test, inline=1, fontsize=10,colors = 'blue')
    
#plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
#    plt.contour(to_np(lons), to_np(lats), to_np(pm25a), contour_levels, colors="black",
#           transform=crs.PlateCarree())
#    pmcontours = plt.contourf(to_np(lons), to_np(lats), to_np(pm25a*.50),levels=clevels,
#             transform=crs.PlateCarree(),
#             colors=colorlist)
    
    pmcontours = plt.contourf(to_np(lons), to_np(lats), to_np(pm10a*.50),levels=clevels,
             transform=crs.PlateCarree(),
             colors=colorlist)
#    pmcontours = plt.contourf(to_np(lons), to_np(lats), to_np(pm25a),level=clevels,
#             transform=crs.PlateCarree(),
#             cmap=get_cmap("YlOrBr"))
    plt.colorbar(pmcontours, ax=ax, shrink=.98)
#plt.colorbar(PM2_5_DRY, ax=ax, shrink=.98)
# Add the 500 hPa wind barbs, only plotting every 125th data point.
    plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
          to_np(ua[::10, ::10]), to_np(va[::10, ::10]),
          transform=crs.PlateCarree(), length=6)
# Add a color bar
#    plt.colorbar(ax=ax, shrink=.98)

# Set the map bounds
#    ax.set_xlim(cartopy_xlim(pm25a))
#    ax.set_ylim(cartopy_ylim(pm25a))
    
# Set the map bounds
    ax.set_xlim(cartopy_xlim(pm10a))
    ax.set_ylim(cartopy_ylim(pm10a))

# Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title("Init 1200 28 Apr Surface PM$_{2.5}$, V3 SLP and barbs (kt) at hour = %s" %timename[i])
    plt.savefig("surpm{}.jpg".format(i))
    plt.show()
