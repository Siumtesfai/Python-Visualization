#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 17:59:07 2020

@author: gregoryjenkins
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 08:30:05 2020

@author: gregoryjenkins
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
from netCDF4 import Dataset

from wrf import (getvar, to_np, vertcross, smooth2d, CoordPair, GeoBounds,
                 get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim,ALL_TIMES)

# Open the NetCDF file

ncfile = Dataset('E:\Greg-Sium-WRF-CHEM\Runs\R2\wrfout2')

# Get sthe sea level pressure

time =getvar(ncfile,'times',timeidx=ALL_TIMES)
ntimes = np.size(time)
timename = time.astype(str)
for i in range(0,ntimes,1):

	pm25 = getvar(ncfile, 'PM2_5_DRY',timeidx=i)
	z = getvar(ncfile, 'z',timeidx=i)

# Set the start point and end point for the cross section
	start_point = CoordPair(lat=6., lon=-17.4)
	end_point = CoordPair(lat=20., lon=-17.4)
# Set the start point and end point for the cross section


# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section in the metadata by setting latlon
# to True.

	pm25_cross = vertcross(pm25, z, wrfin=ncfile, start_point=start_point,
                    end_point=end_point, latlon=True, meta=True)

# Get the lat/lon points
	lats, lons = latlon_coords(pm25)

# Get the cartopy projection object
	cart_proj = get_cartopy(pm25)


# Create a figure that will have 3 subplots
	fig = plt.figure(figsize=(12,9))

	ax_va = fig.add_subplot()

# Download and create the states, land, and oceans using cartopy features
	states = cfeature.NaturalEarthFeature(category='cultural', scale='50m',
                                      facecolor='none',
                                      name='admin_1_states_provinces_shp')
	land = cfeature.NaturalEarthFeature(category='physical', name='land',
                                    scale='50m',
                                    facecolor=cfeature.COLORS['land'])
	ocean = cfeature.NaturalEarthFeature(category='physical', name='ocean',
                                     scale='50m',
                                     facecolor=cfeature.COLORS['water'])


#levels = [10, 55, 150,250,350,500,1000]
# Make the contour plot for wind speed

	levels = [5,20, 35,55,85,150]


	pm25_contours = ax_va.contourf(to_np(pm25_cross),levels=levels, cmap=get_cmap("YlOrBr"))
	cb_v = fig.colorbar(pm25_contours, ax=ax_va)
	cb_v.ax.tick_params(labelsize=8)

# Make the contour plot for dbz



# Set the x-ticks to use latitude and longitude labels
	coord_pairs = to_np(pm25_cross.coords["xy_loc"])
	x_ticks = np.arange(coord_pairs.shape[0])
	x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]


	ax_va.set_xticks(x_ticks[::6])
	ax_va.set_xticklabels(x_labels[::6], rotation=45, fontsize=6) 


#ax_v.set_xticks(x_ticks[::6])
#ax_v.set_xticklabels(x_labels[::6], rotation=45, fontsize=8) 
# Set the y-ticks to be height
	vert_vals = to_np(pm25_cross.coords["vertical"])
	vert_vals1= vert_vals[:6]
#print(vert_vals)
#print(vert_vals1)
#vert_vals = to_np(wspd_cross.coords[10])
#print(vert_vals)
#v_ticks = np.arange(vert_vals.shape[0])
	v_ticks = np.arange(vert_vals.shape[0])
#print(v_ticks)

	ax_va.set_yticks(v_ticks[::2])
	ax_va.set_yticklabels(vert_vals[::2], fontsize=8)

#ax_v.set_yticks(v_ticks[::20])
#ax_v.set_yticklabels(vert_vals[::20], fontsize=4)
#ax_v.set_ybound(lower = 0,upper=10000)
# Set the x-axis and  y-axis labels

	ax_va.set_xlabel("Latitude, Longitude", fontsize=8)
	ax_va.set_ylabel("Height (m)", fontsize=12)
#        plt.title("Initial 1200 12 Mar Surface PM$_{2.5}$,scaled by (.50)  SLP and barbs (kt) at hour = %s" %timename[i])
	plt.title("Initial 1200 28 Apr vertical cross PM2.5 (mg/m3) at hour =  %s" %timename[i])
	plt.savefig("vertcross_pm25_17w{}.jpg".format(i))
#ax_pm25.set_ylim(0,10000)
#

