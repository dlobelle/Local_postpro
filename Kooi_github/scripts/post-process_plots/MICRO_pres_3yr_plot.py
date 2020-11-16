#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 08:30:48 2020

@author: Lobel001
"""

# For MICRO presentation: plotting NPacific Ts and trajectories of 10x10 grid run for 3 years
# Input here comes from: plot_Ts_Zs.py and plot_NPac3yrs_Ts_2004_allsizes.py
#ONLY for 920 for now, but could save and run for 30 and 840 too: 

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np 
import cartopy
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from numpy import *
import cmocean.cm as cmo
import os 
import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")


dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/'
dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'

with open(dirread+'90day_global_Ts_for03and07.pickle', 'rb') as f:
    lons,lats,Ts03,Ts07 = pickle.load(f)

# these sims were not run with updated Kooi kernel so results are wrong
# with open(dirread+'3yr_20to28N_Ts_trajectories_for03.pickle','rb') as f:   #and07
#     lons_traj03_S,lats_traj03_S,lons_sink03_S,lats_sink03_S, Ts03_S = pickle.load(f)

with open(dirread+'3yr_28to36N_Ts_trajectories_for03.pickle','rb') as f:  
    lons_all03_N, lats_all03_N, lons_traj03_N,lats_traj03_N,lons_sink03_N,lats_sink03_N, Ts03_N  = pickle.load(f)

  
sizes = [3] #(3,7)

''' prepare cartopy projection and gridspec'''

projection = cartopy.crs.PlateCarree(central_longitude=72+180)
plt.rc('font', size = 18)
fig_w = 10 
fig_h = 10 

fig0 = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
gs0 = fig0.add_gridspec(figure = fig0, nrows = 2, ncols = 1, height_ratios=[5,1])


fig1 = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
gs1 = fig1.add_gridspec(figure = fig1, nrows = 2, ncols = 1, height_ratios=[5,1])


col = 0 
row = 0
si = sizes[row]
cmap0 = plt.cm.get_cmap('magma_r', 10)
#cmap = plt.cm.get_cmap('haline',10) #'magma_r', 10)
cmap = cmo.haline
if col == 0: 
    loc = 'N'
    latmin = 28
    lonmin = -143  
else:
    loc = 'S'   
    latmin = 20
    lonmin = -148

Ts = eval(f'Ts0{si}') 
Ts_3yr = eval(f'Ts0{si}_{loc}')
lons_traj = eval(f'lons_traj0{si}_{loc}')
lats_traj = eval(f'lats_traj0{si}_{loc}')
lon_sink = eval(f'lons_sink0{si}_{loc}')
lat_sink = eval(f'lats_sink0{si}_{loc}')
lons_all = eval(f'lons_all0{si}_{loc}')
lats_all = eval(f'lats_all0{si}_{loc}')


time = np.arange(0,lons_traj.shape[1]/2,0.5)
time_all = np.tile(time,(lons_traj.shape[0],1))    

size_name = '1 mm'
#letter = (chr(ord('a') + row))

c = col
r = row #int(si)-3

#%%
#  Entire Pacific
ax0 = fig0.add_subplot(gs0[r,c], projection=projection)
ax0.coastlines(resolution='50m',zorder=3)
ax0.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
ax0.set_extent([100, -90, -70, 60])
ax0.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap0, alpha = 0.8, vmin = 0, vmax = 90, s = 100,zorder=1,transform = cartopy.crs.PlateCarree()) 
rect0= plt.Rectangle((lonmin,latmin), 10, 10, ec = 'lightblue', lw = 2, fc = 'none', zorder = 3, transform = cartopy.crs.PlateCarree())
ax0.add_patch(rect0)

#%%
Ts[Ts<100]= np.nan # In order for all values below 100 to be white, and highlight the black region
Ts[Ts>=100]= 1000.
# Box with initial locations 
ax = fig1.add_subplot(gs1[r,c], projection=projection)
ax.coastlines(resolution='50m',zorder=3)
ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
ax.set_extent([179, -109, 10, 40])

scat = ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap0, alpha = 0.8, vmin = 0, vmax = 1000, s = 1500,zorder=0,transform = cartopy.crs.PlateCarree()) 
rect= plt.Rectangle((lonmin-2,latmin-2), 12, 12, ec = 'lightblue', lw = 6, fc = 'none', zorder = 3, transform = cartopy.crs.PlateCarree())
ax.add_patch(rect)
ax.scatter(lons_all[:,0],lats_all[:,0], marker='.', color = 'lightblue', cmap = cmap, vmin = 0, vmax = 1000, s = 500,zorder=2,transform = cartopy.crs.PlateCarree())
ax.title.set_text('radius = 1mm')

cbaxes = fig1.add_axes([0.1, 0.3, 0.8, 0.01])
cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='neither')
cbar.set_label(label='Age [days]', size=18)

plt.rc('font', size = 18)
gl1 = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=False, linewidth=0.5,
                      color='gray', alpha=0.5, linestyle='--')
gl1.xlocator = mticker.FixedLocator(np.arange(-180, -110, 30))
gl1.ylocator = mticker.FixedLocator(np.arange(10, 640, 10))
gl1.xlabels_bottom = True
gl1.xformatter = LONGITUDE_FORMATTER
gl1.ylabels_left = True
gl1.yformatter = LATITUDE_FORMATTER

lons_all[22,:] = np.nan
lons_all[23,:] = np.nan
lats_all[22,:] = np.nan
lats_all[23,:] = np.nan
#%%
# Trajectories

fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
gs = fig.add_gridspec(figure = fig, nrows = 2, ncols = 1, height_ratios=[5,1])

ax = fig.add_subplot(gs[r,c], projection=projection)
ax.coastlines(resolution='50m',zorder=3)
ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
ax.set_extent([179, -109, 10, 40])

ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap0, alpha = 0.8, vmin = 0, vmax = 90, s = 1500,zorder=0,transform = cartopy.crs.PlateCarree()) 
if col == 0:
    ax.scatter(lons_all,lats_all, marker='.', c = time_all, cmap = cmap, alpha = 0.8, vmin = 0, vmax = 90, s = 1,zorder=1,transform = cartopy.crs.PlateCarree())   
scat = ax.scatter(lons_traj,lats_traj, c = time_all, cmap = cmap, vmin = 0, vmax = 1000,  marker='.', s = 1, zorder = 2, transform = cartopy.crs.PlateCarree())
rect= plt.Rectangle((lonmin-2,latmin-2), 12, 12, ec = 'lightblue', lw = 6, fc = 'none', zorder = 3, transform = cartopy.crs.PlateCarree())
ax.add_patch(rect)
#ax.scatter(lon_sink, lat_sink, marker= '.', c = Ts_3yr, cmap = cmap, linewidth = 3, vmin = 0, vmax = 1000, s = 300,zorder=4, transform = cartopy.crs.PlateCarree()) #"$\u25A1$"
ax.title.set_text('radius = 1mm')

cbaxes = fig.add_axes([0.1, 0.3, 0.8, 0.01])
cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='neither')
cbar.set_label(label='Age [days]', size=18)

plt.rc('font', size = 18)
gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=False, linewidth=0.5,
                      color='gray', alpha=0.5, linestyle='--')
gl.xlocator = mticker.FixedLocator(np.arange(-180, -110, 30))
gl.ylocator = mticker.FixedLocator(np.arange(10, 640, 10))
gl.xlabels_bottom = True
gl.xformatter = LONGITUDE_FORMATTER
gl.ylabels_left = True
gl.yformatter = LATITUDE_FORMATTER

#%%
# Particles that sink
fig2 = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
gs2 = fig2.add_gridspec(figure = fig2, nrows = 2, ncols = 1, height_ratios=[5,1])
ax = fig2.add_subplot(gs2[r,c], projection=projection)
ax.coastlines(resolution='50m',zorder=3)
ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
ax.set_extent([179, -109, 10, 40])

ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap0, alpha = 0.5, vmin = 0, vmax = 90, s = 1500,zorder=0,transform = cartopy.crs.PlateCarree()) 
# if col == 0:
#     ax.scatter(lons_all,lats_all, marker='.', c = time_all, cmap = cmap, alpha = 0.8, vmin = 0, vmax = 90, s = 0.1,zorder=1,transform = cartopy.crs.PlateCarree())   
scat = ax.scatter(lons_traj,lats_traj, c = time_all, cmap = cmap, vmin = 0, vmax = 1000,  marker='.', s = 50, zorder = 2, transform = cartopy.crs.PlateCarree())
rect= plt.Rectangle((lonmin-2,latmin-2), 12, 12, ec = 'lightblue', lw = 6, fc = 'none', zorder = 3, transform = cartopy.crs.PlateCarree())
ax.add_patch(rect)
ax.scatter(lon_sink, lat_sink, marker= '.', c = Ts_3yr,  cmap = cmap, linewidth = 5, vmin = 0, vmax = 1000, s = 2000,zorder=4, transform = cartopy.crs.PlateCarree()) #"$\u25A1$"
ax.title.set_text('radius = 1mm')

cbaxes = fig2.add_axes([0.1, 0.3, 0.8, 0.01])
cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='neither')
cbar.set_label(label='Age [days]', size=18)

plt.rc('font', size = 18)
gl2 = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=False, linewidth=0.5,
                      color='gray', alpha=0.5, linestyle='--')
gl2.xlocator = mticker.FixedLocator(np.arange(-180, -110, 30))
gl2.ylocator = mticker.FixedLocator(np.arange(10, 640, 10))
gl2.xlabels_bottom = True
gl2.xformatter = LONGITUDE_FORMATTER
gl2.ylabels_left = True
gl2.yformatter = LATITUDE_FORMATTER