#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 01:36:59 2020

@author: Lobel001
"""

# Based on Erik's script: plotting_NEMO.ipynb

import numpy as np
import xarray as xr
from netCDF4 import Dataset
#import pandas as pd
import matplotlib.pyplot as plt
import cartopy
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from cartopy.util import add_cyclic_point
#import matplotlib.colors as colors
#import matplotlib.gridspec as gridspec
#import matplotlib.patches as mpatches
#import matplotlib.path as mpath
#import matplotlib.ticker as mticker
#from matplotlib.patches import Polygon, Ellipse
#from glob import glob
import cmocean.cm as cmo
import os
from numpy import *
#import gsw
#import pickle

# Added this (from a forum) as a temporary fix to error I was getting regarding 'GeoAxes not having a _hold function'
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

#CHOOSE SIZE TO PLOT TS
size = 'r1e-04'
rho = '920'

dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/NEMO_phys_params/'
dirread_Ts = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/allrho/res_2x2/allr/'
## To plot the MLD and algal conc (converted from MEDUSA phytoplankton conc) subplots in one figure: 2 columns 

#from time import sleep
#from progressbar import progressbar

''' Preparing the projection and subplots'''
plt.rc('font', size = 11)
projection = cartopy.crs.PlateCarree(central_longitude=72+180) #PlateCarree() #central_longitude=0)
fig_w = 10 
fig_h = 11 
fig = plt.figure(figsize=(fig_w,fig_h)) #10,10.5)) #, constrained_layout=True) #True) # 10,5 #(18,25)
gs = fig.add_gridspec(figure = fig, nrows = 5, ncols = 3, height_ratios=[6,6,6,6,1])#,hspace = 1, wspace = 1)

''' Extract all files within a season into one variable '''

seasnames = ['DJF', 'MAM', 'JJA', 'SON']
monvals = ['12','01','02','03','04','05','06','07','08','09','10','11']

    
    
for row in range(4): #progressbar(range(4)): # seasons
    #sleep(0.02)
    for col in range(3): # params
        i = row*3 + col
        seas = seasnames[row]
        

        # Lat, Lon = mask.variables['nav_lat'], mask.variables['nav_lon']
        # latvals = Lat[:]; lonvals = Lon[:] # extract lat/lon values to numpy arrays
        
        
        # with open(dirread+'MLD_2004_'+seas+'.pickle','rb') as f:  
        #     M0 = pickle.load(f)
        yr0 = '2004'
        yr1 = '2004'

        ai = fig.add_subplot(gs[row, col], projection=projection)
        ai.coastlines(resolution='50m',zorder=3)
        ai.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
        ai.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())

        if row == 0: 
            yr0 = '2003'
            mons = monvals[0:3]
        if row == 1: 
            mons = monvals[3:6]
        if row == 2: 
            mons = monvals[6:9]
        if row == 3: 
            mons = monvals[9:12]

        if i ==0 or i == 3 or i ==6 or i ==9:   #i == 2 or i ==4 or i ==6:
            
            fname = dirread+'MLD_2004_'+seas+'.nc'
            M0 = xr.open_dataset(fname).mldr10_1  
            M0 = M0.assign_coords(nav_lat=M0.nav_lat)
            M0 = M0.assign_coords(nav_lon=M0.nav_lon) 

            a = M0.plot(ax=ai, x='nav_lon', y='nav_lat', add_colorbar=False, vmin=0, vmax=200, rasterized=True, cmap=cmo.deep, zorder = 1, transform=cartopy.crs.PlateCarree()) 
            title = 'MLD in %s' % seas

                
        if i ==1 or i == 4 or i ==7 or i ==10:   #i == 3 or i ==5 or i ==7:
 
            fname = dirread+'algal_conc_2004_'+seas+'.nc'
            M = xr.open_dataset(fname).__xarray_dataarray_variable__  
            M = M.assign_coords(nav_lat=M.nav_lat)
            M = M.assign_coords(nav_lon=M.nav_lon) 

            a2 = M.plot(ax=ai, x='nav_lon', y='nav_lat', add_colorbar=False, vmin=0, vmax=4e7, rasterized=True, cmap=cmo.algae, zorder = 1, transform=cartopy.crs.PlateCarree()) 
            title = 'Algal conc. in %s' % seas
            

        #ai.set_title('%s) %s ' % (chr(ord('a') + row*2 + col), title))

        
        if i ==2 or i == 5 or i ==8 or i ==11:
            cmap = plt.cm.get_cmap('magma_r', 9) 
            fname = 'global_'+seas+'_2004_3D_grid2x2_allrho_allr_90days_30dtsecs_12hrsoutdt.nc' 
            
            if not os.path.isfile(dirread_Ts+fname):
                print('%s not found' %fname)
            else:          
                data = Dataset(dirread_Ts+fname,'r')  
                time = data.variables['time'][0,:]/86400
                rho2 = float(rho) #rho_all[ii]
                size2 = float(size[1:len(size)]) #size_all[ii]
                rho_output=data.variables['rho_pl'][:]
                r_output=data.variables['r_pl'][:]
               
                inds = np.where((rho_output==rho2) & (r_output.data==size2))[0]
               
                lons=data.variables['lon'][inds]
                lats=data.variables['lat'][inds] 
                depths =data.variables['z'][inds]
                vs = data.variables['vs'][inds]
                vs_init= data.variables['vs_init'][inds]
                w = data.variables['w'][inds]
            time = time - time[0]
            
            ''' Defining Ts using depth Vs + w > 0 m/s (particle's velocity + vertical advection is downward)'  '''
            vs_init2 = []
            w2 = []
            w_vs = []
            w_ind = []
            z_set2 = []
            t_set = np.zeros(depths.shape[0]) 
            t_set[:] = np.nan 
            boo = np.zeros((depths.shape[0],depths.shape[1]))
            
            for i in range(depths.shape[0]): #9620: number of particles 
                vs_init2 = vs_init[i,:]
                w2 = w[i,:]
                w_vs = vs_init2 + w2 
                w_ind = np.where(w_vs>0)
                
                z_set2 = []
                if w_ind[0].any():
                    z_set2 = w_ind[0][0]-1   
                    
                    t_set[i] = time[z_set2]
            
            t_set[isnan(t_set)] = 100.
            # if idy == 0:
            #     letter = '(a)'
            #     r_1s = 0
            #     c_1s = 0
            # if idy == 1:
            #     letter = '(b)'
            #     r_1s = 0
            #     c_1s = 1
            # if idy == 2:
            #     letter = '(c)'
            #     r_1s = 1
            #     c_1s = 0
            # if idy == 3:
            #     letter = '(d)'
            #     r_1s = 1
            #     c_1s = 1

            # ax = fig.add_subplot(gs_1s[r_1s,c_1s], projection=projection)
            # ax.coastlines(resolution='50m', zorder=3)
            # ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
            # ax.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())          
            a3 = ai.scatter(lons[:,0], lats[:,0], marker='.', c=t_set,cmap = cmap, vmin = 0, vmax = 90, s = 5, zorder=1,transform = cartopy.crs.PlateCarree()) #,crs=cartopy.crs.PlateCarree()) #scat = 

            title = '$T_s$ in %s' % seas
            

        ai.set_title('%s) %s ' % (chr(ord('a') + row*2 + col), title))

cbaxes = fig.add_axes([0.02, 0.045, 0.3, 0.015]) # defines the x, y, w, h of the colorbar 
plt.colorbar(a, cax=cbaxes, orientation="horizontal", aspect=100, extend='max', label='[m]', use_gridspec=True) #, fontsize = 22)
#plt.rc('font', size = 14)

cbaxes2 = fig.add_axes([0.35, 0.045, 0.3, 0.015]) # defines the x, y, w, h of the colorbar 
plt.colorbar(a2, cax=cbaxes2, orientation="horizontal", aspect=100, extend='max', label='[no. m$^{-3}$]', use_gridspec=True) #, fontsize = 22)

cbaxes3 = fig.add_axes([0.68, 0.045, 0.3, 0.015]) # defines the x, y, w, h of the colorbar 
plt.colorbar(a3, cax=cbaxes3, orientation="horizontal", aspect=100, extend='max', label='[days]', use_gridspec=True) #, fontsize = 22)
#plt.rc('font', size = 14)

fig.canvas.draw()
plt.tight_layout()
#plt.savefig('/home/dlobelle/Kooi_figures/non_parcels_output/NEMO_MLD_algae_allseas_2004.pdf')