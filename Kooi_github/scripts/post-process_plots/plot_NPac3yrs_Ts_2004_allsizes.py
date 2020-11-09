#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 09:46:24 2020

@author: Lobel001
"""

# Adapted from plot_Ts_Vs_2004_allsizes.py but only for particles released in N Pacific and run for 3 yrs

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np 
import cartopy
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
#import math
from numpy import *
#import cmocean
import os #, fnmatch 
#import cartopy
import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")


rho =  '920' #np.array([920])  # [kgm-3]: density of the plastic 
#res = '2x2' # [deg]: resolution of the global release of particles
#loc = 'global'


fname = '3yr_NPac_3D_grid2x2_allrho_allr_1089days_60dtsecs_12hrsoutdt.nc'
dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/allrho/res_2x2/allr/'
dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/'

plt.close("all")

#%%

yr = 2004
numyrs = 1
#seas = ['DJF', 'MAM', 'JJA', 'SON']

cmap = plt.cm.get_cmap('magma_r', 10) #RdBu #Blues #vidiris
# row = 1 
# col = 1

''' prepare cartopy projection and gridspec'''

projection = cartopy.crs.PlateCarree(central_longitude=180) #central_longitude=72+180)
plt.rc('font', size = 30)

    
# size_st = 3 # this is the superscript for 1e-0+'size1', hence starts on 1e-02
# numsizes = 5
sizes = np.arange(3,8)
fig_w = 25 #8 #30# 10 #30 # 30 25 #28
fig_h = 20 #7 #25#9 #25 #30 30 #18
fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
gs = fig.add_gridspec(figure = fig, nrows = 6, ncols = 2, height_ratios=[6,6,6,6,6,1],hspace = 0.5) #, wspace = 0.5) #, height_ratios=[10,1,10,1], wspace = 0.05, hspace = 1) # gridspec.GridSpec
 
# fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
# gs2 = fig.add_gridspec(figure = fig, nrows = 6, ncols = 1, height_ratios=[6,6,6,6,6,1],hspace = 0.5)   
#%%  
for row in range(5): #for ids in range(size_st,size_st+numsizes):#8): #sizes
    for col in range(2): # params
        i = row*2 + col
        si = sizes[row]
        size = 'r1e-0'+str(si) #size = 'r1e-0'+str(ids) #'r1e-02'
    
        print(size)
    
        if not os.path.isfile(dirread+fname):
            print('%s not found' %fname)
        else:          
            data = Dataset(dirread+fname,'r')  
            time = data.variables['time'][0,:]/86400
            time_all = (data.variables['time'][:,:]/86400)#-data.variables['time'][0,0]
            rho2 = float(rho) #rho_all[ii]
            size2 = float(size[1:len(size)]) #size_all[ii]
            rho_output=data.variables['rho_pl'][:]
            r_output=data.variables['r_pl'][:]
           
            inds = np.where((rho_output==rho2) & (r_output.data==size2))[0]
           
            lons=data.variables['lon'][inds] #[:]
            lats=data.variables['lat'][inds] 
            depths =data.variables['z'][inds]
            
            ''' find the first time and max depth where particles sink (below 1 m, since that's the initial release depth) '''      
            time = time - time[0]
            Ts = np.zeros(depths.shape[0])  #t_set = 
            Ts[:] = np.nan 
    
            z = []   
            depths2 = []
            lon_sink = np.zeros(depths.shape[0])
            lat_sink = np.zeros(depths.shape[0])
            for i in range(depths.shape[0]): #9620: number of particles 
                depths2 = depths[i,:] 
                #max_lon[i] = np.max(lons[i,:])
                #min_lon[i] = np.min(lons[i,:])
                z = np.where(depths2 > 1.)
                
                z2 = [0,] * depths.shape[1]
                zind = []
                z_set = []
                if z[0].any():
                    z_set = z[0][0]-1
        
                    # for ii in range(z_set,depths.shape[1]): #(2,depths.shape[1]):
                    #     z2[ii-1] = depths[i,ii]-depths[i,ii-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
                    # zf = np.where(np.array(z2) < 0.)[:]           
                    # zind = zf[0][0] if np.array(zf).any() else [] #np.nan
                
                    
                    Ts[i] = time[z_set] #t_set = 
                    lon_sink[i] = lons[i,z_set]
                    lat_sink[i] = lats[i,z_set]
                    #Ts[Ts==0]= np.nan
        # print('For '+size+' min lon = '+str(np.min(min_lon)))
        # print('For '+size+' max lon = '+str(np.max(max_lon)))
        #     t_set = np.array(t_set)
        #     t_set[t_set==0]= np.nan
    
           
        
        #     t_set_all[:,idy] = t_set   
        #     #vs_max_all[:,idy] = vs_max
            
        
        # t_set_all[t_set_all==0] = np.nan        
        # Ts = np.nanmean(t_set_all,axis=1) #[:,idx]
        # #Vs = np.nanmean(vs_max_all,axis=1)
    
        Ts[isnan(Ts)] = 1500. # to that     
    #%%    
        c = 1
        r = int(si)-3
    
        if si == 3:
            size_name = '1 mm'
        if si == 4:
            size_name = '0.1 mm'
        if si == 5:
            size_name = '10 \u03bcm'
        if si == 6:
            size_name = '1 \u03bcm'
        if si == 7:
            size_name = '0.1 \u03bcm'
        
        letter = (chr(ord('a') + row))
        # if ids == 2:
        #     letter = '(a)'
        #     r = 0
        #     c = 0
        #     Ts02 = Ts
        #     size_name = '10 mm'
        # if ids == 3:
        #     letter = '(b)'
        #     r = 0
        #     c = 1
        #     Ts03 = Ts
        #     size_name = '1 mm'
        # if ids == 4:
        #     letter = '(c)'
        #     r = 1
        #     c = 0
        #     Ts04 = Ts
        #     size_name = '0.1 mm'
        # if ids == 5:
        #     letter = '(d)'
        #     r = 1
        #     c = 1
        #     Ts05 = Ts
        #     size_name = '10 \u03bcm'
        # if ids == 6:
        #     letter = '(e)'
        #     r = 2
        #     c = 0
        #     Ts06 = Ts
        #     size_name = '1 \u03bcm'
        # if ids == 7:
        #     letter = '(f)'
        #     r = 2
        #     c = 1
        #     Ts07 = Ts
        #     size_name = '0.1 \u03bcm'
        
        if col == 0:
            plt.rc('font', size = 30)
            ax = fig.add_subplot(gs[row,col], projection=projection)
            ax.coastlines(resolution='50m',zorder=3)
            ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
            #ax.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())
            ax.set_extent([-150, -138, 18, 30]) #ax.set_extent([-180, -80, 0, 40])
            
            scat = ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap, vmin = 0, vmax = 1000, s = 1000,zorder=1,transform = cartopy.crs.PlateCarree()) 
            ax.title.set_text(letter+ ' radius = '+size_name)  #plt.title(letter+ ' radius = '+size_name)#, fontsize = 26) #str(yr)+' all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)
                
            cbaxes = fig.add_axes([0.1, 0.06, 0.8, 0.01])
            cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max')
            cbar.set_label(label='[days]') #, size=18)
            plt.rc('font', size = 20)
            gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=False, linewidth=0.5,
                                  color='gray', alpha=0.5, linestyle='--')
            gl.xlocator = mticker.FixedLocator(np.arange(-152, -132, 4))
            gl.ylocator = mticker.FixedLocator(np.arange(18, 32, 2))
            gl.xlabels_bottom = True
            gl.xformatter = LONGITUDE_FORMATTER
            gl.ylabels_left = True
            gl.yformatter = LATITUDE_FORMATTER
        
        # test = lons[np.where(lons==np.nan)]
        # time_all[test]==np.nan
        if col == 1:
            plt.rc('font', size = 30)
            if si == 2: #ids ==2:
                fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
                gs = fig.add_gridspec(figure = fig, nrows = 4, ncols = 2, height_ratios=[6,6,6,1],hspace = 0.5) 
            ax = fig.add_subplot(gs[row,col], projection=projection)
            ax.coastlines(resolution='50m',zorder=3)
            ax.add_feature(cartopy.feature.LAND, color='lightgrey') #, zorder=2)
            ax.set_extent([110, -80, 0, 50]) #,crs=cartopy.crs.PlateCarree())
            ax.set_ylim([0, 50])
            scat1 = ax.scatter(lons[:,0], lats[:,0], c = 'k', marker='.', vmin = 0, vmax = 1000, s = 50,zorder=2, transform = cartopy.crs.PlateCarree())
            scat2 = ax.scatter(lons, lats, marker='.',  vmin = 0, vmax = 1000, s = 50,zorder=1, transform = cartopy.crs.PlateCarree()) # c = time_all, cmap = cmap,
            scat3 = ax.scatter(lon_sink, lat_sink, marker='.',  c=Ts,cmap = cmap, vmin = 0, vmax = 1000, s = 500,zorder=1, transform = cartopy.crs.PlateCarree())
            ax.title.set_text(letter+ ' radius = '+size_name) #plt.title
            
            cbaxes2 = fig.add_axes([0.1, 0.06, 0.8, 0.01])
            cbar2 = plt.colorbar(scat, cax=cbaxes2, orientation="horizontal", aspect=50, extend='max')
            cbar2.set_label(label='[days]') #, size=18)
            plt.rc('font', size = 20)
            gl2 = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=False, linewidth=0.5,
                                  color='gray', alpha=0.5, linestyle='--')
            # gl2.xlocator = mticker.FixedLocator(np.arange(110, -80, 4))
            # gl2.ylocator = mticker.FixedLocator(np.arange(0, 50, 2))
            gl2.xlabels_bottom = False #True
            gl2.xformatter = LONGITUDE_FORMATTER
            gl2.ylabels_left = False #True
            gl2.yformatter = LATITUDE_FORMATTER
        
        
        #cbar.ax.tick_params(labelsize=18)
        #plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max', label='[days]')
                  
    
        #plt.savefig(dirwritefigs + 'Ts.pdf')
    
    # with open(dirwrite+'Ts2004_allsizes_for_diffplot_with_advection.pickle', 'wb') as f:
    #     pickle.dump([lons,lats,Ts02,Ts03,Ts04,Ts05,Ts06,Ts07], f)
    
    # ''' plotting Ts against Vs which could be linear relationship''' - but it isn't! 
    # fig3 = plt.figure(figsize =(20,10))
    # plt.scatter(Ts,Vs) 
    # plt.xlabel('Ts')
    # plt.ylabel('Vs')
    #plt.ylim(bottom=0, top =0.05) 

