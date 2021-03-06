#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 07:36:23 2020

@author: Lobel001
"""

# Sinking velocities for 920 kgm-3 data, all sizes (1e-2 to 1e-7) and only for 2004
    # output plots: global maps of sinking timescales and velocities only here (the year defined as "most typical") 

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np 
import cartopy
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
from numpy import *
#import math
#from pylab import *
import cmocean.cm as cmo
import os #, fnmatch 
#import cartopy
import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")

''' Choose what to plot '''
plots = 'allsizes' #'allsizes' # '1size_allseas'


rho =  '920' #np.array([920])  # [kgm-3]: density of the plastic 
# size = 'r1e-02' #np.array([1e-7]) # [m]: size of plastic
res = '2x2' # [deg]: resolution of the global release of particles
loc = 'global'


plt.close("all")

#%%

yr = 2004
numyrs = 1
seas = ['DJF', 'MAM', 'JJA', 'SON']

cmap = plt.cm.get_cmap('magma_r', 9) #RdBu #Blues #vidiris cmap = plt.cm.get_cmap(cmo.matter_r,9) #
row = 1 
col = 1

''' prepare cartopy projection and gridspec'''

projection = cartopy.crs.PlateCarree() #central_longitude=72+180)
plt.rc('font', size = 30)

if plots == 'allsizes':
    size_st = 2 # this is the superscript for 1e-0+'size1', hence starts on 1e-02
    numsizes = 6
    fig_w = 24 #8 #30# 10 #30 # 30 25 #28
    fig_h = 21 #7 #25#9 #25 #30 30 #18
    fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) # 10,5
    gs = fig.add_gridspec(figure = fig, nrows = 4, ncols = 2, height_ratios=[7,7,7,1]) #,hspace = 0.5, wspace = 0.5) #, height_ratios=[10,1,10,1], wspace = 0.05, hspace = 1) # gridspec.GridSpec
elif plots == '1size_allseas':
    #plt.rc('font', size = 20)
    size_st = 4  # this is the superscript for 1e-0+'size1', hence it's only for 1e-02
    numsizes = 1
    fig_w = 15 #30 
    fig_h = 8 #17
    fig_1s = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) #True) # 10,5 #
    gs_1s = fig_1s.add_gridspec(figure = fig_1s, nrows = 3, ncols = 2, height_ratios=[6.5,6.5,1], hspace = 0.5)
    
     
for ids in range(size_st,size_st+numsizes):#8): #sizes
    size = 'r1e-0'+str(ids) #'r1e-02'
    if size == 'r1e-02':
        size_fn = 'r0.01'
    else:
        size_fn = size #'r1e-02'
        
    dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
    dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/'
    
    if size == 'r1e-05' or size == 'r1e-07':
        num_part = 9582
    else:
        num_part = 9620
        
    t_set_all = np.zeros((num_part,len(seas))) # just keeping it to have the same num of variables in the pickle. 
    vs_max_all = np.zeros((num_part,len(seas)))
    
    
    for idy,s in enumerate(seas):
    
        print(size, s, yr)
        
        if size == 'r1e-05' or size == 'r1e-07': # for those sizes that have 3 latitudinal bands 
            loc = ['south_global','eq_global','north_global']
            for idz, lo in enumerate(loc):
                fname = '%s_%s_%s_3D_grid2x2_rho%s_%s_90days_60dtsecs_12hrsoutdt.nc' % (lo, s, yr, rho, size)
                if not os.path.isfile(dirread+fname):
                    print('%s not found' %fname)
                else:          
                    data = Dataset(dirread+fname,'r') 
                    time = data.variables['time'][0,:]/86400
                    if idz == 0:
                        lons = data.variables['lon'][:]#+72+180
                        lats = data.variables['lat'][:]
                        depths =data.variables['z'][:]
                        vs = data.variables['vs'][:]
                    else:
                        lons = np.vstack((lons,data.variables['lon'][:]))#+72+180 
                        lats= np.vstack((lats,data.variables['lat'][:]))
                        depths =np.vstack((depths,data.variables['z'][:]))
                        vs = np.vstack((vs,data.variables['vs'][:]))
        elif size == 'r1e-03' or size == 'r1e-06':
            fname = 'global_%s_%s_3D_grid2x2_allrho_allr_90days_60dtsecs_12hrsoutdt.nc' % (s, yr)
            dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/allrho/res_'+res+'/allr/'
            if not os.path.isfile(dirread+fname):
                print('%s not found' %fname)
            else:          
                data = Dataset(dirread+fname,'r')  
                time = data.variables['time'][0,:]/86400
                rho2 = float(rho) #rho_all[ii]
                size2 = float(size[1:len(size)]) #size_all[ii]
                rho_output=data.variables['rho_pl'][:]
                r_output=data.variables['r_pl'][:]
               
                inds = np.where((rho_output==rho2) & (r_output.data==size2))[0]
               
                lons=data.variables['lon'][inds]#+72+180
                lats=data.variables['lat'][inds] 
                depths =data.variables['z'][inds]
                vs = data.variables['vs'][inds]
                
        else:     
            fname = 'global_%s_%s_3D_grid2x2_rho%s_%s_90days_60dtsecs_12hrsoutdt.nc' % (s, yr, rho, size_fn)
            if not os.path.isfile(dirread+fname):
                print('%s not found' %fname)
            else:          
                data = Dataset(dirread+fname,'r') 
                time = data.variables['time'][0,:]/86400
                lons=data.variables['lon'][:]#+72+180
                lats=data.variables['lat'][:] 
                depths =data.variables['z'][:]
                vs = data.variables['vs'][:]
        
        ''' find the first time and max depth where particles sink (below 1 m, since that's the initial release depth) '''      
        time = time - time[0]
        #t_set = [0, ] * depths.shape[0]
        t_set = np.zeros(depths.shape[0]) 
        t_set[:] = np.nan 
        #vs_max = [np.nan, ] * depths.shape[0]

        z = []    
        for i in range(depths.shape[0]): #9620: number of particles 
            depths2 = depths[i,:] 
            z = np.where(depths2 > 1.)
            
            z2 = [0,] * depths.shape[1]
            zind = []
            z_set = []
            if z[0].any():
                z_set = z[0][0]-1
    
                for ii in range(z_set,depths.shape[1]): #(2,depths.shape[1]):
                    z2[ii-1] = depths[i,ii]-depths[i,ii-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
                zf = np.where(np.array(z2) < 0.)[:]           
                zind = zf[0][0] if np.array(zf).any() else [] #np.nan
            
                ''' finding the largest (positive) vs_max to plot max sinking velocity''' 
                # if np.array(zind).any():
                #     v0 = vs[i,z_set:zind]
                #     v_f = np.where(v0>0.)
                #     vs_max[i] = np.max(v0[v_f[0]]) if np.array(v_f).any() else np.nan #vs[i,z_set-1]
                # else: 
                #     vs_max[i] = np.nan  #np.max(vs[i,:]) #[i,0:zind]) if np.array(zind).any() else np.nan
                
                t_set[i] = time[z_set]

        t_set = np.array(t_set)
        #t_set[t_set==0]= np.nan

        if plots == '1size_allseas':
            t_set[isnan(t_set)] = 100.
            if idy == 0:
                letter = '(a)'
                r_1s = 0
                c_1s = 0
            if idy == 1:
                letter = '(b)'
                r_1s = 0
                c_1s = 1
            if idy == 2:
                letter = '(c)'
                r_1s = 1
                c_1s = 0
            if idy == 3:
                letter = '(d)'
                r_1s = 1
                c_1s = 1

            #plt.rc('font', size = 20)
            ax_1s = fig_1s.add_subplot(gs_1s[r_1s,c_1s], projection=projection)
            ax_1s.coastlines(resolution='50m', zorder=3)
            ax_1s.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
            ax_1s.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())
            #ax_1s.set_extent([-28, -20, 45, 55]) 
            
            scat = ax_1s.scatter(lons[:,0], lats[:,0], marker='.', c=t_set,cmap = cmap, vmin = 0, vmax = 90, s = 20, zorder=1) #,crs=cartopy.crs.PlateCarree()) #scat = 
            # cbar = plt.colorbar(scat, label = 'days') #,ax=axs[idx])  
            # cbar.set_label(label='days', size='20')   
            # cbar.ax.tick_params(labelsize=20) 
            plt.title(letter+ ' '+str(s), fontsize = 20)
    
        t_set_all[:,idy] = t_set   
        #vs_max_all[:,idy] = vs_max
     
    
    #t_set_all[t_set_all==0] = np.nan        
    Ts = np.nanmean(t_set_all,axis=1) #[:,idx]
    #Vs = np.nanmean(vs_max_all,axis=1)
    Ts[isnan(Ts)] = 100.

#%%    
    if plots == 'allsizes':    
        
        #fig1 = plt.figure(figsize=(20,10))
        if ids == 2:
            letter = '(a)'
            r = 0
            c = 0
            Ts02 = Ts
            size_name = '10 mm'
        if ids == 3:
            letter = '(b)'
            r = 0
            c = 1
            Ts03 = Ts
            size_name = '1 mm'
        if ids == 4:
            letter = '(c)'
            r = 1
            c = 0
            Ts04 = Ts
            size_name = '0.1 mm'
        if ids == 5:
            letter = '(d)'
            r = 1
            c = 1
            Ts05 = Ts
            size_name = '10 \u03bcm'
        if ids == 6:
            letter = '(e)'
            r = 2
            c = 0
            Ts06 = Ts
            size_name = '1 \u03bcm'
        if ids == 7:
            letter = '(f)'
            r = 2
            c = 1
            Ts07 = Ts
            size_name = '0.1 \u03bcm'
        
        
        ax = fig.add_subplot(gs[r,c], projection=projection)
        ax.coastlines(resolution='50m',zorder=3)
        ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
        ax.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())
        #ax.set_extent([73, 72, -70, 60])
        #ax.set_extent([-150, -140, 45, 55])

        
        scat = ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap, vmin = 0, vmax = 90, s = 50,zorder=1) #scat = 
        # cbar = plt.colorbar(scat, label = 'days') #,ax=axs[idx])  
        # cbar.set_label(label='days', size='20')   
        # cbar.ax.tick_params(labelsize=20) 
        plt.title(letter+ ' radius = '+size_name)#, fontsize = 26) #str(yr)+' all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)
            
        cbaxes = fig.add_axes([0.1, 0.06, 0.8, 0.01])
        cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max')
        cbar.set_label(label='[days]') #, size=18)
        #cbar.ax.tick_params(labelsize=18)
        #plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max', label='[days]')
              

if plots == '1size_allseas':  
    plt.rc('font', size = 18)
    #cbar = plt.colorbar(scat, label = 'days') 
    cbaxes = fig_1s.add_axes([0.1, 0.085, 0.8, 0.01])
    cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max')
    cbar.set_label(label='[days]') #, size=18)
    #plt.rc('font', size = 18)

#%% UNCOMMENT TO SAVE FIGS OR VARIABLES  - NB now I've converted NaN (no sinking within 90 days) to 100 days.       
    # with open(dirwrite+'Ts_2004_4seasons'+size+'.pickle', 'wb') as f:
    #     pickle.dump([lons,lats,Ts02,Ts03,Ts04,Ts05,Ts06,Ts07], f)

    #plt.savefig(dirwritefigs + 'Ts_2004_4seasons'+size+'.pdf')
#else: 
#    plt.savefig(dirwritefigs + 'Ts.jpg') #'pdf')

# with open(dirwrite+'Ts2004_allsizes_for_diffplot_with_advection.pickle', 'wb') as f:
#     pickle.dump([lons,lats,Ts02,Ts03,Ts04,Ts05,Ts06,Ts07], f)

# ''' plotting Ts against Vs which could be linear relationship''' - but it isn't! 
# fig3 = plt.figure(figsize =(20,10))
# plt.scatter(Ts,Vs) 
# plt.xlabel('Ts')
# plt.ylabel('Vs')
#plt.ylim(bottom=0, top =0.05) 

