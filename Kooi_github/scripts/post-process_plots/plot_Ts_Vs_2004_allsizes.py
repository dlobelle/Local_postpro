#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 07:36:23 2020

@author: Lobel001
"""

# Sinking velocities for 920 kgm-3 data, all sizes (1e-2 to 1e-7) and only for 2004
    # output plots: global maps of settling velocities only here (the year defined as "most typical") 

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np 
import cartopy
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# from cartopy.util import add_cyclic_point
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

#import math
#from pylab import *
import cmocean
import os #, fnmatch 
#import cartopy
import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")

rho =  '920' #np.array([920])  # [kgm-3]: density of the plastic 
# size = 'r1e-02' #np.array([1e-7]) # [m]: size of plastic
res = '2x2' # [deg]: resolution of the global release of particles
loc = 'global'

plt.close("all")

#%%
#list_aux = os.listdir(dirread) # complete list of all files in dirread

yr = 2004
numyrs = 1
seas = ['DJF', 'MAM', 'JJA', 'SON']

cmap = plt.cm.get_cmap('magma_r', 11) #RdBu #Blues #vidiris
row = 1 
col = 1

projection = cartopy.crs.PlateCarree() #PlateCarree() #central_longitude=0)
fig = plt.figure(figsize=(30,17), constrained_layout=True) #True) # 10,5
gs = fig.add_gridspec(figure = fig, nrows = 4, ncols = 2) #, height_ratios=[10,1,10,1], wspace = 0.05, hspace = 1) # gridspec.GridSpec

for ids in range(2,8): #sizes
    size = 'r1e-0'+str(ids) #'r1e-02'
    if size == 'r1e-02':
        size_fn = 'r0.01'
    else:
        size_fn = size #'r1e-02'
        
    dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
    dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    
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
                        lons = data.variables['lon'][:]
                        lats = data.variables['lat'][:]
                        depths =data.variables['z'][:]
                        vs = data.variables['vs'][:]
                    else:
                        lons = np.vstack((lons,data.variables['lon'][:])) 
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
               
                lons=data.variables['lon'][inds] #[:]
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
                lons=data.variables['lon'][:] 
                lats=data.variables['lat'][:] 
                depths =data.variables['z'][:]
                vs = data.variables['vs'][:]
        
            ''' find the first time and max depth where particles sink (below 1 m, since that's the initial release depth) '''      
        time = time - time[0]
        #z_set = [0, ] * depths.shape[0]
        t_set = [0, ] * depths.shape[0]
        vs_max = [np.nan, ] * depths.shape[0]
        #vs_max[:] = np.nan
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
                    #if depths2[ii]>0.: #1.
                    z2[ii-1] = depths[i,ii]-depths[i,ii-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
                zf = np.where(np.array(z2) < 0.)[:]           
                zind = zf[0][0] if np.array(zf).any() else [] #np.nan
            
                ''' finding the largest (positive) vs to plot max sinking velocity''' 
                if np.array(zind).any():
                    v0 = vs[i,z_set:zind]
                    v_f = np.where(v0>0.)
                    vs_max[i] = np.max(v0[v_f[0]]) if np.array(v_f).any() else np.nan #vs[i,z_set-1]
                else: 
                    vs_max[i] = np.nan 
                # vs_max[i] = np.max(vs[i,:]) #[i,0:zind]) if np.array(zind).any() else np.nan
                t_set[i] = time[z_set]
        
        
        if size == 0.01:
            sz_fn = 'r1e-02'
        else:
            sz_fn = size
       
        # fig = plt.figure(figsize=(20,10))
    
        # cmap = plt.cm.get_cmap('magma_r', 11) #RdBu #Blues #vidiris
        
        # m = Basemap(projection='robin',lon_0=-180,resolution='l') #, ax=axs[idx])
        # m.drawparallels(np.array([-60,-30,0,30,60]), labels=[True, False, False, True], linewidth=2.0, size=20,zorder=0)
        # m.drawmeridians(np.array([50,150,250,350]), labels=[False, False, False, True], linewidth=2.0, size=20,zorder=0)
        # m.drawcoastlines()
        # m.fillcontinents(color='lightgrey',zorder=3)
        # xs, ys = m(lons[:,0], lats[:,0])
        # scat = m.scatter(xs, ys, marker='.', c=vs_max,cmap = cmap, s = 100,zorder=1) #vmin = 0, vmax = 1e-6, 
        # cbar = m.colorbar(label = 'm/s') #,ax=axs[idx])  
        # cbar.set_label(label='[m/s]', size='20')   
        # cbar.ax.tick_params(labelsize=20) 
        # plt.title(str(yr)+' '+str(s)+' average maximum settling velocity [m/s], rho ='+str(rho)+', size ='+str(size) ,size = 20)
           
        t_set_all[:,idy] = t_set   
        vs_max_all[:,idy] = vs_max
        
    
    t_set_all[t_set_all==0] = np.nan        
    Ts = np.nanmean(t_set_all,axis=1) #[:,idx]
    Vs = np.nanmean(vs_max_all,axis=1)
    
    # with open(dirwrite+str(rho)+str(sz_fn)+'JJA2001_av_Ts_Zmax.pickle', 'wb') as f:
    #     pickle.dump([lons,lats,lon_set,lat_set,t_set,z_max,z_max2,i_noset,prob_set], f)

#%%    
    
#fig1 = plt.figure(figsize=(20,10))
    if ids == 2:
        letter = '(a)'
        r = 0
        c = 0
    if ids == 3:
        letter = '(b)'
        r = 0
        c = 1
    if ids == 4:
        letter = '(c)'
        r = 1
        c = 0
    if ids == 5:
        letter = '(d)'
        r = 1
        c = 1
    if ids == 6:
        letter = '(e)'
        r = 2
        c = 0
    if ids == 7:
        letter = '(f)'
        r = 2
        c = 1
    
    
    ax = fig.add_subplot(gs[r,c], projection=projection)
    ax.coastlines(resolution='50m',zorder=3)
    ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
    ax.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())
    #ax.set_extent([73, 72, -70, 60])
    
    scat = ax.scatter(lons[:,0], lats[:,0], marker='.', c=Ts,cmap = cmap, vmin = 0, vmax = 90, s = 100,zorder=1) #scat = 
    # cbar = plt.colorbar(scat, label = 'days') #,ax=axs[idx])  
    # cbar.set_label(label='days', size='20')   
    # cbar.ax.tick_params(labelsize=20) 
    plt.title(letter+ ' radius = '+str(size)+'m', fontsize = 26) #str(yr)+' all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)
    

#cbar = plt.colorbar(scat, label = 'days') #,ax=axs[idx])
cbaxes = fig.add_axes([0.1, 0.15, 0.8, 0.01])
plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='both', label='[m]')

plt.savefig('/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/Ts.pdf')

#fig1 = plt.figure(figsize=(20,10))
#m = Basemap(projection='robin',lon_0=-180,resolution='l') #, ax=axs[idx])
# m.drawparallels(np.array([-60,-30,0,30,60]), labels=[True, False, False, True], linewidth=2.0, size=20,zorder=0)
# m.drawmeridians(np.array([50,150,250,350]), labels=[False, False, False, True], linewidth=2.0, size=20,zorder=0)
# m.drawcoastlines()
# m.fillcontinents(color='lightgrey',zorder=3)
# xs, ys = m(lons[:,0], lats[:,0])
# scat = m.scatter(xs, ys, marker='.', c=Ts,cmap = cmap, vmin = 0, vmax = 90, s = 100,zorder=1) #scat = 
# cbar = m.colorbar(label = 'days') #,ax=axs[idx])  
# cbar.set_label(label='days', size='20')   
# cbar.ax.tick_params(labelsize=20) 
# plt.title(str(yr)+' all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)


# fig2 = plt.figure(figsize=(20,10))

# cmap = plt.cm.get_cmap('magma_r', 11) #RdBu #Blues #vidiris

# m = Basemap(projection='robin',lon_0=-180,resolution='l') #, ax=axs[idx])
# m.drawparallels(np.array([-60,-30,0,30,60]), labels=[True, False, False, True], linewidth=2.0, size=20,zorder=0)
# m.drawmeridians(np.array([50,150,250,350]), labels=[False, False, False, True], linewidth=2.0, size=20,zorder=0)
# m.drawcoastlines()
# m.fillcontinents(color='lightgrey',zorder=3)
# xs, ys = m(lons[:,0], lats[:,0])
# scat = m.scatter(xs, ys, marker='.', c=Vs,cmap = cmap, s = 100,zorder=1) #vmin = 0, vmax = 1e-6,
# cbar = m.colorbar(label = 'm/s') #,ax=axs[idx])  
# cbar.set_label(label='[m/s]', size='20')   
# cbar.ax.tick_params(labelsize=20) 
# plt.title(str(yr)+' all seas average maximum settling velocity [m/s], rho ='+str(rho)+', size ='+str(size) ,size = 20)


# ''' plotting Ts against Vs which could be linear relationship''' 
# fig3 = plt.figure(figsize =(20,10))
# plt.scatter(Ts,Vs) 
# plt.xlabel('Ts')
# plt.ylabel('Vs')
#plt.ylim(bottom=0, top =0.05) 

#%%
# v900 = vs[1000,:]
# d900 = depths[1000,:]

# fig4 = plt.figure(figsize =(20,10))
# plt.plot(time[0:6],v900[0:6])
# plt.ylim(bottom=-0.1, top =0) 
# plt.title('settling velocity for 1 trajectory rho ='+str(rho)+', size ='+str(size))

# fig5 = plt.figure(figsize =(20,10))

# plt.plot(time[0:6],d900[0:6]*-1)
# plt.title('depth for 1 trajectory rho ='+str(rho)+', size ='+str(size))