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
size = 'r1e-06' #np.array([1e-7]) # [m]: size of plastic
res = '2x2' # [deg]: resolution of the global release of particles
loc = 'global'

if size == 'r1e-02':
    size_fn = 'r0.01'
else:
    size_fn = size #'r1e-02'

dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'


plt.close("all")

#%%
#list_aux = os.listdir(dirread) # complete list of all files in dirread

yr = 2004
numyrs = 1
seas = ['DJF', 'MAM', 'JJA', 'SON']

t_set_all = np.zeros((9620,len(seas))) # just keeping it to have the same num of variables in the pickle. 
vs_max_all = np.zeros((9620,len(seas)))

for idy,s in enumerate(seas):

    print(s, yr)
    
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
                    vs = np.vstack((depths,data.variables['vs'][:]))
    elif size == 'r1e-03' or size == 'r1e-06':
        fname = 'noadv_global_%s_%s_3D_grid2x2_allrho_allr_90days_60dtsecs_12hrsoutdt.nc' % (s, yr)
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
        fname = 'noadv_global_%s_%s_3D_grid2x2_rho%s_%s_90days_60dtsecs_12hrsoutdt.nc' % (s, yr, rho, size_fn)
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
    z_set = [0, ] * depths.shape[0]
    vs_max = [0, ] * depths.shape[0]
    for i in range(depths.shape[0]): #9620: number of particles 
        z = np.where(depths[i,:] > 1.)[0]
        z_set[i] = z[0] if z.any() else 0
              
        z2 = [0,] * depths.shape[1]
        for ii in range(2,depths.shape[1]):
            z2[ii-1] = depths[i,ii]-depths[i,ii-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
        zf = np.where(np.array(z2) < 0.)[:]
        zind = zf[0][0] if np.array(zf).any() else []     
        vs_max[i] = np.max(vs[i,0:zind]) if np.array(zind).any() else np.nan
            
    t_set = time[z_set]
    t_set[t_set == 0] = np.nan 
    
    
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
    
    # t_set_all = np.hstack(t_set)     #t_set_all[:,idy] = t_set #
    # vs_max_all = np.hstack(np.array([vs_max]).T) #.append(np.array(vs_max).T)[:] #vs_max_all[:,idy] = vs_max
    # vs_max_all = np.concatenate((vs_max_all,np.array([vs_max]).T)) #,axis=1)
        
Ts = np.nanmean(t_set_all,axis=1) #[:,idx]
Vs = np.nanmean(vs_max_all,axis=1)
    
    # with open(dirwrite+str(rho)+str(sz_fn)+'JJA2001_av_Ts_Zmax.pickle', 'wb') as f:
    #     pickle.dump([lons,lats,lon_set,lat_set,t_set,z_max,z_max2,i_noset,prob_set], f)

#%%    
    
fig1 = plt.figure(figsize=(20,10))

cmap = plt.cm.get_cmap('magma_r', 11) #RdBu #Blues #vidiris

m = Basemap(projection='robin',lon_0=-180,resolution='l') #, ax=axs[idx])
m.drawparallels(np.array([-60,-30,0,30,60]), labels=[True, False, False, True], linewidth=2.0, size=20,zorder=0)
m.drawmeridians(np.array([50,150,250,350]), labels=[False, False, False, True], linewidth=2.0, size=20,zorder=0)
m.drawcoastlines()
m.fillcontinents(color='lightgrey',zorder=3)
xs, ys = m(lons[:,0], lats[:,0])
scat = m.scatter(xs, ys, marker='.', c=Ts,cmap = cmap, vmin = 0, vmax = 90, s = 100,zorder=1) #scat = 
cbar = m.colorbar(label = 'days') #,ax=axs[idx])  
cbar.set_label(label='days', size='20')   
cbar.ax.tick_params(labelsize=20) 
plt.title(str(yr)+' no advection all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)
    

fig2 = plt.figure(figsize=(20,10))

cmap = plt.cm.get_cmap('magma_r', 11) #RdBu #Blues #vidiris

m = Basemap(projection='robin',lon_0=-180,resolution='l') #, ax=axs[idx])
m.drawparallels(np.array([-60,-30,0,30,60]), labels=[True, False, False, True], linewidth=2.0, size=20,zorder=0)
m.drawmeridians(np.array([50,150,250,350]), labels=[False, False, False, True], linewidth=2.0, size=20,zorder=0)
m.drawcoastlines()
m.fillcontinents(color='lightgrey',zorder=3)
xs, ys = m(lons[:,0], lats[:,0])
scat = m.scatter(xs, ys, marker='.', c=Vs,cmap = cmap, s = 100,zorder=1) #vmin = 0, vmax = 1e-6,
cbar = m.colorbar(label = 'm/s') #,ax=axs[idx])  
cbar.set_label(label='[m/s]', size='20')   
cbar.ax.tick_params(labelsize=20) 
plt.title(str(yr)+' all seas average maximum settling velocity [m/s], rho ='+str(rho)+', size ='+str(size) ,size = 20)
