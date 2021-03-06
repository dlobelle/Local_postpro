#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 18:15:42 2020

@author: Lobel001
"""
# 30/07/20- This has been written for 920 kg m-3 only. 

#from math import sin, cos, sqrt, atan2, radians
#from geopy import distance

import matplotlib 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.matlib 
#import math
#from pylab import *
#import cmocean
import os #, fnmatch 
#import cartopy
#import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")

rho = '920' # [kgm-3]: density of the plastic 
res = '2x2' # [deg]: resolution of the global release of particles
# size = 'r1e-07' # [m]: size of plastic
loc = 'global'

# if size == 'r1e-02':
#     size_fn = 'r0.01'
# else:
#     size_fn = size 

# dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
# dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
# dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'

plt.close("all")

#%%

# yr_st = 2004
# numyrs = 1
# yrs = np.arange(yr_st,yr_st+numyrs)
seas = 'MAM' #'DJF', 'MAM', 'JJA', 'SON']
s = seas
yr = 2004

fig = plt.figure(figsize=(20,20), constrained_layout=True)
gs = fig.add_gridspec(figure = fig, nrows = 4, ncols = 2, height_ratios = [6,6,6,1])
plt.rc('font', size = 24)

for ids in range(2,8): #sizes
    size = 'r1e-0'+str(ids) #'r1e-02'
    if size == 'r1e-02':
        size_fn = 'r0.01'
    else:
        size_fn = size #'r1e-02'
    
    print(size)
    dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
    dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    
    if size == 'r1e-05' or size == 'r1e-07':
        num_part = 9582
    else:
        num_part = 9620

        
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
                else:
                    lons = np.vstack((lons,data.variables['lon'][:])) 
                    lats= np.vstack((lats,data.variables['lat'][:]))
                    depths =np.vstack((depths,data.variables['z'][:]))
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
    
    time = time - time[0]


    """ To find the horiz. distance until the first depth it sinks to """
    # time_p = np.zeros((depths.shape[0],depths.shape[1])) #,time.shape[0]-1)) #lons.shape[0]
    # time_p[:] = np.nan  
    
    # depth_p = np.zeros((depths.shape[0],depths.shape[1])) #,time.shape[0]-1)) #lons.shape[0]
    # depth_p[:] = np.nan  
    #zind = np.zeros((depths.shape[0]))

    boo = np.zeros((depths.shape[0],depths.shape[1]))
    depths2 = []
    zind = []
    for i in range(depths.shape[0]): #9582: number of particles 
        
        z2 = [0,] * depths.shape[1]
        depths2 = depths[i,:]
        if any(depths2>1.):
            f = np.where(depths2>1.)
            for ii in range(f[0][0],depths.shape[1]): #(2,depths.shape[1]):
                if depths2[ii]>0.: #1.
                    z2[ii-1] = depths[i,ii]-depths[i,ii-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
            zf = np.where(np.array(z2) < 0.)[:]           
            zind = zf[0][0] if np.array(zf).any() else [] #np.nan
            if np.array(zind).any():
                j = np.array(zind)
                boo[i,:j+1] = 1.
            # if np.array(zind).any(): 
            #     time_p[i,:] = time[0:np.array(zind)]
            #     depth_p[i,:] = depths[i,0:np.array(zind)]
            #     for t in range(np.array(zind)): #time.shape[0]-1): 
            #         #sinks = np.where(depths[0,:]>1.0)
            #         if lats[i,t+1].any(): #and np.array(sinks).any():
            #             loc1 = (lats[i,t],lons[i,t])
            #             loc2 = (lats[i,t+1],lons[i,t+1])                        
            #             dist[i,t] = distance.geodesic(loc1,loc2).km 
            #         else:
            #             dist[i,t] = np.nan                   
                        
 #%%                   
               
    # distcum = np.cumsum(distnan,axis = 1)
    # n = np.array(([[np.nan, ]] * distcum.shape[0]))
    # dist_p = np.append(distcum,n, axis= 1)
    
    #lats_p = np.tile(lats[:,0],(2,1)).T
    
    ''' median depth and distance'''
    #z_med = np.nanmedian(depths[depths>1].ravel().data)
    # dist_med = np.nanmedian(dist_p.ravel().data)
    
    ''' to get the cmap for line plot below, need to get it from scatterplot '''
    plot_idx = np.random.permutation(depths.shape[0])
    
    time_p = np.tile(time,(depths.shape[0],1))#.T
    lats_p = np.tile(lats[plot_idx,0],(lats.shape[1],1)).T
    depths_p = depths[plot_idx,:]
    boo2 = np.array(boo[plot_idx,:],dtype='bool')

    
    fig1 = plt.figure(figsize=(15,10))
    cmap = plt.cm.get_cmap('coolwarm',7)
    scat = plt.scatter(time_p[boo2],(depths_p[boo2]*-1), vmin = -70, vmax = 70, c = lats_p[boo2], cmap = cmap, alpha = 0.6) #, cbarlabel = 'initial latitude') # 
    #ax = plt.gca()
    #ax.set_facecolor('darkgrey') 
    plt.colorbar()
    plt.ylim(top=0, bottom =-240) 
    plt.xlim(left=0, right = 90)
    ax = plt.gca()
    #ax.set_facecolor('lightgrey')   
    plt.ylabel('Depth [m]', size = 20)
    plt.xlabel('Time [days]', size = 20)
    font = {'size'   : 20} ##'family' : 'normal',
            #'weight' : 'bold',
    plt.title(str(s)+' First sinking depth reached within 90 days by rho rho ='+str(rho)+', size ='+str(size) ,size = 20)
    
    matplotlib.rc('font', **font)
    #%% 20/07/20- Line plot with marker depth vs horizontal distance 

    #fig2 = plt.figure(figsize=(15,10))
    
    
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

    ax = fig.add_subplot(gs[r,c])
    
    
    ''' using colorbar above to separate colours by initial release latitudinal bins'''
    #plot_idx = np.random.permutation(dist_p.shape[0])
    #lats_p = lats
    
    time_save = np.zeros(plot_idx.shape[0])
    time_save[:] = np.nan
    depths_save = np.zeros(plot_idx.shape[0])
    depths_save[:] = np.nan  
    lats_save = np.zeros(plot_idx.shape[0])
    d = np.zeros((depths.shape[0],depths.shape[1]))
    d[:] = np.nan 
    for ii in range(plot_idx.shape[0]):
        # ii = plot_idx[i]
        boo_p = boo2[ii,:]
        
        if lats_p[ii,0]<-50.:
            rgb = cmap(0)[:3]
        elif lats_p[ii,0]<-30. and lats_p[ii,0]>=-50.:
            rgb = cmap(1)[:3]
        elif lats_p[ii,0]<-10. and lats_p[ii,0]>=-30.:
            rgb = cmap(2)[:3]
        elif lats_p[ii,0]<10. and lats_p[ii,0]>=-10.:
            rgb = cmap(3)[:3]            
        elif lats_p[ii,0]<30. and lats_p[ii,0]>=10.:
            rgb = cmap(4)[:3]
        elif lats_p[ii,0]<50. and lats_p[ii,0]>=30.:
            rgb = cmap(5)[:3]
        elif lats_p[ii,0]>=50.:
            rgb = cmap(6)[:3]
  
        """To plot a scatter dot for the first max depth"""
        # d0 = depths_p[ii,boo_p]*-1 
        # d[ii,0:d0.shape[0]] = d0 if d0.shape[0]>0. else 0.
        
        # d2 = np.zeros((depths.shape[0],depths.shape[1]))
        # for zz in range(1,depths.shape[1]):
        #     d2[ii,zz] = d[ii,zz]-d[ii,zz-1] #np.where(depths[i,ii]<depths[i,ii-1])[0]
        #     zf2 = np.where(np.array(d2) > 0.)[:] 
            
            
        ax.plot(time_p[ii,boo_p],(depths_p[ii,boo_p]*-1), c = rgb, linewidth=2, alpha = 0.6)# 0.6) #alpha = 0.6, 
        ind_nonan = np.where(depths_p[ii,boo_p]>1.)
        if np.array(ind_nonan).any():
            last_ind = ind_nonan[0][-1]
            ax.plot(time_p[ii,last_ind],(depths_p[ii,last_ind]*-1), marker = 'o', c = rgb, markersize=10, linewidth = 0.5, markeredgecolor='black',alpha = 0.7) 
        
            time_save[ii] = time_p[ii,last_ind]
            depths_save[ii] = depths_p[ii,last_ind]
            lats_save[ii] = lats_p[ii,0]
    
    ''' median depth and distance'''
    z_med = np.nanmedian(depths_save.ravel().data)
    time_med = np.nanmedian(time_save.ravel().data)
        
    ax.axvline(x=time_med, color = 'k', linewidth = 4)
    ax.axhline(y=z_med*-1, color = 'k', linewidth = 4)
    ax.set_ylim(top=0, bottom =-240) 
    ax.set_xlim(left=0, right = 90) #2500)
    #ax = plt.gca()  
    if ids == 4:
        ax.set_ylabel('Depth [m]', size = 26)
    if ids == 6 or ids == 7:
        ax.set_xlabel('Time [days]', size = 26)
    # font = {'size'   : 20} ##'family' : 'normal',
            #'weight' : 'bold',
    
    ax.set_title(letter+ ' radius = '+str(size)+'m') #, fontsize = 26)
    #ax.set_title(str(s)+' First sinking depth reached within 90 days by rho ='+str(rho)+', size ='+str(size) ,size = 20)
    
    ax = plt.gca()
    plt.rc('font', size = 24)

cbaxes = fig.add_axes([0.1, 0.04, 0.8, 0.01])
cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='both', label = 'initial latitude')
plt.rc('font', size = 24)

fig.savefig('/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/Zs.pdf')

# with open('/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/Matias_test.pickle', 'wb') as f:
#     pickle.dump([time_save,depths_save,lats_save], f)

#fig.text(0.5, 0.04, 'Time [days]', ha='center', va='center', fontsize = 24)
#fig.text(0.06, 0.5, 'Depth [m]', ha='center', va='center', rotation='vertical', fontsize = 24)

#%% 30/07/20- No longer want scatterplot
# plot_idx = np.random.permutation(dist_p.shape[0])

# fig = plt.figure(figsize=(15,10))
# cmap = plt.cm.get_cmap('coolwarm',7) #cmocean.cm.delta,8)
# # cm = plt.cm.get_cmap('viridis',len(years))

# s = plt.scatter(dist_p[plot_idx,:],(depths[plot_idx,:]*-1), vmin = -70, vmax = 70, c = lats_p[plot_idx], cmap = cmap) #, alpha = 0.4) #, cbarlabel = 'initial latitude') # 
# #plt.plot(dist_p[plot_idx,:],(depths[plot_idx,:]*-1),'k')
# plt.axvline(x=dist_med, color = 'k', linewidth = 3)
# plt.axhline(y=z_med*-1, color = 'k', linewidth = 3)
# plt.ylim(top=0) 
# plt.xlim(left=0)
# plt.ylabel('Depth [m]', size = 18)
# plt.xlabel('Horizontal distance [km]', size = 18)
# plt.colorbar(label = 'initial latitude')    
# #ax = plt.gca()
# #ax.set_facecolor('grey')    

  
#%% 30/07/20- This worked only if the initial latitudes were equidistant, however there are more particles/ocean area in southern hemisphere
# ind_sortlat = numpy.argsort(lats[:,0])
# cmap_sort= cmocean.cm.delta(np.linspace(0, 1, 9620))

# fig2 = plt.figure(figsize=(15,10))
# plt.colorbar(s)
# for ii in range(dist2.shape[0]):
#     order = ind_sortlat[ii]
#     plt.plot(dist_p[order,:],(depths[order,:]*-1), '-o', c = cmap_sort[ii])

# plt.ylim(top=0) 
# plt.xlim(left=0)

#fig2.colorbar(cmap_sort)
#%%

# for             
#     plt.plot(dist_p[500],(depths[500]*-1),'-o')            

# for idx,s in enumerate(seas): #, 'JJA', 'SON']): # the 4 seasons 
#     for idy,yr in enumerate(yrs): #yrs): #2011):
#         print(s, yr)
#         fname = 'global_%s_%s_3D_grid2x2_rho%s_%s_90days_60dtsecs_12hrsoutdt.nc' % (s, yr, rho, size_fn)
#         if not os.path.isfile(dirread+fname):
#             print('%s not found' %fname)
#         else:          
#             data = Dataset(dirread+fname,'r') 
#             time = data.variables['time'][0,:]/86400
#             lons=data.variables['lon'][:] 
#             lats=data.variables['lat'][:] 
#             depths =data.variables['z'][:]
    
#             time = time - time[0]
            
#             dist = np.zeros((depths.shape[0],time.shape[0]-1)) #lons.shape[0]
            
            
#             for i in range(lons.shape[0]):
#                 for t in range(time.shape[0]-1):                   
#                     if lats[i,t+1].any():
#                         loc1 = (lats[i,t],lons[i,t])
#                         loc2 = (lats[i,t+1],lons[i,t+1])                        
#                         dist[i,t] = distance.geodesic(loc1,loc2).km 
#                     else:
#                         dist[i,t] = np.nan