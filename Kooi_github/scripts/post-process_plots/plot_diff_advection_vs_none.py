#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 22:25:00 2020

@author: Lobel001
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np 
import cartopy
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

#import math
#from pylab import *
#import cmocean
import os #, fnmatch 
#import cartopy
import pickle 
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore", "Mean of empty slice")



rho =  '920' #np.array([920])  # [kgm-3]: density of the plastic 
res = '2x2' # [deg]: resolution of the global release of particles
loc = 'global'

plt.close("all")

# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return idx # array[idx]

def getclosest_ij(lats,lons,latpt,lonpt):     
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_sq = (lats-latpt)**2 + (lons-lonpt)**2                 # find squared distance of every point on grid
    minindex_flattened = dist_sq.argmin() 
    #return np.unravel_index(minindex_flattened, lats.shape) 
    return minindex_flattened
#%%

yr = 2004
numyrs = 1
seas = ['DJF', 'MAM', 'JJA', 'SON']

cmap = plt.cm.get_cmap('magma_r', 9) #RdBu #Blues #vidiris
row = 1 
col = 1

''' pickle from plot_Ts_Vs_2004_allsizes.py'''

with open('/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/Ts2004_allsizes_for_diffplot_with_advection.pickle','rb') as f:  
    lons_adv,lats_adv,Ts02,Ts03,Ts04,Ts05,Ts06,Ts07 = pickle.load(f)
    
    
''' prepare cartopy projection and gridspec'''

projection = cartopy.crs.PlateCarree() #central_longitude=72+180)

size_st = 2 # this is the superscript for 1e-0+'size1', hence starts on 1e-02
numsizes = 6
fig_w = 8 #10 #30 # 30 25 #28
fig_h = 7#8 #25 #30 30 #18

#plt.rc('font', size = 12)
fig = plt.figure(figsize=(fig_w,fig_h), constrained_layout=True) #True) # 10,5
gs = fig.add_gridspec(figure = fig, nrows = 4, ncols = 2, height_ratios=[7,7,7,1]) #,hspace = 0.5, wspace = 0.5) #, height_ratios=[10,1,10,1], wspace = 0.05, hspace = 1) # gridspec.GridSpec

dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/post_pro/'
 
for ids in range(size_st,size_st+numsizes):#8): #sizes
    size = 'r1e-0'+str(ids) #'r1e-02'
    if size == 'r1e-02':
        size_fn = 'r0.01'
    else:
        size_fn = size #'r1e-02'
        
    #dirread = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/data_output/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    #dirwrite = '/Users/Lobel001/Desktop/Local_postpro/Kooi_data/post_pro_data/'
    #dirwritefigs = '/Users/Lobel001/Desktop/Local_postpro/Kooi_figures/rho_'+rho+'kgm-3/res_'+res+'/'+size+'/'
    
    # if size == 'r1e-05' or size == 'r1e-07':
    #     num_part = 9582
    # else:
    #     num_part = 9620
    
    Ts_adv = eval('Ts0'+str(ids))
    num_part = 9620   
    t_set_all = np.zeros((num_part,len(seas))) # just keeping it to have the same num of variables in the pickle. 
    vs_max_all = np.zeros((num_part,len(seas)))
    
    for idy,s in enumerate(seas):
    
        print(size, s, yr)
            

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

                
        # vs = data.variables['vs'][inds]
        time = time - time[0]
        ''' find the first time and max depth where particles sink (below 1 m, since that's the initial release depth) '''      
        t_set = np.zeros(depths.shape[0]) 
        t_set[:] = np.nan 
        
        #t_set = [0, ] * depths.shape[0]
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
            
                t_set[i] = time[z_set] 
                            
                ''' removed section on finding vs since no longer in paper '''
            
        t_set = np.array(t_set)
        #t_set[t_set==0]= np.nan
       
        t_set_all[:,idy] = t_set  
        
    
    #t_set_all[t_set_all==0] = np.nan  

             
    Ts = np.nanmean(t_set_all,axis=1) #[:,idx]

    ''' Slicing needed: sims with adv for these 2 sizes (those run with 3 lon bands) are 9582 particles and lat mask is from -70:78 N, whereas all others are 9620 particles and lat mask -68:80 N'''
    if ids == 5 or ids == 7: 
        ind_lat_noadv = np.argwhere(lats[:,0] <=78.1)
        ind_lat_adv = np.argwhere(lats_adv[:,0]>=-68)
        
        Ts_noadv = Ts[ind_lat_noadv[:,0]]
        Ts_adv0 = Ts_adv[ind_lat_adv[:,0]]        
        #Ts = Ts_adv0 #np.flip(Ts_adv0)
        
        lats_noadv0 = lats[ind_lat_noadv[:,0],0]
        lons_noadv0 = lons[ind_lat_noadv[:,0],0]
        
        lats_adv0 = lats_adv[ind_lat_adv[:,0],0] #np.flip(
        lons_adv0 = lons_adv[ind_lat_adv[:,0],0] #np.flip(
        
        ind2 = np.zeros(lats_noadv0.shape)
        for ii in range(np.array(lats_noadv0.shape[0])):
            latpt = lats_noadv0[ii]
            lonpt = lons_noadv0[ii]
            
            #ind = np.argwhere(np.min(lats_adv0-latpt)) and np.min(lons_adv0-lonpt)) 
            #ind2[ii] = find_nearest(lats_adv0, latpt) and find_nearest(lons_adv0, lonpt) 
            ind2[ii] = getclosest_ij(lats_adv0,lons_adv0,latpt,lonpt)
        ind1 = ind2.astype(int)
        
        Ts_diff = abs(Ts_noadv-Ts_adv0[ind1])
        lats_plot = lats_adv0[ind1] #[ind_lat_noadv[:,0],0]
        lons_plot = lons_adv0[ind1] #[ind_lat_noadv[:,0],0]
        
        
        
        #Ts_diff = abs(Ts_noadv-Ts_adv0)
    else:
        lats_plot = lats[:,0]
        lons_plot = lons[:,0]
        Ts 
        Ts_diff = abs(Ts-Ts_adv)
     
#%%    
    if ids == 2:
        letter = '(a)'
        r = 0
        c = 0
        size_name = '10 mm'
    if ids == 3:
        letter = '(b)'
        r = 0
        c = 1
        size_name = '1 mm'
    if ids == 4:
        letter = '(c)'
        r = 1
        c = 0
        size_name = '0.1 mm'
    if ids == 5:
        letter = '(d)'
        r = 1
        c = 1
        size_name = '10 \u03bcm'
    if ids == 6:
        letter = '(e)'
        r = 2
        c = 0
        size_name = '1 \u03bcm'
    if ids == 7:
        letter = '(f)'
        r = 2
        c = 1
        size_name = '0.1 \u03bcm'
        
    #plt.rc('font', size = 12)   
    ax = fig.add_subplot(gs[r,c], projection=projection)
    ax.coastlines(resolution='50m',zorder=3)
    ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=2)
    ax.set_ylim([-70, 80]) #,crs=cartopy.crs.PlateCarree())
    #ax.set_extent([73, 72, -70, 60])
    
    scat = ax.scatter(lons_plot, lats_plot, marker='.', c=Ts_diff,cmap = cmap, vmin = 0, vmax = 90, s = 4,zorder=1) #scat = 
    # cbar = plt.colorbar(scat, label = 'days') #,ax=axs[idx])  
    # cbar.set_label(label='days', size='20')   
    # cbar.ax.tick_params(labelsize=20) 
    plt.title(letter+ ' radius = '+size_name)#, fontsize = 26) #str(yr)+' all seas average settling onset time [days], rho ='+str(rho)+', size ='+str(size) ,size = 20)
    
    #plt.rc('font', size = 10)   
    cbaxes = fig.add_axes([0.1, 0.06, 0.8, 0.01])
    cbar = plt.colorbar(scat, cax=cbaxes, orientation="horizontal", aspect=50, extend='max')
    cbar.set_label(label='[days]') #, size=18)    
    
    print(size + ' median diff between adv and no adv = '+str(np.nanmedian(Ts_diff)))
    #print(size + 'max diff between adv and no adv = '+str(np.nanmax(Ts_diff)))
    
plt.savefig(dirwritefigs + 'Ts_diff_adv_noadv.pdf')