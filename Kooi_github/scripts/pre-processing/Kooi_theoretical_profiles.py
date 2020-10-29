#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 12:10:22 2020

@author: Lobel001
"""

# V3 COMBINES ALL PROFILES THAT ARE TIME DEPENDENT INTO ONE LOOP TO SAVE AS PICKLE 

#------------ Kooi et al. 2017 kernals --------------
## TOTAL EQUATION FOR VERTICAL VELOCITY: Vs = -(((rho_tot - rho_sw)/rho_sw) * g * omega * upsilon_sw)**(1/3)

#from IPython import get_ipython
#get_ipython().magic('reset -sf') 

#%matplotlib qt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math as math
#import pandas as pd 
#import operator
import pickle
from numpy import *
#import scipy.linalg
#import scipy.integrate as integrate
#from math import e

plt.close("all")

def cube(x):
    if 0<=x: return x**(1./3.)
    return -(-x)**(1./3.)

# depth profile used in NEMO [metres]
nemo_z = [0, 1.023907, 2.10319, 3.251309, 4.485053, 5.825238, 7.297443, 
    8.932686, 10.7679, 12.84599, 15.21527, 17.92792, 21.03757, 24.59599, 
    28.64965, 33.23697, 38.3871, 44.12101, 50.45447, 57.40257, 64.9846, 
    73.2287, 82.17556, 91.88141, 102.4202, 113.8852, 126.3909, 140.074, 
    155.095, 171.6402, 189.9228, 210.1845, 232.697, 257.7629, 285.7158, 
    316.9199, 351.768, 390.6786, 434.0905, 482.4563, 536.2332, 595.8721, 
    661.8052, 734.4321, 814.1057, 901.118, 995.6885, 1097.954, 1207.963, 
    1325.672, 1450.95, 1583.582, 1723.28, 1869.693, 2022.425, 2181.044, 
    2345.101, 2514.137, 2687.699, 2865.347, 3046.659, 3231.24, 3418.723, 
    3608.769, 3801.072, 3995.354, 4191.367, 4388.89, 4587.726, 4787.702, 
    4988.667, 5190.488, 5393.049, 5596.249, 5800]

# CHOOSE THESE: 

rho_pl = 1050.                 # density of plastic (kg m-3): 840, 920, 940, 1050, 1380 (last 2 are initially non-buoyant)
r_pl = 10.**(-4)               # radius of plastic (m): ranges from 10 mm to 0.1 um or 10-2 to 10-7 m

# CONSTANTS
#rho_tot = rho_pl + rho_bf     # total density of plastic (kg m-3)
g = 7.32e10#/(86400**2)         # gravitational acceleration (m d-2), now [s-2]
k = 1.0306E-13#/(86400**2)      # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)
D_n = 2*r_pl                   # equivalent spherical diameter [m]


# INITIAL CONDITIONS
A_0 = 0                        # initial number of attached algae [no. m-2] - changes dynamically with Runge-Kutta integration
z = 0                          # initial depth [m]
t = 0                          # time? [days?] 

""" Temperature, salinity and density profiles: rho_sw [kg m-3] """ # using a Hill function (Hill et al. 1910)

# ------Temperature [C] (Eq. 22)------
surf_z = 0
tot_z = -4000       # total depth of ocean
z_step = -1
T_surf = 25.
T_bot = 1.5 
p = 2               # steepness temp decrease
z_c = -300          # depth of thermocline (m)


T_z = []
for z_i in range(len(nemo_z)): #range (surf_z,tot_z,z_step):
    z = nemo_z[z_i]
    T_z.append(T_surf + ((T_bot - T_surf)  * (z**p/(z**p + z_c**p))))
                     
#y = int(-tot_z+surf_z/-z_step)
depth = nemo_z #range(0,y)

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=1, ncols=5) #, sharex=True, sharey=True)
ax1.scatter(T_z,depth)
ax1.title.set_text('Temperature [C]')
ax1.set_ylabel('Depth [m]')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_xlim([min(T_z), max(T_z)])

    
    
# ------Salinity [g kg-1] (Eq. 24)------
# Using North Pacific coefficients for salinity: used in paper (Table S4)

c1 = 9.9979979767E-17
c2 = 1.0536246487E-12
c3 = 3.9968286066E-09
c4 = 6.5411526250E-06
c5 = 4.1954014008E-03
c6 = 3.5172984035E+01
s_fix = 34.6  # [g kg-1]  constant salinity from z > zfix
z_fix = -1000 # NOTE: depth must be negative for equation to work



#ind = (enumerate(nemo_z)*-1>z_fix, key=operator.itemgetter(1))

ind = nonzero(np.array(nemo_z)*-1>np.array(z_fix))

z_i = []
S_z_g = []
for z_i in range (size(ind)): #surf_z,z_fix, z_step):
    z = nemo_z[z_i]*-1
    S_z_g.append((c1*(z**5)) + (c2*(z**4)) + (c3*(z**3)) + (c4*(z**2)) + (c5*z) + c6)


# to add a linear fit from 1000 to 2000 m 

idx = (np.array(nemo_z)<2000)*(np.array(nemo_z)>1000)
n = np.where(idx)
S_end = S_z_g[-1]
s_ = np.linspace(S_end,s_fix,size(n)+2)  

# to add the fixed salinity below 2000m (in paper it is 1000 but in excel, linear interp to 2000m)
S_rest = size(depth)- size(s_) - size(S_z_g)+1                                                     
s = np.array([s_fix] * S_rest) 
S_z_g = np.concatenate((S_z_g[0:-1],s_,s))

#y2 = int(-z_fix/-z_step)
#depth2 = range(0,y2)

ax2.scatter(S_z_g,depth) # depth2
ax2.title.set_text('Salinity [g kg-1]')
ax2.set_ylim(ax2.get_ylim()[::-1])
ax2.set_xlim([min(S_z_g), max(S_z_g)])


#------ Density profile [kg m-3] (Eq. 23)------
a1 = 9.999E2
a2 = 2.034E-2
a3 = -6.162E-3
a4 = 2.261E-5
a5 = -4.657E-8

b1 = 8.020E2
b2 = -2.001
b3 = 1.677E-2
b4 = -3.060E-5 #2.261E-5 # THIS WAS A TYPO IN TABLE S1 AND KOOI CODE
b5 = -4.657E-5

rho_sw = []

S_z = S_z_g/1000 # NOTE: salinity must be in kg/kg instead of g/kg for equation to work

for i in range(len(depth)):
    rho_sw.append((a1 + (a2*T_z[i]) + (a3*(T_z[i]**2)) + (a4*(T_z[i]**3)) + 
                   (a5*(T_z[i]**4))) + ((b1*S_z[i]) + (b2*S_z[i]*T_z[i]) + 
                   (b3*S_z[i]*(T_z[i]**2)) + (b4*S_z[i]*(T_z[i]**3)) + (b5*(S_z[i]**2)*(T_z[i]**2))))
    
ax3.scatter(rho_sw,depth) # depth2
ax3.title.set_text('Density [kg m-3]')
ax3.set_ylim(ax3.get_ylim()[::-1])
ax3.set_xlim([min(rho_sw), max(rho_sw)])

mu_w = []  # dynamic viscosity [kg m-1 s-1]
A = []
B = []
mu_sw = []
upsilon_sw = []



""" Kinematic viscosity: upsilon_sw [m2 s-1]: Eq. 25 to 29"""
for ii in range(len(depth)):
    mu_w.append(4.2844E-5 + (1/((0.157*(T_z[ii] + 64.993)**2)-91.296))) # kg m-1 s-1
    A.append(1.541 + 1.998E-2*T_z[ii]- 9.52E-5*T_z[ii]**2)              # 
    B.append(7.974 - 7.561E-2*T_z[ii] + 4.724E-4*T_z[ii]**2)
    mu_sw.append(mu_w[ii]*(1+ A[ii]*S_z[ii] + B[ii]*S_z[ii]**2))
    upsilon_sw.append(mu_sw[ii]/rho_sw[ii])


ax4.scatter(upsilon_sw,depth) # depth2
ax4.title.set_text('Kinematic viscosity [m2 s-1]')
ax4.set_ylim(ax4.get_ylim()[::-1])
ax4.set_xlim([min(upsilon_sw), max(upsilon_sw)])
ax4.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


ax5.scatter(mu_sw,depth) # depth2
ax5.title.set_text('Seawater viscosity [m2 s-1]')
ax5.set_ylim(ax5.get_ylim()[::-1])
ax5.set_xlim([min(mu_sw), max(mu_sw)])
ax5.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
plt.show()


 # TO SAVE FIRST PROFILES IN SAME PICKLE
#with open('profiles.pickle', 'wb') as f:
#    pickle.dump([depth,T_z,S_z,rho_sw,upsilon_sw,mu_sw], f)
    
    
""" Algal growth: Eq. 3 to 21"""

#------ Algae properties [Table S1]------
rho_bf = 1388.           # density of biofilm (kg m-3)
V_A = 2.0E-16            # Volume of 1 algal cell [m-3]
R_A = 0.1                # respiration rate [d-1] 
Q10 = 2.                 # temperature coefficient respiration
mu_max = 1.85#/86400      # maximum growth rate algae [d-1], now [s-1]
alpha = 0.12#/86400       # initial slope [d-1], now [s-1]
T_min = 0.2              # min temp algal growth [oC]
T_opt = 26.7             # optimal algal growth temp [oC]
T_max = 33.3             # max temp algal growth [oC]
I_opt = 1.75392E13#/86400 # optimal light intensity algal growth [uE m-2 d-1], now [s-1]
gamma = 1.728E5#/86400    # shear [d-1], now [s-1]


#------ Light profile [Table S1]: I tried different light extinction coefficients to check why the euphotic zone was so shall in Kooi------
I_m = 1.2E8#/86400        # Surface light intensity at noon [uE m-2 d-1], now [s-1]
k_w = 0.2 #0.1    #0.2            # extinction coefficient water [m-1]
k_p = 0.02   #0.12            # extinction coefficient algae [m-1 mg-1 chl L-1]


#------ Variables computed below using equations 3 to 21------
V_pl = []           # volume of plastic [m3]
theta_pl = []       # surface area of plastic particle [m2]
r_tot = []          # total radius of plastic + biofilm thickness
r_A = []            # radius of individual particle [m]
D_pl = []           # diffusivity of plastic particle [m2 s-1]
D_A = []            # diffusivity of algae [m2 s-1]
beta_Abrown = []    # Brownian motion [m3 s-1]
beta_Aset = []      # differential settling [m3 s-1]
beta_Ashear = []    # advective shear [m3 s-1]
deltaA = []         # attached algal growth
V_bf = []           # volume of biofilm [m3]
V_tot = []          # volume of total [m3]
t_bf = []           # biofilm thickness [m]
rho_tot = []        # total density [kg m-3]
Dstar = []          # dimensionless particle diameter
omega_star = []     # dimensionless settling velocity
epsilon = []        # light extinction coefficient [m-1]
I_0 = []            # light availability at the sea surface [uE m-2 s-1]
I_z = []            # light intensity at depth z [uE m-2 s-1]
phi = []            # temperature influence on algal growth 
mu = []             # algal growth [s-1]


# ----------- Volumes -----------------
V_pl = (4/3)*math.pi*r_pl**3 # [m3]

theta_pl = 4*math.pi*r_pl**2 # surface area of plastic particle [m2]

V_bf = (V_A*A_0)*theta_pl # [m3]

V_tot = V_bf + V_pl # [m3]

t_bf = cube(V_tot*(3/(4*math.pi)))-r_pl # [m] #V_tot*(3/(4*math.pi))**(1/3) - r_pl

r_tot = r_pl + t_bf #[m]

rho_tot = ((r_pl**3) * rho_pl + ((r_pl + t_bf)**3 - (r_pl**3))*rho_bf)/((r_pl + t_bf)**3) # [kg m-3]

theta_tot = 4*math.pi*r_tot**2 #[m2]

#Dstar = ((rho_tot - rho_sw)*g*D_n**3)/(rho_sw*upsilon_sw**2) 

#omega_star = 1.74E-4*Dstar**2 


# --------------- Chlorophyll profile (same as Uitz et al. 2006), proxy for A_A (ambient algal concentrations [no. m-3]): 
# using chl surf conc of 0.04-0.08 mg m-3 as default from table S2 -----------

Chla_Zbase = 0.151 # mg m-3
Cb = 0.533         # -
Cmax = 1.194       # - 
zmax = 92.01       # m
deltaz = 43.46     # m
s2 = 1.72E-03        # m-1
euph_z = 120 
 
C = np.zeros(len(depth))
depth2 = []
e_z2 = euph_z*2
d2 = np.array(depth)<e_z2
depth2 = np.array(depth)[d2]
C2 = []

for zz in range (len(depth2)):  
    z = depth2[zz]
    C[zz] = ((Cb + Cmax*math.exp(-1*((z-zmax)/deltaz)**2)- s2*z)*Chla_Zbase)
    C2.append((Cb + Cmax*math.exp(-1*((z-zmax)/deltaz)**2)- s2*z)*Chla_Zbase)


plt.figure(2)
plt.scatter(C,depth)
plt.xlabel('Kooi code Chla')
plt.ylabel('Depth [m]')
plt.ylabel('Depth [m]')
plt.gca().invert_yaxis()



#------------Profiles: light, carbon, ambient algae, algal growth (t,z): using 24 hours and 75 depth levels--------------

shift = 0 #0.5 #I_m/2 #0.5 # in the excel for Merel's test of changing num of hrs of light in a day
days = 1 #3
hours = 24 # mins = 24*60*days
time = np.linspace(0,days,hours) #mins)

# Defining the euphotic layer depth (1% of surf irradiance at midday)
epsilon = k_w + (k_p*0.05) # Merel only uses surface concentration Chla 
Iz_max = []
for zz in range (len(depth)):
    z = depth[zz]
    Iz_max.append(I_m*math.exp(-1*epsilon*z))

idx = nonzero(np.array(Iz_max)<I_m*0.01)
euph_z_true = depth[idx[0][0]]

z =[]
zz = []
t = []
tt = []
I_0 = [] 
#I_z_t = []
#Carbon_t = []
#A_A_t = [] 
#mu_A_t = []

Carbon = np.zeros((len(time),len(depth)))
A_A = np.zeros((len(time),len(depth)))
phi = np.zeros((len(time),len(depth)))
mu_opt = np.zeros((len(time),len(depth)))
mu_A = np.zeros((len(time),len(depth)))
I_z = np.zeros((len(time),len(depth)))


for tt in range (0,hours): #mins): 
    t = time[tt]    
    e2 = 2*math.pi*t
    if I_m*(math.sin(e2)+shift)<0:
        I_0.append(0)
    else:
        I_0.append(I_m*(math.sin(2*math.pi*t)+shift))
                
    for zz in range (len(depth2)): #0,euph_z[0]*2): #(len(depth)):
        #Chla_int.append(np.trapz(C[0:z])) # np.cumsum(C) # 
        #epsilon.append(k_w + (k_p*Chla_int[z])) #C[z]))
        z = depth2[zz]
        I_z[tt,zz] = (I_0[tt]*math.exp(-1*epsilon*z)) 
        Carbon[tt,zz] = (C[zz]/(0.003+1.0154*math.exp(0.050*T_z[zz])*math.exp(-0.059*I_z[tt,zz]/1E6)))  
        if Carbon[tt,zz]/(2726*1E-9)<0:
            A_A[tt,zz] = 0
        else:
            A_A[tt,zz] = (Carbon[tt,zz]/(2726*1E-9))
            
        phi[tt,zz] = (((T_z[zz] - T_max)*(T_z[zz] - T_min)**2)/((T_opt-T_min)*((T_opt-T_min)*
                     (T_z[zz]-T_opt)-(T_opt-T_max)*(T_opt+T_min-(2*T_z[zz])))))
        
        if T_z[zz]<T_min:
            mu_maxT = 0
        elif T_z[zz]>T_max:
            mu_maxT = 0
        else:
            mu_maxT = mu_max 
#        
        if z > euph_z_true: #Iz_max[zz]<I_m*0.01: # since euphotic zone is defined as depth of 1% of midday surf light
            mu_opt[tt,zz] = 0 
        else:                             
            mu_opt[tt,zz] = (mu_maxT*(I_z[tt,zz]/(I_z[tt,zz] + (mu_maxT/alpha)*((I_z[tt,zz]/I_opt)-1)**2)))
        
        mu_A[tt,zz] = (mu_opt[tt,zz]*phi[tt,zz])

        
    plt.figure(4)
    plt.scatter(I_z[tt,:],depth)
    
    plt.figure(5)
    plt.scatter(Carbon[tt,:],depth)
    
    plt.figure(6)
    plt.scatter(A_A[tt,:],depth)
    
    plt.figure(7)
    plt.scatter(mu_A[tt,:]/86400,depth)
    
    
#    a2.scatter(I_z[tt,0:len(depth2)],depth2)
#    a3.scatter(A_A[tt,0:len(depth2)],depth2)
#    a4.scatter(mu_A[tt,0:len(depth2)],depth2)
#    
    #I_z_t.append(I_z)
#    Carbon_t.append(Carbon)
#    A_A_t.append(A_A)
#    mu_A_t.append(mu_A)
    
#    I_z_t2[tt,] = I_z


plt.figure(3)
plt.scatter(time,I_0) 
plt.xlabel('Time [days]')
plt.ylabel('Light availability [uE m-2 min-1]') # double check units 
plt.title('Light availability at the surface over 1 day [uE m-2 d-1]')

plt.figure(4)
plt.xlabel('Light availability [uE m-2 d-1]')
plt.ylabel('Depth [m]')
plt.title('Light availability for every surface I_0 [uE m-2 d-1]')
plt.gca().invert_yaxis()

plt.figure(5)
plt.xlabel('Carbon [mg C m-3]')
plt.ylabel('Depth [m]')
plt.title('Carbon')
plt.gca().invert_yaxis()

plt.figure(6)
plt.xlabel('Ambient algae [no. m-3]')
plt.ylabel('Depth [m]')
plt.title('Ambient algae')
plt.gca().invert_yaxis()

plt.figure(7)    
#plt.scatter(mu_A,depth)
plt.xlabel('Algal growth [s-1]') #[d-1]
plt.ylabel('Depth [m]')
plt.title('Algal growth')
plt.gca().invert_yaxis()

A_A_t = A_A
mu_A_t = mu_A

#fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=4) #, sharex=True, sharey=True)
#C2 = np.array(C)[d2]
#a1.scatter(list(C2),list(depth2))

#ax1.scatter(T_z,depth)
#ax1.title.set_text('Temperature [C]')
#ax1.set_ylabel('Depth [m]')
#ax1.set_ylim(ax1.get_ylim()[::-1])
#ax1.set_xlim([min(T_z), max(T_z)])



#with open('profiles_t.pickle', 'wb') as p:
#    pickle.dump([depth,time,A_A_t,mu_A_t], p)

####### STOPPED FROM HERE ON SINCE REALISED THAT THE TIME DERVIATIVE WILL BE EASIER TO COMPUTE FROM WITHIN PARCELS 

