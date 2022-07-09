#!/usr/bin/env python
# coding: utf-8

# # O2 Budget
# 

# This notebook analyzes the O$_2$ budget from POP output. Dissolved oxygen in the ocean's interior is simulated following:
# 
# $$
# \frac{\partial{O_2}}{\partial{t}}= \underbrace{- \frac{\partial{U.O_2}}{\partial{x}} -\frac{\partial{V.O_2}}{\partial{y}}}_\text{Lateral Advection}
# - \overbrace{\frac{\partial{W.O_2}}{\partial{z}}}^\text{Vertical Advection}
# + \underbrace{A_h.\nabla^2{O_2}}_\text{Lateral Mixing}
# +\overbrace{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}^\text{Vertical Mixing}
# + \underbrace{ J(O_2)  }_\text{Sources - Sinks}
# $$
# 
# 

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
get_ipython().run_line_magic('matplotlib', 'inline')

import warnings
warnings.simplefilter("ignore") # Silence warnings

import xarray as xr
import numpy as np
from tqdm import tqdm
import xgcm 
import pop_tools

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import cartopy
import cartopy.crs as ccrs
import proplot 
import matplotlib.patches as patches

from utils import *


# In[18]:


C=CLSTR(1,15,200,10)


# In[19]:


C


# In[4]:


dsh=xr.open_mfdataset({'/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/O2_Budget_Mean.nc',
                      '/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/O2_Budget_Eddy_Mean.nc',})

du=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/CLM/[UVW]VEL.CLM.nc')
dshm=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/[UVW]_dO2m_Mean.nc')
dshd=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/[UVW]_dO2_Mean.nc')
dshe=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/DIVde_Mean.nc')

# dshm=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/DIVd_Mean.nc')
# dshd=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/DIVdm_Mean.nc')
# dshe=xr.open_mfdataset('/glade/scratch/yeddebba/Mesoscale/HR/TPAC/final_budget/DIVde_Mean.nc')


# In[77]:


fig = plt.figure(figsize=(25,7),dpi=100,)
gs = fig.add_gridspec(nrows=2, ncols=3, left=0.05, right=0.95, wspace=0.2, hspace=0.4)

fs=12; fsx=5; 
SS=-20; NN=20
tops=0. ; bots=500.; 
coords=[[bots,tops],[SS,NN]]

om_clr='blue'

var= [dsh.TEND_O2.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),
      dsh.DIV.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),
      dsh.DIVm.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),
      dsh.DIVp.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),
      dsh.VDIF.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),
      dsh.J_O2.isel(nlon=slice(W,E)).sum('nlon')*nmol_mol/(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon')*cm3_m3),]

cntf= True ; cnt= False 
units=['mol.m$^{-3}$.s$^{-1}$']
title=[r'$\overline{\frac{\partial{O_2}}{\partial{t}}}$ ', 
     r'$-\nabla{\overline{(U.O_2)}}$ ',
     r'$-\nabla{(\overline{U}.\overline{O_2})}$ ',
     r"$-\nabla{(\overline{U'.O_2'})}$",
     r'$\overline{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}$ ',
     r'$\overline{J(O_2)}$ ']

clr=[clr_o2]
lw=0.5; alp=1;cb_on=True

lo = -3e-9; hi = -lo; hi2 = -lo; dsc = -lo*0.5e-1
y=dsh.TLAT[:,E]
z= dsh.z_t*1e-2

i=[0,0,0,1,1,1]; j=[0,1,2,2,1,0]
for k in tqdm(np.arange(0,6)):#len(VAR))):
    ax = fig.add_subplot(gs[i[k],j[k]])
    cax= plot_depth_section(ax,y,z,var[k],cntf,np.arange(lo, hi+dsc, dsc),
                            clr[0],units[0],title[k],coords,fs,fsx,lw,alp,cb_on,'lon')
    pch= cax.contour(y, z,dsh.O2.isel(nlon=slice(W,E)).mean('nlon'),
                     colors='black', alpha=0.6,linewidths=0.3, levels = np.arange(0,300.1,10))
    pch= cax.contour(y, z,dsh.O2.isel(nlon=slice(W,E)).mean('nlon'),
                     colors=om_clr, alpha=1, linewidths=1, levels = np.arange(200,200.1,0.1))
    pch= cax.contour(y, z,dsh.O2.isel(nlon=slice(W,E)).mean('nlon'),
                     colors=omz_clr, alpha=1, linewidths=1, levels = np.arange(60,60.1,0.1))
    pch= cax.contour(y, z,dsh.O2.isel(nlon=slice(W,E)).mean('nlon'),
                     colors=odz_clr, alpha=1, linewidths=1, levels = np.arange(20,20.1,0.1))
    
    
#     cs2= cax.contour(y, z,(dsp.PD.isel(nlon=slice(W,E)).mean('nlon').squeeze()-1)*1e3, colors='orange', linewidths=1.0, levels=np.arange(20,27,0.5))
#     ax2.clabel(cs2,np.arange(20,27,0.5), fontsize=10)


fig.tight_layout


# In[107]:


fig, ax = plt.subplots(1,2, figsize=(20,4), sharey=False)

VARS=[dshb.TEND_O2,
      dshb.DIV,
#       dshb.DIVm,
#       dshb.DIVp,
      dshb.VDIF,
      dshb.J_O2,]
title=[r' $\overline{\frac{\partial{O_2}}{\partial{t}}}$ ', 
     r' $-\nabla{\overline{(U.O_2)}}$ ',
#      r' $-\nabla{(\overline{U}.\overline{O_2})}$ ',
#      r" $-\nabla{(\overline{U'.O_2'})}$",
     r' $\overline{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}$ ',
     r' $\overline{J(O_2)}$ ',]

clr=['grey','black','orange','green'] #'dodgerblue','cyan',

i=0
for VAR in VARS:
    ax[0].plot(dsh.TLAT[:,0],VAR * nmol_mol/VOL_hr,
               label=str(title[i]), color=clr[i],lw=1)
    i+=1
ax[0].legend()
ax[0].set_title('O$_2$ budget in CESM-HR',size=16)# integrated over '+str(abs(Wi))+'W-'+str(abs(Ei))+'W '+str(top)+'-'+str(bot)+'m')
ax[0].set_ylabel('O$_2$ Tendency (mol.m-3.s-1)',fontsize=16)
ax[0].set_xlabel('Latitude (ºN)',fontsize=16)
ax[0].set_xlim(-15,15)
ax[0].set_ylim(-1.2e-8,1.0e-8)


# Plot Advective budget terms for both resolutions     
VARS=[dshb.TEND_O2,
      dshb.DIV,
      dshb.DIVm,
      dshb.DIVp,
      dshb.VDIF,
      dshb.J_O2,]
title=[r' $\overline{\frac{\partial{O_2}}{\partial{t}}}$ ', 
     r' $-\nabla{\overline{(U.O_2)}}$ ',
     r' $-\nabla{(\overline{U}.\overline{O_2})}$ ',
     r" $-\nabla{(\overline{U'.O_2'})}$",
     r' $\overline{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}$ ',
     r' $\overline{J(O_2)}$ ',]

clr=['grey','black','dodgerblue','cyan','orange','green']

i=1
for VAR in VARS[1:4]:
    ax[1].plot(dsh.TLAT[:,0],VAR * nmol_mol/VOL_hr,
               label=str(title[i]), color=clr[i],lw=1.5)
    i+=1
ax[1].legend()
ax[1].set_title('Advective O$_2$ budget in CESM-HR vs CESM-LR',size=16)# integrated over '+str(abs(Wi))+'W-'+str(abs(Ei))+'W '+str(top)+'-'+str(bot)+'m',size=20)
ax[1].set_ylabel('O$_2$ Tendency (mol.m-3.s-1)',fontsize=16)
ax[1].set_xlabel('Latitude (ºN)',fontsize=16)
ax[1].set_xlim(-15,15)
ax[1].set_ylim(-1.0e-8,0.20e-8)


# In[9]:


N,W=find_indices(dsh.TLONG,dsh.TLAT,360+Wi,Ni)
S,E=find_indices(dsh.TLONG,dsh.TLAT,360+Ei,Si)
# volume integral
dshb=(dsh.isel(nlon=slice(W,E)).sum('nlon').isel(nlat=slice(S,N)).sum('nlat').sel(z_t=slice(top*1e2,bot*1e2)).sum(dim='z_t')).compute()
dshb
VOL_hr=(dsh.VOL.isel(nlon=slice(W,E)).sum('nlon').isel(nlat=slice(S,N)).sum('nlat').sel(z_t=slice(top*1e2,bot*1e2)).sum(dim='z_t')*cm3_m3).compute()


# In[10]:


fig = plt.figure(figsize=(8,4),dpi=100,)
VARS=[dshb.TEND_O2,
      dshb.DIV,
      dshb.DIVm,
      dshb.DIVp,
      dshb.HDIF,
      dshb.VDIF,
      dshb.J_O2,]
title=[r' $\overline{\frac{\partial{O_2}}{\partial{t}}}$ ', 
     r' $-\nabla{\overline{(U.O_2)}}$ ',
     r' $-\nabla{(\overline{U}.\overline{O_2})}$ ',
     r" $-\nabla{(\overline{U'.O_2'})}$",
      r' $\overline{A_h.\nabla_H^2{({O_2})}}$',
     r' $\overline{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}$ ',
     r' $\overline{J(O_2)}$ ',]
clr=['grey','blue','dodgerblue','cyan','purple','orange','green']
x=np.arange(0,len(VARS))
plt.bar(x,[V * nmol_mol/VOL_hr for V in VARS],width=0.85,align='center',color=clr,alpha=0.5,lw=1.0,edgecolor='grey')

plt.xticks(x, title,size=11)
plt.hlines(0,-0.5,6.5,lw=0.2,color='k')
plt.grid(b=None)
plt.xlim(-0.5,6.5)
# plt.ylim(-6e5,6e5)
plt.ylabel(r"O$_2$ tendency (mol/s)",size=14)
plt.text(5.5,5.0e5,'CESM-HR',color="dodgerblue",size=18)


# In[99]:


VARS=[dshb.TEND_O2,
      dshb.DIV,
      dshb.DIVm,
      dshb.DIVp,
      dshb.VDIF,
      dshb.J_O2,]
VARC=[r' $\overline{\frac{\partial{O_2}}{\partial{t}}}$ ', 
     r' $-\nabla{\overline{(U.O_2)}}$ ',
     r' $-\nabla{(\overline{U}.\overline{O_2})}$ ',
     r" $-\nabla{(\overline{U'.O_2'})}$",
     r' $\overline{\frac{\partial{}}{\partial{z}}k.\frac{\partial{O_2}}{\partial{z}}}$ ',
     r' $\overline{J(O_2)}$ ',]

clr=['grey','black','dodgerblue','cyan','orange','green']
dpt=500e2
conv_mmol_m3_mo=nmol_mol


i=0
fig, ax = plt.subplots(1,2, figsize=(10,7), sharey=False)
for VAR in VARS:
    (VAR*nmol_mol/VOL_hr).plot(y="z_t", ylim=(dpt, 0), label=str(VARC[i]), ax=ax[0],color=clr[i],lw=1)
    i+=1
ax[0].legend()
ax[0].set_title('O$_2$ budget in CESM-LR integrated over '+str(Si)+'N-'+str(Ni)+'N '+str(abs(Wi))+'W-'+str(abs(Ei))+'W ')
ax[0].set_ylabel('Depth (cm)')
ax[0].set_xlabel('Tendency (mmol/m3/yr)')
ax[0].set_xlim(-3e-9,3e-9)
ax[0].set_ylim(500e2,50e2)

