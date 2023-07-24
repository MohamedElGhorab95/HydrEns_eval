# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:06:29 2023

@author: M Elghorab
"""

from engine.fr_entities_tools import *
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
rad = Observation("C:/netCDFs/fertig/radRW_ICO.nc").gen_observation_field()

mu = rad.extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec().average
ma = rad.extract_by_shp("shp/Mandau/egp6741491_Zittau_5_TEZG_neu.shp").aggr_temporal(3).avg_areal_prec().average
we = rad.extract_by_shp("shp/WeiseElster/OelsnitzTEZG_DHDN.shp").aggr_temporal(3).avg_areal_prec().average

ado = rad.extract_by_shp("shp/Adorf/Adorf.shp").aggr_temporal(3).avg_areal_prec().average
bad = rad.extract_by_shp("shp/BadElster/BadElster.shp").aggr_temporal(3).avg_areal_prec().average
geis = rad.extract_by_shp("shp/Geising/Geising.shp").aggr_temporal(3).avg_areal_prec().average
gros = rad.extract_by_shp("shp/Grossschoenau/Grossschoenau.shp").aggr_temporal(3).avg_areal_prec().average
lau  = rad.extract_by_shp("shp/Lauenstein4/Lauenstein4.shp").aggr_temporal(3).avg_areal_prec().average
nied = rad.extract_by_shp("shp/Niederoderwitz/Niederoderwitz.shp").aggr_temporal(3).avg_areal_prec().average
seif = rad.extract_by_shp("shp/Seifennersdorf/Seifennersdorf.shp").aggr_temporal(3).avg_areal_prec().average

# mu.plot()

# ma.plot()

# we.plot()


tot = xr.concat([mu, ma, we, ado, bad, geis, gros, lau, nied, seif ], dim='time')

# tot.plot()

vals = [i  for i in tot.values if i >0]

plt.plot([i for i in range(len(vals))], vals)
plt.plot([0,len(vals)],[3.5,3.5])

# def get_exceeded_value(values):
#     sorted_values = sorted(values)
#     threshold_index = int(len(sorted_values) *(1- 0.1))
#     exceeded_value = sorted_values[threshold_index]
#     return exceeded_value

# get_exceeded_value(vals)


sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(dpi=750)

# plotting the results
plt.plot([i for i in range(len(vals))], vals,'o', mfc='none',label='Concatinated Events' )
plt.plot([0,len(vals)], [2.5,2.5], 'r--',linewidth=1.75,label='90% Quantile' )

ax.legend()
ax.set_xlim(0,3000)
ax.set_ylim(0,20)
# ax.set_title("Skill score of ICOND2 forecast performance \n{} catchment")
ax.set_xlabel('Event instances')
ax.set_ylabel('Rainfall depth (mm/3hrs)')





