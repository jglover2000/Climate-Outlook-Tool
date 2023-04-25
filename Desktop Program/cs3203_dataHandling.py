# -*- coding: utf-8 -*-
"""
Created on Sat Apr 08 18:58:36 2023

@author: James
"""

from netCDF4 import Dataset, num2date, MFDataset
from cartopy.util import add_cyclic_point
from scipy.ndimage import gaussian_filter
from metpy.units import units

import netCDF4
import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import metpy.calc as mpcalc
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

temp_df = Dataset('air.2m.mon.mean.nc', 'r')
solar_df = Dataset('dswrf.sfc.mon.mean.nc', 'r')
longwave_df = Dataset('ulwrf.sfc.mon.mean.nc', 'r')
precip_df = Dataset('prate.sfc.mon.mean.nc', 'r')
shum_df = Dataset('shum.2m.mon.mean.nc', 'r')

start_year = int(input("Enter start year (between 1948 and 2023): "))
start_month = int(input("Enter start month as an integer (January = 1, February = 2, etc.): "))
end_year = int(input("Enter end year (between 1948 and 2023): "))
end_month = int(input("Enter end month as an integer (January = 1, February = 2, etc.): "))

start = datetime.datetime(start_year, start_month, 1, 0, 0)
stop = datetime.datetime(end_year, end_month, 1, 0, 0)

initial_index = netCDF4.date2index(start, temp_df.variables['time'], select = 'nearest')
final_index = netCDF4.date2index(stop, temp_df.variables['time'], select = 'nearest')

temperature = temp_df.variables['air'][initial_index:final_index, :, :]
solar_rad = solar_df.variables['dswrf'][initial_index:final_index, :, :]
longwave_rad = longwave_df.variables['ulwrf'][initial_index:final_index, :, :]
precipitation = precip_df.variables['prate'][initial_index:final_index, :, :]
specific_hum = shum_df.variables['shum'][initial_index:final_index, :, :]

temperature = np.where(temperature.mask, np.nan, temperature.data)
solar_rad = np.where(solar_rad.mask, np.nan, solar_rad.data)
longwave_rad = np.where(longwave_rad.mask, np.nan, longwave_rad.data)
precipitation = np.where(precipitation.mask, np.nan, precipitation.data)
specific_hum = np.where(specific_hum.mask, np.nan, specific_hum.data)

temp_mean = np.nanmean(temperature, axis = 0)
solar_mean = np.nanmean(solar_rad, axis = 0)
longwave_mean = np.nanmean(longwave_rad, axis = 0)
precip_mean = np.nanmean(precipitation, axis = 0)
shum_mean = np.nanmean(specific_hum, axis = 0)

temp_mean = (temp_mean - 273.15) * (9/5) + 32
precip_mean = precip_mean * 60 * 60 * 24

print("Average global temperature: " + "{:.2f}".format(np.nanmean(temp_mean)) + " \N{DEGREE FAHRENHEIT}")
print("Average global solar radiation flux (downward): " + "{:.2f}".format(np.nanmean(solar_mean)) + " W/m^2")
print("Average global longwave radiation flux (upward): " + "{:.2f}".format(np.nanmean(longwave_mean)) + " W/m^2")
print("Average global precipitation rate: " + "{:.2f}".format(np.nanmean(precip_mean)) + " mm/day")
print("Average global specific humidity: " + "{:.2f}".format(np.nanmean(shum_mean)) + " g/kg")

longitude = temp_df.variables['lon'][:]
latitude = temp_df.variables['lat'][:]

cmin = 0.0
cmax = 8.1
cint = 0.5
clevs = np.round(np.arange(cmin, cmax, cint), 1)
colormap = plt.cm.Spectral_r

# Plot for average temperature

plot_temp, longitude2 = add_cyclic_point(temp_mean, longitude)
plotLon, plotLat = np.meshgrid(longitude2, latitude)

projection = ccrs.Robinson(central_longitude = 180)
data_crs = ccrs.PlateCarree()

temp_fig = plt.figure(figsize = (12, 8))
temp_ax = plt.subplot(111, projection = projection)
temp_ax.add_feature(cfeature.COASTLINE, linewidth = 1, zorder = 1)
temp_ax.add_feature(cfeature.BORDERS, linewidth = 0.5, zorder = 1)

cs = temp_ax.contourf(plotLon, plotLat, plot_temp, np.arange(-105.5, 105.5, 0.5), cmap = colormap, transform = data_crs, zorder = 0)
cbar = temp_fig.colorbar(cs, orientation = 'horizontal', extend = 'max', aspect = 40, shrink = 0.8, pad = 0.05, extendrect = 'True')
cbar.set_label('Climatological Average Temperature (\N{DEGREE FAHRENHEIT})', size = 'large')
cbar.solids.set_edgecolor('face')

plt.title('Average Temperature (Global Mean = ' + '{:.2f}'.format(np.nanmean(temp_mean)) + ' \N{DEGREE FAHRENHEIT})', size = 'x-large')
plt.tight_layout()
plt.savefig('Global Avg Temp.png', dpi = 500)

# Plot for average precipitation

plot_precip, longitude2 = add_cyclic_point(precip_mean, longitude)

precip_fig = plt.figure(figsize = (12, 8))
precip_ax = plt.subplot(111, projection = projection)
precip_ax.add_feature(cfeature.COASTLINE, linewidth = 1, zorder = 1)
precip_ax.add_feature(cfeature.BORDERS, linewidth = 0.5, zorder = 1)

cs = precip_ax.contourf(plotLon, plotLat, plot_precip, np.arange(0, 27.5, 0.5), cmap = colormap, transform = data_crs, zorder = 0)
cbar = precip_fig.colorbar(cs, orientation = 'horizontal', extend = 'max', aspect = 40, shrink = 0.8, pad = 0.05, extendrect = 'True')
cbar.set_label('Climatological Average Precipitation Rate (mm/day)', size = 'large')
cbar.solids.set_edgecolor('face')

plt.title('Average Precipitation Rate (Global Mean = ' + '{:.2f}'.format(np.nanmean(precip_mean)) + ' mm/day)', size = 'x-large')
plt.tight_layout()
plt.savefig('Global Avg Precip Rate.png', dpi = 500)