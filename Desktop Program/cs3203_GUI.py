# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 17:21:01 2023

@author: James
"""

# Import libraries
import getopt, sys
import os.path
import subprocess
import xml.etree.ElementTree as ET
import PySimpleGUI as psg

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

def calcAverages():
    
    temp_df = Dataset('air.2m.mon.mean.nc', 'r')
    solar_df = Dataset('dswrf.sfc.mon.mean.nc', 'r')
    longwave_df = Dataset('ulwrf.sfc.mon.mean.nc', 'r')
    precip_df = Dataset('prate.sfc.mon.mean.nc', 'r')
    shum_df = Dataset('shum.2m.mon.mean.nc', 'r')
    
    start = datetime.datetime(start_year, startMonth_idx, 1, 0, 0)
    stop = datetime.datetime(end_year, endMonth_idx, 1, 0, 0)
    
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
    
    temp_str = "Surface temperature: " + "{:.2f}".format(np.nanmean(temp_mean)) + " degF"
    solarRad_str = "Downward solar radiation flux: " + "{:.2f}".format(np.nanmean(solar_mean)) + " W/m^2"
    longwaveRad_str = "Upward longwave radiation flux: " + "{:.2f}".format(np.nanmean(longwave_mean)) + " W/m^2"
    precip_str = "Precipitation rate: " + "{:.2f}".format(np.nanmean(precip_mean)) + " mm/day"
    shum_str = "Specific humidity: " + "{:.2f}".format(np.nanmean(shum_mean)) + " g/kg"
    
    longitude = temp_df.variables['lon'][:]
    latitude = temp_df.variables['lat'][:]
    
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
    cbar.set_label('Climatological Average Temperature (\N{DEGREE FAHRENHEIT})', size = 'x-large')
    cbar.solids.set_edgecolor('face')
    
    plt.title('Average Temperature (Global Mean = ' + '{:.2f}'.format(np.nanmean(temp_mean)) + ' \N{DEGREE FAHRENHEIT})', size = 'xx-large')
    plt.tight_layout()
    plt.savefig('Global Avg Temp.png', dpi = 50)
    
    # Plot for average precipitation
    
    plot_precip, longitude2 = add_cyclic_point(precip_mean, longitude)
    
    precip_fig = plt.figure(figsize = (12, 8))
    precip_ax = plt.subplot(111, projection = projection)
    precip_ax.add_feature(cfeature.COASTLINE, linewidth = 1, zorder = 1)
    precip_ax.add_feature(cfeature.BORDERS, linewidth = 0.5, zorder = 1)
    
    cs = precip_ax.contourf(plotLon, plotLat, plot_precip, np.arange(0, 27.5, 0.5), cmap = colormap, transform = data_crs, zorder = 0)
    cbar = precip_fig.colorbar(cs, orientation = 'horizontal', extend = 'max', aspect = 40, shrink = 0.8, pad = 0.05, extendrect = 'True')
    cbar.set_label('Climatological Average Precipitation Rate (mm/day)', size = 'x-large')
    cbar.solids.set_edgecolor('face')
    
    plt.title('Average Precipitation Rate (Global Mean = ' + '{:.2f}'.format(np.nanmean(precip_mean)) + ' mm/day)', size = 'xx-large')
    plt.tight_layout()
    plt.savefig('Global Avg Precip Rate.png', dpi = 50)
    
    layout = [[psg.Text('Average Global Climatological Variables', text_color = 'white', font = ('Lucinda', 18, 'bold'), justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(start_month + ' 1, ' + str(start_year) + ' - ' + end_month + ' 1, ' + str(end_year), font = ('Lucinda', 16, 'bold'), justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(temp_str, font = 'Lucinda', justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(solarRad_str, font = 'Lucinda', justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(longwaveRad_str, font = 'Lucinda', justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(precip_str, font = 'Lucinda', justification = 'center', expand_x = True, expand_y = True)],
              [psg.Text(shum_str, font = 'Lucinda', justification = 'center', expand_x = True, expand_y = True)],
              [psg.Image('Global Avg Temp.png', expand_x = True, expand_y = True), psg.Image('Global Avg Precip Rate.png', expand_x = True, expand_y = True)],
              [psg.Push(), psg.Button('Exit', font = ('Times New Roman', 12, 'bold'), key = 'exit', expand_x = False, expand_y = False)]]
    
    win = psg.Window('Results', layout, modal = True)
    
    while True:
        e,v = win.read()
        if e in ('exit', None):
            win.close()
            break

#################################################################################################################################################################################################################################################################################

debugflag = True

month_options = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
year_options = [1948, 1949, 1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970,
                1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993,
                1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
                2017, 2018, 2019, 2020, 2021, 2022, 2023]

# Set up the GUI
psg.theme('DarkBlue8')

layout = [[psg.Text('Username: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.InputText(size = (10, 1), font = 'Lucinda', key = 'username', expand_x = True, expand_y = True)],
          [psg.Text('Password: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.InputText(size = (10, 1), font = 'Lucinda', key = 'password', enable_events = True, expand_x = True, expand_y = True)],
          [psg.Text(size = (30, 1), font = 'Lucinda', justification = 'left', key = 'greeting', expand_x = True, expand_y = True)],
          [psg.Text()],
          [psg.Text('Climate Outlook Tool', text_color = 'white', font = ('Lucinda', 18, 'bold'), justification = 'center', expand_x = True, expand_y = True)],
          [psg.Text(font = ('Lucinda', 12, 'bold'), justification = 'center', key = 'prompt', expand_x = True, expand_y = True)],
          [psg.Text('Enter start year: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.Spin(year_options, initial_value = 1948, readonly = False, size = (9, 1), disabled = True, key = 'Start_Year', expand_x = True, expand_y = True)],
          [psg.Text('Enter start month: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.Spin(month_options, initial_value = 'January', readonly = False, size = (9, 1), disabled = True, key = 'Start_Month', expand_x = True, expand_y = True)],
          [psg.Text('Enter end year: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.Spin(year_options, initial_value = 1948, readonly = False, size = (9, 1), disabled = True, key = 'End_Year', expand_x = True, expand_y = True)],
          [psg.Text('Enter end month: ', size = (20, 1), font = 'Lucinda', justification = 'left'), psg.Spin(month_options, initial_value = 'January', readonly = False, size = (9, 1), disabled = True, key = 'End_Month', expand_x = True, expand_y = True)],
          [psg.Button('Calculate Averages', font = ('Times New Roman', 12, 'bold'), disabled = True, key = 'calcBtn', expand_x = True, expand_y = True), psg.Push(), psg.Button('Exit', font = ('Times New Roman', 12, 'bold'), key = 'exit', expand_x = True, expand_y = True)]]

win = psg.Window('Climate Outlook', layout)

# Make some functions
def buttonDisable():
    win['Start_Year'].update(disabled = True)
    win['Start_Month'].update(disabled = True)
    win['End_Year'].update(disabled = True)
    win['End_Month'].update(disabled = True)
    win['calcBtn'].update(disabled = True)
    win.Refresh()

def buttonEnable():
    win['Start_Year'].update(disabled = False)
    win['Start_Month'].update(disabled = False)
    win['End_Year'].update(disabled = False)
    win['End_Month'].update(disabled = False)
    win['calcBtn'].update(disabled = False)
    win.Refresh()

# Now run the working loop
while True:
    e, v = win.read()
    if e in ('exit', None):
        win.close()
        break
    if e == 'password':
        if v['password'] == '4774' and v['username'] != '':
            buttonEnable()
            win['greeting'].update('Welcome to Climate Outlook, ' + v['username'] + '!')
            win['prompt'].update('Please enter a range of dates below!')
            win['username'].update('')
            win['password'].update('')
        elif v['password'] == '4774' and v['username'] == '':
            win['password'].update('')
            win['greeting'].update('Enter your username and try again!')
    if e == 'calcBtn':
        buttonDisable()
        if debugflag:
            print('Configure Pressed')
        start_month = v['Start_Month']
        startMonth_idx = month_options.index(v['Start_Month']) + 1
        start_year = v['Start_Year']
        end_month = v['End_Month']
        endMonth_idx = month_options.index(v['End_Month']) + 1
        end_year = v['End_Year']
        calcAverages()
        if debugflag:
            print("config is: ", 'config')
        buttonEnable()

win.close()