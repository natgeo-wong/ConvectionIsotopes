#!/usr/bin/env python
import cdsapi
import os

datadir = '/n/holyscratch01/kuang_lab/nwong/ColumbiaIsotope/data/reanalysis/'

c = cdsapi.Client()

if not os.path.exists(datadir):
    os.makedirs(datadir)

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'vertical_velocity',
        'pressure_level': 500,
        'year': [
            2001, 2002, 2003, 2004, 2005,
            2006, 2007, 2008, 2009, 2010,
            2011, 2012, 2013, 2014, 2015,
            2016, 2017, 2018, 2019, 2020,
        ],
        'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        'area': [30.0, -100.0, -15.0, -55.0],
        'time': '00:00',
    },
    datadir + 'era5m-OTRECx0.25-w-500hPa.nc')
