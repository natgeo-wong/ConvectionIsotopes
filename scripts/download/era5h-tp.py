#!/usr/bin/env python
import cdsapi
import os

datadir = '/n/holyscratch01/kuang_lab/nwong/ColumbiaIsotope/data/reanalysis/tp/'

c = cdsapi.Client()

if not os.path.exists(datadir):
    os.makedirs(datadir)

for mo in range(1,13):
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'format': 'netcdf',
            'product_type': 'reanalysis',
            'variable': 'total_precipitation',
            'year': 2020,
            'month': mo,
            'day': [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                22, 23, 24, 25, 26, 27, 28, 29, 30, 31
            ],
            'time': [
                '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00', '21:00', '22:00', '23:00',
            ],
            'area': [30.0, -100.0, -15.0, -55.0],
        },
        datadir + 'era5h-OTRECx0.25-tp-' + str(2020) + str(mo).zfill(2) + '.nc')

for mo in range(1,7):
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'format': 'netcdf',
            'product_type': 'reanalysis',
            'variable': 'total_precipitation',
            'year': 2021,
            'month': mo,
            'day': [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                22, 23, 24, 25, 26, 27, 28, 29, 30, 31
            ],
            'time': [
                '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00', '21:00', '22:00', '23:00',
            ],
            'area': [30.0, -100.0, -15.0, -55.0],
        },
        datadir + 'era5h-OTRECx0.25-tp-' + str(2021) + str(mo).zfill(2) + '.nc')
