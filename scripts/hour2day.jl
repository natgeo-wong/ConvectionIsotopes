using DrWatson
@quickactivate "ConvectionIsotopes"

using ERA5Reanalysis

e5ds = ERA5Hourly(start=Date(2019,1),stop=Date(2021,12),path=datadir())
egeo = ERA5Region("OTREC",resolution=0.25)
epre = era5Pressures(); epre = epre[epre.>=10]

for pre in epre
    evar = PressureVariable("w",hPa=pre)
    hourly2daily(e5ds,evar,egeo)
end