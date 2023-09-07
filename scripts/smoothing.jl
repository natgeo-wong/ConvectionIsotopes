using DrWatson
@quickactivate "ColombiaIsotope"

using ERA5Reanalysis

e5ds = ERA5Daily(start=Date(2019,7),stop=Date(2021,6),path=datadir())
egeo = ERA5Region("OTREC",resolution=0.25)
epre = era5Pressures(); epre = epre[epre.>=10]

for pre in epre
    evar = PressureVariable("w",hPa=pre)
    smoothing(e5ds,evar,egeo,days=30)
end