using DrWatson
@quickactivate "ConvectionIsotopes"
using ERA5Reanalysis

addGeoRegions(srcdir("OTRECGeoRegions.txt"))

e5ds = ERA5Hourly(start=Date(2019,1,1),stop=Date(2021,12,31),path=datadir())
ereg = ERA5Region("OTREC")
evar_sp = SingleVariable("sp"); download(e5ds,evar_sp,ereg,overwrite=true)
evar_w  = PressureVariable("w",hPa=0); download(e5ds,evar_w,ereg,pall=true,ptop=10,overwrite=true)
