using DrWatson
@quickactivate "ConvectionIsotopes"
using ERA5Reanalysis

addGeoRegions(srcdir("OTRECGeoRegions.txt"))

e5ds = ERA5Monthly(dtbeg=Date(2019,1,1),dtend=Date(2021,12,31),eroot="/n/kuangdss01/lab")
ereg = ERA5Region(GeoRegion("OTREC"))
evar_sp = SingleVariable("sp")
download(e5ds,evar_sp,ereg,ispy=true)

for lvl in era5Pressures()
    evar_w = PressureVariable("w",hPa=lvl)
    download(e5ds,evar_w,ereg,ispy=true)
end

getLandSea(e5ds,ereg,returnlsd=false)
