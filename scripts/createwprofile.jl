using DrWatson
@quickactivate "ColombiaIsotope"
using ERA5Reanalysis

include(srcdir("wwgtpre_profile.jl"))

addGeoRegions(srcdir("OTRECGeoRegions.txt"))
addERA5Variables(srcdir("templatesc.txt"))

e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,1,31),path="/n/kuangdss01/lab")

wprofile(e5ds,ERA5Region(GeoRegion("OTREC")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_SAN")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_CLB")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_MNT")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_INL")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_PAC")))
wprofile(e5ds,ERA5Region(GeoRegion("CIS_CRB")))