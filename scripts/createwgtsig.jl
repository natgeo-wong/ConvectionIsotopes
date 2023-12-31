using DrWatson
@quickactivate "ConvectionIsotopes"
using ERA5Reanalysis

include(srcdir("wwgtpre_create.jl"))

addGeoRegions(srcdir("OTRECGeoRegions.txt"))
addERA5Variables(srcdir("templatesc.txt"))

e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,1,31),path="/n/kuangdss01/lab")

create_wσ(e5ds,ERA5Region(GeoRegion("OTREC")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_SAN")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_CLB")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_MNT")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_INL")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_PAC")))
create_wσ(e5ds,ERA5Region(GeoRegion("CIS_CRB")))