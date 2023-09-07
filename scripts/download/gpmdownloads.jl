using DrWatson
@quickactivate "ColombiaIsotope"
using NASAPrecipitation

addGeoRegions(srcdir("OTRECGeoRegions.txt"))

npd = IMERGFinalHH(dtbeg=Date(2019,7,1),dtend=Date(2021,6,30),sroot=datadir())
geo = GeoRegion("OTREC")
download(npd,geo)