using DrWatson
@quickactivate "ConvectionIsotopes"
using NASAPrecipitation

npd = IMERGFinalHH(start=Date(2019,7,1),stop=Date(2021,6,30),path=datadir())
geo = GeoRegion("OTREC",path=srcdir())
download(npd,geo)