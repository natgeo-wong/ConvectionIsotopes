using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using GeoRegions
using Printf

include(srcdir("wrfdqdp.jl"))

for iiso = ["HDO","O18"], istn = 1 : 12
    stngeo = GeoRegion("OTREC_STN$(@sprintf("%02d",istn))")
    wrfdqdp(stngeo,iso=iiso,start=Date(2019,8),stop=Date(2019,12,30),days=7)
    wrfdqdp(stngeo,iso=iiso,start=Date(2019,8),stop=Date(2019,12,30),days=30)
    for ibox = 1 : 4
        boxgeo = GeoRegion("OTREC_STN$(@sprintf("%02d",istn))_$(@sprintf("%02d",ibox))")
        wrfdqdp(boxgeo,iso=iiso,start=Date(2019,8),stop=Date(2019,12,30),days=7)
        wrfdqdp(boxgeo,iso=iiso,start=Date(2019,8),stop=Date(2019,12,30),days=30)
    end
end