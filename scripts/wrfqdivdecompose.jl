using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using Distributed
using Printf

include(srcdir("wrfbudget.jl"))

for istn = 1 : 12, ibox = 1 : 16
    stnstr = "OTREC_STN$(@sprintf("%02d",istn))_$(@sprintf("%02d",ibox))"
    wrfqdivdecompose(GeoRegion(stnstr),start=Date(2019,8),stop=Date(2019,12,30))
end