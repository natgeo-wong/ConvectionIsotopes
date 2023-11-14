using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using Distributed
using Printf

addprocs(15)

@everywhere include(srcdir("wrfbudget.jl"))

@distributed for i in 1 : (12 * 16)
    boxii = @sprintf("%02d",mod(i,16)+1)
    stnii = @sprintf("%02d",ceil(i/16))
    stnstr = "OTREC_STN$(stnii)_$(boxii)"
    wrfqdivdecompose(GeoRegion(stnstr),start=Date(2019,8),stop=Date(2019,12,30))
end