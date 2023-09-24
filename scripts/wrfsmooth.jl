using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfsmooth.jl"))

wrf3Dsmooth("QVAPOR",GeoRegion("OTREC_SAN"))