using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfbudget.jl"))

wrfqdiv("QVAPOR",GeoRegion("OTREC_SAN"),start=Date(2019,8),stop=Date(2019,12,30))