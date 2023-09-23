using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfbudget.jl"))

wrfqbudget("QVAPOR",GeoRegion("OTREC_SAN"))