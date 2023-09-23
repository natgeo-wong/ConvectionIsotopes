using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfbudget.jl"))

wrfbudget("QVAPOR",GeoRegion("OTREC_SAN"))