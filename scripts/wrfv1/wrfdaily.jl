using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfdaily.jl"))

wrf2Ddaily("RAINNC",start=Date(2019,8),stop=Date(2019,12,31))