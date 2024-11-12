using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using GeoRegions
using Printf

include(srcdir("wrfregridtoera5.jl"))

wrfnewregridera5("RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=7)
wrfnewregridera5("RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=30)

wrfnewregridera5("HDO_RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=7)
wrfnewregridera5("HDO_RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=30)

wrfnewregridera5("O18_RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=7)
wrfnewregridera5("O18_RAINNC",start=Date(2019,8,1),stop=Date(2020,12,30),days=30)

wrfoldregridera5("RAINNC",days=7)
wrfoldregridera5("RAINNC",days=30)

wrfoldregridera5("HDO_RAINNC",days=7)
wrfoldregridera5("HDO_RAINNC",days=30)

wrfoldregridera5("O18_RAINNC",days=7)
wrfoldregridera5("O18_RAINNC",days=30)