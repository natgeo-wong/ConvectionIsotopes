using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates

include(srcdir("wrfsmooth.jl"))

wrfwwgtpre(GeoRegion("OTREC_STN01"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN02"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN03"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN04"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN05"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN06"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN07"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN08"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN09"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN10"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN11"),smooth=true,smoothtime=7)
wrfwwgtpre(GeoRegion("OTREC_STN12"),smooth=true,smoothtime=7)

wrfwwgtpre(GeoRegion("OTREC_STN01"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN02"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN03"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN04"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN05"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN06"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN07"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN08"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN09"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN10"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN11"),smooth=true,smoothtime=30)
wrfwwgtpre(GeoRegion("OTREC_STN12"),smooth=true,smoothtime=30)