using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using Printf

include(srcdir("wrfwwgtpre.jl"))

for istn = 1 : 12

    geoID = "OTREC_STN$(@sprintf("%02d",istn))"
    wrfwwgtpre(GeoRegion(geoID),start=Date(2019,8),stop=Date(2019,12,30),days=7)
    wrfwwgtpre(GeoRegion(geoID),start=Date(2019,8),stop=Date(2019,12,30),days=30)

    for ibox = 1 : 4
        boxID = "OTREC_STN$(@sprintf("%02d",istn))_$(@sprintf("%02d",ibox))"
        wrfwwgtpre(GeoRegion(boxID),start=Date(2019,8),stop=Date(2019,12,30),days=7)
        wrfwwgtpre(GeoRegion(boxID),start=Date(2019,8),stop=Date(2019,12,30),days=30)
    end

end