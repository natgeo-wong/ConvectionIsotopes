using Dates
using NASAPrecipitation
using Printf
using Statistics

function gpmrain(
    npd  :: IMERGDataset,
    geo  :: GeoRegion;
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    
    olsd = getLandSea(npd,GeoRegion("OTREC",path=srcdir()))
    nlsd = getLandSea(npd,geo)
    ggrd = RegionGrid(geo,olsd.lon,olsd.lat)

    dtvec = npd.start : Day(1) : npd.stop; ndt = length(dtvec) * 24

    tmp1 = zeros(Float32,length(olsd.lon),length(olsd.lat),48)
    tmp2 = zeros(Float32,length(nlsd.lon),length(nlsd.lat),48)
    pvec = zeros(ndt*48)
    ii = 0
	for idt = npd.start : Day(1) : npd.stop
		ds = read(npd,GeoRegion("OTREC",path=srcdir()),idt)
		NCDatasets.load!(ds["precipitation"].var,tmp,:,:,:)
		close(ds)

        extract!(tmp2,tmp1,ggrd)

        for ihr = 1 : 48
            ii += 1
            iprcp = @views tmp2[:,:,ihr]
            pvec[ii] = mean(iprcp[.!isnan.(iprcp)])
        end
	end

    pvec = dropdims(mean(reshape(pvec,2,:),dims=1),dims=1)

    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","$(geo.ID)-gpmrain-$timestr.nc")

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt * 24

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncrain = defVar(ds,"precipitation",Float32,("date",),attrib=Dict(
        "units"     => "mm day**-1",
    ))

    nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
    ncrain[:] = pvec

    close(ds)

end