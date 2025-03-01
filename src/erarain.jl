using Dates
using ERA5Reanalysis
using Printf
using Statistics

function erarain(
    e5ds :: ERA5Reanalysis.ERA5Hourly,
    ngeo :: ERA5Region;
)

    dtbegstr = Dates.format(e5ds.start,dateformat"yyyymmdd")
    dtbegend = Dates.format(e5ds.stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"

    evar = SingleVariable("tp")
    
    ogeo = GeoRegion("OTREC_wrf_d02",path=srcdir())
    egeo = ERA5Region(ogeo)
    olsd = getLandSea(e5ds,egeo)
    ggrd = RegionGrid(ngeo,olsd.lon,olsd.lat)

    dtvec = e5ds.start : Month(1) : e5ds.stop; ndt = length(dtvec)

    tmp  = zeros(Float32,length(ggrd.lon),length(ggrd.lat),24)
    pvec = zeros(ndt*24)
    ii = 0
	for idt = e5ds.start : Day(1) : e5ds.stop
		ds = read(e5ds,evar,egeo,idt)
		tprcp = nomissing(ds[evar.ID][:,:,:])
		close(ds)
        nt = size(tprcp,3)

        extract!(tmp,tprcp,ggrd)

        for ihr = 1 : nt
            ii += 1
            iprcp = @views tmp2[:,:,ihr]
            pvec[ii] = mean(iprcp[.!isnan.(iprcp)])
        end
	end

    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","$(ngeo.geo.ID)-erarain-$timestr.nc")

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt * 24

    nctime = defVar(ds,"time",Float64,("date",),attrib=Dict(
        "units"     => "hours since $(e5ds.start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncrain = defVar(ds,"precipitation",Float32,("date",),attrib=Dict(
        "units"     => "mm day**-1",
    ))

    nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
    ncrain[:] = pvec * 1000

    close(ds)

end