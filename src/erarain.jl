using Dates
using ERA5Reanalysis
using Printf
using Statistics

function erarain(
    e5ds :: ERA5Reanalysis.ERA5Hourly,
    geo  :: GeoRegion;
)

    dtbegstr = Dates.format(npd.start,dateformat"yyyymmdd")
    dtbegend = Dates.format(npd.stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"

    evar = SingleVariable("tp")
    
    ogeo = GeoRegion("OTREC_wrf_d02",path=srcdir())
    olsd = getLandSea(e5ds,ERA5Region(ogeo))
    ggrd = RegionGrid(geo,olsd.lon,olsd.lat)

    dtvec = e5ds.start : Month(1) : e5ds.stop; ndt = length(dtvec)

    tmp  = zeros(Float32,length(ggrd.lon),length(ggrd.lat),24)
    pvec = zeros(ndt*24)
    ii = 0
	for idt = npd.start : Day(1) : npd.stop
		ds = read(e5ds,evar,ogeo,idt)
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
    fnc = datadir("wrf3","processed","$(geo.ID)-erarain-$timestr.nc")

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt * 24

    nctime = defVar(ds,"time",Float64,("date",),attrib=Dict(
        "units"     => "hours since $(npd.start) 00:00:00.0",
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