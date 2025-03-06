using Dates
using GeoRegions
using RegionGrids
using NCDatasets
using Printf
using Statistics
using Trapz

include(srcdir("backend.jl"))

function wrfwwgtpre2(
    geo  :: GeoRegion;
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    close(ds)

    ggrd = RegionGrid(geo,Point2.(lon,lat))
    lon1 = minimum(ggrd.ilon); lon2 = maximum(ggrd.ilon)
    lat1 = minimum(ggrd.ilat); lat2 = maximum(ggrd.ilat)

    if iszero(days)
        ds = NCDataset(datadir("wrf3","2D","p_wwgt-daily-$timestr.nc"))
    else
        ds = NCDataset(datadir("wrf3","2D","p_wwgt-daily-$timestr-$smthstr.nc"))
    end

    pw = ds["p_wwgt"][lon1:lon2,lat1:lat2,:]; ndt = size(pw,3)
    σw = ds["σ_wwgt"][lon1:lon2,lat1:lat2,:]

    pw2 = zeros(ndt)
    σw2 = zeros(ndt)

    for it = 1 : ndt
        ipw = @views pw[:,:,it]
        iσw = @views σw[:,:,it]
        pw2[it] = mean(.!isnan.(ipw))
        σw2[it] = mean(.!isnan.(iσw))
    end

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt2-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt2-daily-$timestr-$smthstr.nc")
    end

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncpwgt = defVar(ds,"p_wwgt",Float32,("date",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "Pa",
    ))

    ncσwgt = defVar(ds,"σ_wwgt",Float32,("date",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "0-1",
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpwgt[:] = pw2
    ncσwgt[:] = σw2

    close(ds)

end