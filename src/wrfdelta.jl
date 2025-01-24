using Dates
using GeoRegions
using NCDatasets
using Printf
using RegionGrids
using Statistics
using Trapz

include(srcdir("backend.jl"))

function wrfdelta(
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

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    ndt  = length(dtvec)

    tmpRAIN = zeros(Float32,nlon,nlat)
    tmpHDO  = zeros(Float32,nlon,nlat)
    tmpO18  = zeros(Float32,nlon,nlat)

    rain = zeros(Float32,ndt)
    HDOr = zeros(Float32,ndt)
    O18r = zeros(Float32,ndt)

    if iszero(days)
        ds1 = NCDataset(datadir("wrf3","2D","RAINNC-daily-$timestr.nc"))
        ds2 = NCDataset(datadir("wrf3","2D","HDO_RAINNC-daily-$timestr.nc"))
        ds3 = NCDataset(datadir("wrf3","2D","O18_RAINNC-daily-$timestr.nc"))
    else
        ds1 = NCDataset(datadir("wrf3","2D","RAINNC-daily-$timestr-$smthstr.nc"))
        ds2 = NCDataset(datadir("wrf3","2D","HDO_RAINNC-daily-$timestr-$smthstr.nc"))
        ds3 = NCDataset(datadir("wrf3","2D","O18_RAINNC-daily-$timestr-$smthstr.nc"))
    end

    for ii in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for Day $ii of $ndt"
        flush(stderr)

        NCDatasets.load!(ds1["RAINNC"].var,tmpRAIN,lon1:lon2,lat1:lat2,ii)
        NCDatasets.load!(ds2["HDO_RAINNC"].var,tmpHDO,lon1:lon2,lat1:lat2,ii)
        NCDatasets.load!(ds3["O18_RAINNC"].var,tmpO18,lon1:lon2,lat1:lat2,ii)

        rain[ii] = mean(tmpRAIN)
        HDOr[ii] = mean(tmpHDO)
        O18r[ii] = mean(tmpO18)
        

    end

    close(ds1)
    close(ds2)
    close(ds3)

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-rain-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-rain-daily-$timestr-$smthstr.nc")
    end

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncrain = defVar(ds,"RAINNC",Float32,("date",),attrib=Dict(
        # "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        # "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "mm day**-1",
    ))

    ncHDO = defVar(ds,"HDO_RAINNC",Float32,("date",),attrib=Dict(
        # "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        # "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "mm day**-1",
    ))

    ncO18 = defVar(ds,"O18_RAINNC",Float32,("date",),attrib=Dict(
        # "long_name" => "pressure",
        # "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "mm day**-1",
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncrain[:] = rain
    ncHDO[:]  = HDOr
    ncO18[:]  = O18r

    close(ds)

end