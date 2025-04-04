using GeoRegions
using RegionGrids
using NCDatasets
using StatsBase

include(srcdir("backend.jl"))

function wrfdqdp(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)
    
    ds  = NCDataset(datadir("wrf3","grid.nc"))
    lon = ds["longitude"][:,:]
    lat = ds["latitude"][:,:]
    pbs = ds["pressure_base"][:,:,:]
    close(ds)

    ggrd = RegionGrid(geo,Point2.(lon,lat))
    lon1 = minimum(ggrd.ilon); lon2 = maximum(ggrd.ilon)
    lat1 = minimum(ggrd.ilat); lat2 = maximum(ggrd.ilat)

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtm = sum(wgts)

    if iso != ""; iso = "$(iso)_" end
    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    tmpq = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmpp = zeros(Float32,nlon,nlat,nlvl,ndt)

    pvec  = zeros(Float32,nlvl,ndt)
    dqdp  = zeros(Float32,nlvl,ndt)

    pbs = dropdims(sum(pbs[lon1:lon2,lat1:lat2,:] .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm

    if iszero(days)
        dsq = NCDataset(datadir("wrf3","3D","$(iso)QVAPOR-daily-$timestr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
    else
        dsq = NCDataset(datadir("wrf3","3D","$(iso)QVAPOR-daily-$timestr-$smthstr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
    end
    NCDatasets.load!(dsq["$(iso)QVAPOR"].var,tmpq,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dsp["P"].var,tmpp,lon1:lon2,lat1:lat2,:,:)

    close(dsq)
    close(dsp)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)
        
        iiq = @view tmpq[:,:,:,idt]
        iip = @view tmpp[:,:,:,idt]

        q = dropdims(sum(iiq .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
        p = dropdims(sum(iip .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm .+ pbs

        for ilvl = 2 : (nlvl-1)
            dqdp[ilvl,idt] = (q[ilvl+1] - q[ilvl-1]) / (p[ilvl+1] - p[ilvl-1])
        end
        dqdp[1,idt] = (q[2] - q[1]) / (p[2] - p[1])
        dqdp[nlvl,idt] = q[nlvl-1] / p[nlvl-1]
        pvec[:,idt] = p

    end

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)dqdp-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)dqdp-daily-$timestr-$smthstr.nc")
    end
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["level"] = nlvl
    ds.dim["date"]  = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncpres = defVar(ds,"P",Float32,("level","date",),attrib=Dict(
        "units" => "Pa",
        "long_name" => "Pressure"
    ))

    ncdqdp = defVar(ds,"$(iso)dqdp",Float32,("level","date",),attrib=Dict(
        "long_name" => "Gradient of $(iso)QVAPOR (relative to SMOW) against pressure"
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpres[:,:]   = pvec
    ncdqdp[:,:]   = dqdp

    close(ds)

end

function wrfdhqdp(
    geo   :: GeoRegion;
    iso   :: AbstractString,
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)
    
    ds  = NCDataset(datadir("wrf3","grid.nc"))
    lon = ds["longitude"][:,:]
    lat = ds["latitude"][:,:]
    pbs = ds["pressure_base"][:,:,:]
    close(ds)

    ggrd = RegionGrid(geo,Point2.(lon,lat))
    lon1 = minimum(ggrd.ilon); lon2 = maximum(ggrd.ilon)
    lat1 = minimum(ggrd.ilat); lat2 = maximum(ggrd.ilat)

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtm = sum(wgts)

    if iso != ""; iso = "$(iso)_" end
    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    tmpq = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmph = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmpp = zeros(Float32,nlon,nlat,nlvl,ndt)
    
    pvec  = zeros(Float32,nlvl,ndt)
    hq    = zeros(Float32,nlvl,ndt)
    dhqdp = zeros(Float32,nlvl,ndt)

    pbs = dropdims(sum(pbs[lon1:lon2,lat1:lat2,:] .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm

    if iszero(days)
        dsh = NCDataset(datadir("wrf3","3D","$(iso)QVAPOR-daily-$timestr.nc"))
        dsq = NCDataset(datadir("wrf3","3D","QVAPOR-daily-$timestr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
    else
        dsh = NCDataset(datadir("wrf3","3D","$(iso)QVAPOR-daily-$timestr-$smthstr.nc"))
        dsq = NCDataset(datadir("wrf3","3D","QVAPOR-daily-$timestr-$smthstr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
    end

    NCDatasets.load!(dsh["$(iso)QVAPOR"].var,tmph,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dsq["QVAPOR"].var,tmpq,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dsp["P"].var,tmpp,lon1:lon2,lat1:lat2,:,:)

	close(dsh)
	close(dsq)
	close(dsp)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)
        
        iiq = @view tmpq[:,:,:,idt]
        iih = @view tmph[:,:,:,idt]
        iip = @view tmpp[:,:,:,idt]

		q            = dropdims(sum(iiq .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
		h            = dropdims(sum(iih .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
		pvec[:,idt] .= dropdims(sum(iip .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm .+ pbs

        hq[:,idt] .= h ./ q

		for ilvl = 2 : (nlvl-1)
			dhqdp[ilvl,idt] = (hq[ilvl+1,idt]-hq[ilvl-1,idt]) ./ (pvec[ilvl+1,idt]-pvec[ilvl-1,idt])
		end
		dhqdp[1,idt] = (hq[2,idt]-hq[1,idt]) / (pvec[2,idt]-pvec[1,idt])
		dhqdp[end,idt] = (hq[end,idt]-hq[end-1,idt]) / (pvec[end,idt]-pvec[end-1,idt])

    end

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-dhqdp-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-dhqdp-daily-$timestr-$smthstr.nc")
    end
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["level"] = nlvl
    ds.dim["date"]  = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncpres = defVar(ds,"P",Float32,("level","date",),attrib=Dict(
        "units" => "Pa",
        "long_name" => "Pressure"
    ))

    nchq = defVar(ds,"$(iso)hq",Float32,("level","date",),attrib=Dict(
        "units" => "‰",
        "long_name" => "Depletion of $(iso)VAPOR relative to SMOW"
    ))

    ncdhqdp = defVar(ds,"$(iso)dhqdp",Float32,("level","date",),attrib=Dict(
        "units" => "‰ Pa**-1",
        "long_name" => "Gradient of $(iso)VAPOR/QVAPOR (relative to SMOW) against pressure"
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpres[:,:] = pvec
    nchq[:,:]   = (hq .- 1) * 1000
    ncdhqdp[:,:] = dhqdp * 1000

    close(ds)

end

function wrfdhdq(
    geo   :: GeoRegion;
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)
    
    ds  = NCDataset(datadir("wrf3","grid.nc"))
    lon = ds["longitude"][:,:]
    lat = ds["latitude"][:,:]
    pbs = ds["pressure_base"][:,:,:]
    close(ds)

    ggrd = RegionGrid(geo,Point2.(lon,lat))
    lon1 = minimum(ggrd.ilon); lon2 = maximum(ggrd.ilon)
    lat1 = minimum(ggrd.ilat); lat2 = maximum(ggrd.ilat)

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtm = sum(wgts)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    tmpq = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmph = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmpo = zeros(Float32,nlon,nlat,nlvl,ndt)
    tmpp = zeros(Float32,nlon,nlat,nlvl,ndt)
    
    pvec  = zeros(Float32,nlvl-1,ndt)
    dhdq  = zeros(Float32,nlvl-1,ndt)
    dodq  = zeros(Float32,nlvl-1,ndt)

    pbs = dropdims(sum(pbs[lon1:lon2,lat1:lat2,:] .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm

    if iszero(days)
        dsh = NCDataset(datadir("wrf3","3D","HDO_QVAPOR-daily-$timestr.nc"))
        dso = NCDataset(datadir("wrf3","3D","O18_QVAPOR-daily-$timestr.nc"))
        dsq = NCDataset(datadir("wrf3","3D","QVAPOR-daily-$timestr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
    else
        dsh = NCDataset(datadir("wrf3","3D","HDO_QVAPOR-daily-$timestr-$smthstr.nc"))
        dso = NCDataset(datadir("wrf3","3D","O18_QVAPOR-daily-$timestr-$smthstr.nc"))
        dsq = NCDataset(datadir("wrf3","3D","QVAPOR-daily-$timestr-$smthstr.nc"))
        dsp = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
    end

    NCDatasets.load!(dsh["HDO_QVAPOR"].var,tmph,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dso["O18_QVAPOR"].var,tmpo,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dsq["QVAPOR"].var,tmpq,lon1:lon2,lat1:lat2,:,:)
    NCDatasets.load!(dsp["P"].var,tmpp,lon1:lon2,lat1:lat2,:,:)

	close(dsh)
	close(dso)
	close(dsq)
	close(dsp)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)
        
        iiq = @view tmpq[:,:,:,idt]
        iih = @view tmph[:,:,:,idt]
        iio = @view tmpo[:,:,:,idt]
        iip = @view tmpp[:,:,:,idt]

		q = dropdims(sum(iiq .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
		h = dropdims(sum(iih .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
	    o = dropdims(sum(iio .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
		p = dropdims(sum(iip .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm .+ pbs

        dhdq[:,idt] .= (h[1:(end-1)] .- h[2:end]) ./ (q[1:(end-1)] .- q[2:end])
        dodq[:,idt] .= (o[1:(end-1)] .- o[2:end]) ./ (q[1:(end-1)] .- q[2:end])
        pvec[:,idt] .= (p[1:(end-1)] .+ p[2:end]) / 2

    end

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-dhodq-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-dhodq-daily-$timestr-$smthstr.nc")
    end
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["level"] = nlvl - 1
    ds.dim["date"]  = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncpres = defVar(ds,"P",Float32,("level","date",),attrib=Dict(
        "units" => "Pa",
        "long_name" => "Pressure"
    ))

    ncdhdq = defVar(ds,"dHDOdH2O",Float32,("level","date",),attrib=Dict(
        "long_name" => "Gradient of HDO_QVAPOR/Gradient of QVAPOR (relative to SMOW) against pressure"
    ))

    ncdodq = defVar(ds,"dO18dH2O",Float32,("level","date",),attrib=Dict(
        "long_name" => "Gradient of O18_QVAPOR/Gradient of QVAPOR (relative to SMOW) against pressure"
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpres[:,:]   = pvec
    ncdhdq[:,:]   = dhdq
    ncdodq[:,:]   = dodq

    close(ds)

end

function wrfcdhdq(
    geo   :: GeoRegion;
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-dhodq-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-dhodq-daily-$timestr-$smthstr.nc")
    end

    ds = NCDataset(fnc)
    time = ds["time"][:]; nt = length(time)
    pvec = ds["P"][:,:]; np = size(pvec,1)
    dhdq = ds["dHDOdH2O"][:,:]
    dodq = ds["dO18dH2O"][:,:]
    close(ds)

    x = ones(np,2)
    cdhdq = zeros(2,nt)
    cdodq = zeros(2,nt)

    for it = 1 : nt

        ip = pvec[:,it]; ii = ip .> 500e2
        iip = ip[ii]
        iidhdq = dhdq[ii,it]
        iidodq = dodq[ii,it]
        ix = @views x[ii,:]
        ix[:,2] .= iip
        cdhdq[:,it] = ix \ iidhdq
        cdodq[:,it] = ix \ iidodq

    end

    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-cdhodq-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-cdhodq-daily-$timestr-$smthstr.nc")
    end
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]  = nt
    ds.dim["coeff"] = 2

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncdhdq = defVar(ds,"cdHDOdH2O",Float32,("coeff","date",),attrib=Dict(
        "long_name" => "Coefficient of Gradient of HDO_QVAPOR/Gradient of QVAPOR (relative to SMOW) against pressure"
    ))

    ncdodq = defVar(ds,"cdO18dH2O",Float32,("coeff","date",),attrib=Dict(
        "long_name" => "Coefficient of Gradient of O18_QVAPOR/Gradient of QVAPOR (relative to SMOW) against pressure"
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncdhdq[:,:]   = cdhdq
    ncdodq[:,:]   = cdodq

    close(ds)

end