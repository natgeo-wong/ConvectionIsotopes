using Dates
using GeoRegions
using RegionGrids
using NCDatasets
using Printf
using Statistics
using Trapz

include(srcdir("backend.jl"))

function wrfwwgtpre(
    geo  :: GeoRegion;
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    Rd = 287.053
    
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
    nlvl = 50
    ndt  = length(dtvec)

    tmp_wmat = zeros(Float32,nlon,nlat,52)
    tmp_pmat = zeros(Float32,nlon,nlat,52)
    tmp_ρmat = zeros(Float32,nlon,nlat,52)

    pwgt = zeros(Float32,ndt)
    σwgt = zeros(Float32,ndt)
    pvec = zeros(Float32,52,ndt)
    wvec = zeros(Float32,52,ndt)

    ds   = NCDataset(datadir("wrf3","grid.nc"))
    pbse = ds["pressure_base"][lon1:lon2,lat1:lat2,:]
    close(ds)

    if iszero(days)
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr.nc"))
    else
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr-$smthstr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr-$smthstr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr-$smthstr.nc"))
    end

    @info "$(now()) - Loading data in the GeoRegion ..."; flush(stderr)
    warr = wds["W"].var[lon1:lon2,lat1:lat2,:,:]
    parr = pds["P"].var[lon1:lon2,lat1:lat2,:,:] .+ pbse
    tarr = tds["T"].var[lon1:lon2,lat1:lat2,:,:] .+ 290
    psfc = sds["PSFC"].var[lon1:lon2,lat1:lat2,:]

    @views @. tarr *= (100000 / parr) ^ (287/1004)

    for it in 1 : ndt

        @info "$(now()) - Calculation for Day $it of $ndt ..."; flush(stderr)

        for ilat = 1 : nlat, ilon = 1 : nlon
            
            tmp_pmat[ilon,ilat,1] = psfc[ilon,ilat,it]

            for ilvl = 1 : nlvl
                tmp_wmat[ilon,ilat,ilvl+1]  = (warr[ilon,ilat,ilvl,it] + warr[ilon,ilat,ilvl+1,it]) / 2
                tmp_ρmat[ilon,ilat,ilvl+1]  =  parr[ilon,ilat,ilvl,it] / Rd / tarr[ilon,ilat,ilvl,it]
                tmp_wmat[ilon,ilat,ilvl+1] *= (tmp_ρmat[ilon,ilat,ilvl+1] * -9.81)
                tmp_pmat[ilon,ilat,ilvl+1]  = parr[ilon,ilat,ilvl,it]
            end

        end

        tmp_wvec = dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2))
        tmp_pvec = dropdims(mean(tmp_pmat,dims=(1,2)),dims=(1,2))
        tmp_psfc = mean(view(psfc,:,:,it))

        calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
        if (calc > 0) & (calc < tmp_psfc)
            pwgt[it] = calc
            σwgt[it] = calc / tmp_psfc
        else
            pwgt[it] = NaN32
            σwgt[it] = NaN32
        end
        wvec[:,it] = reverse(tmp_wvec)
        pvec[:,it] = reverse(tmp_pvec)

    end

    close(wds)
    close(pds)
    close(tds)
    close(sds)

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt-daily-$timestr-$smthstr.nc")
    end

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]   = ndt
    ds.dim["levels"] = 52

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

    ncp = defVar(ds,"P",Float32,("levels","date",),attrib=Dict(
        "long_name" => "pressure",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "Pa",
    ))


    ncw = defVar(ds,"W",Float32,("levels","date",),attrib=Dict(
        "long_name" => "pressure-velocity",
        "full_name" => "Pressure Velocity",
        "units"     => "Pa s**-1",
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpwgt[:] = pwgt
    ncσwgt[:] = σwgt
    ncp[:,:] = pvec
    ncw[:,:] = wvec

    close(ds)

end

function wrfwwgtpre(
    geov  :: Vector{GeoRegion};
    start :: Date,
    stop  :: Date,
	days  :: Int = 0
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]; nlon,nlat = size(lon)
    lat  = ds["latitude"][:,:]
    pbse = ds["pressure_base"][:,:,:]; nlvl = size(pbse,3)
    close(ds)
    
    ngeo  = length(geov)
    lon1v = zeros(Int64,ngeo); lat1v = zeros(Int64,ngeo)
    lon2v = zeros(Int64,ngeo); lat2v = zeros(Int64,ngeo)
    for igeo = 1 : ngeo
        ggrd = RegionGrid(geov[igeo],Point2.(lon,lat))
        lon1v[igeo] = minimum(ggrd.ilon); lon2v[igeo] = maximum(ggrd.ilon)
        lat1v[igeo] = minimum(ggrd.ilat); lat2v[igeo] = maximum(ggrd.ilat)
    end

    dtvec = start : Day(1) : stop; ndt = length(dtvec)

    tmp_wmat = zeros(Float32,nlon,nlat,nlvl+2)
    tmp_pmat = zeros(Float32,nlon,nlat,nlvl+2)
    tmp_ρmat = zeros(Float32,nlon,nlat,nlvl+2)

    pwgt = zeros(Float32,ndt,ngeo)
    σwgt = zeros(Float32,ndt,ngeo)
    pvec = zeros(Float32,nlvl+2,ndt,ngeo)
    wvec = zeros(Float32,nlvl+2,ndt,ngeo)

    if iszero(days)
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr.nc"))
    else
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr-$smthstr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr-$smthstr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr-$smthstr.nc"))
    end

    @info "$(now()) - Loading data in the GeoRegion ..."; flush(stderr)
    warr = wds["W"].var[:,:,:,:]
    parr = pds["P"].var[:,:,:,:]
    tarr = tds["T"].var[:,:,:,:]
    psfc = sds["PSFC"].var[:,:,:]

    close(wds)
    close(pds)
    close(tds)
    close(sds)

    @info "$(now()) - Doing some data adjustment ..."; flush(stderr)
    @views @. parr += pbse
    @views @. tarr += 290
    @views @. tarr *= (100000 / parr) ^ (287/1004)

    for igeo = 1 : ngeo

        lon1 = lon1v[igeo]; lon2 = lon2v[igeo]; nilon = lon2 - lon1 + 1; lonv = lon1:lon2
        lat1 = lat1v[igeo]; lat2 = lat2v[igeo]; nilat = lat2 - lat1 + 1; latv = lat1:lat2

        tmp_wmat = zeros(Float32,nilon,nilat,nlvl+2)
        tmp_pmat = zeros(Float32,nilon,nilat,nlvl+2)
        tmp_ρmat = zeros(Float32,nilon,nilat,nlvl+2)

        for it in 1 : ndt

            @info "$(now()) - Calculation for Day $it of $ndt for $(geov[igeo].ID) ..."; flush(stderr)

            @views @. tmp_pmat[:,:,1] = psfc[lonv,latv,it]
            for ilvl = 1 : nlvl
                @views @. tmp_wmat[:,:,ilvl+1] = (warr[lonv,latv,ilvl,it] + 
                                                  warr[lonv,latv,ilvl+1,it]) / 2
                @views @. tmp_ρmat[:,:,ilvl+1] = parr[lonv,latv,ilvl,it] / Rd / 
                                                 tarr[lonv,latv,ilvl,it]
            end
            @views @. tmp_wmat[:,:,2:(nlvl+1)] *= (tmp_ρmat[lonv,latv,:,it] * -9.81)
            @views @. tmp_pmat[:,:,2:(nlvl+1)]  =      parr[lonv,latv,:,it]

            tmp_wvec = dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2))
            tmp_pvec = dropdims(mean(tmp_pmat,dims=(1,2)),dims=(1,2))
            tmp_psfc = mean(view(psfc,lonv,latv,it))

            calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
            if (calc > 0) & (calc < mean(tmp_psfc))
                pwgt[it,igeo] = calc
                σwgt[it,igeo] = calc / mean(tmp_psfc)
            else
                pwgt[it,igeo] = NaN32
                σwgt[it,igeo] = NaN32
            end
            wvec[:,it,igeo] = reverse(tmp_wvec)
            pvec[:,it,igeo] = reverse(tmp_pvec)

        end

    end

    mkpath(datadir("wrf3","processed"))
    for igeo = 1 : ngeo
        if iszero(days)
            fnc = datadir("wrf3","processed","$(geov[igeo].ID)-p_wwgt-daily-$timestr.nc")
        else
            fnc = datadir("wrf3","processed","$(geov[igeo].ID)-p_wwgt-daily-$timestr-$smthstr.nc")
        end

        if isfile(fnc); rm(fnc,force=true) end

        ds = NCDataset(fnc,"c")
        ds.dim["date"]   = ndt
        ds.dim["levels"] = 52

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

        ncp = defVar(ds,"P",Float32,("levels","date",),attrib=Dict(
            "long_name" => "pressure",
            "full_name" => "Vertical Wind Weighted Column Sigma",
            "units"     => "Pa",
        ))


        ncw = defVar(ds,"W",Float32,("levels","date",),attrib=Dict(
            "long_name" => "pressure-velocity",
            "full_name" => "Pressure Velocity",
            "units"     => "Pa s**-1",
        ))

        nctime.var[:] = collect(0 : (ndt-1))
        ncpwgt[:] = pwgt[:,igeo]
        ncσwgt[:] = σwgt[:,igeo]
        ncp[:,:] = pvec[:,:,igeo]
        ncw[:,:] = wvec[:,:,igeo]

        close(ds)
    end

end

function wrfwwgtpre(;
    start :: Date,
    stop  :: Date,
    days  :: Int = 0
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    pbse = ds["pressure_base"][:,:,:]
    close(ds)

    dtvec = start : Day(1) : stop

    nlon = size(lon,1)
    nlat = size(lat,2)
    nlvl = 50
    ndt  = length(dtvec)

    warr = zeros(Float32,nlon,nlat,nlvl+1)
    parr = zeros(Float32,nlon,nlat,nlvl)
    tarr = zeros(Float32,nlon,nlat,nlvl)

    tmp_wvec = zeros(Float32,52)
    tmp_pvec = zeros(Float32,52)
    tmp_ρvec = zeros(Float32,52)

    pwgt = zeros(Float32,nlon,nlat,ndt)
    σwgt = zeros(Float32,nlon,nlat,ndt)
    psfc = zeros(Float32,nlon,nlat,ndt)

    if iszero(days)
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr.nc"))
    else
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr-$smthstr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr-$smthstr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr-$smthstr.nc"))
    end

    NCDatasets.load!(sds["PSFC"].var,psfc,:,:,:)

    for idt = 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Calculating p_wwgt and σ_wwgt for $idt"
		flush(stderr)

        NCDatasets.load!(wds["W"].var,warr,:,:,:,idt)
        NCDatasets.load!(pds["P"].var,parr,:,:,:,idt)
        NCDatasets.load!(tds["T"].var,tarr,:,:,:,idt)

        for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
            parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            tarr[ilon,ilat,ilvl] += 290
            tarr[ilon,ilat,ilvl]  = tarr[ilon,ilat,ilvl] * (100000 / parr[ilon,ilat,ilvl]) ^ (287/1004)
        end

        for ilat = 1 : nlat, ilon = 1 : nlon
            
            tmp_pvec[1] = psfc[ilon,ilat,idt]

            for ilvl = 1 : nlvl
                tmp_ρvec[ilvl+1]  =  parr[ilon,ilat,ilvl] / Rd / tarr[ilon,ilat,ilvl]
                tmp_wvec[ilvl+1]  = (warr[ilon,ilat,ilvl] + warr[ilon,ilat,ilvl+1]) / 2
                tmp_wvec[ilvl+1] *= (tmp_ρvec[ilvl+1] * -9.81)
                tmp_pvec[ilvl+1]  =  parr[ilon,ilat,ilvl]
            end

            calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
            if (calc > 0) & (calc < psfc[ilon,ilat,idt])
                pwgt[ilon,ilat,idt] = calc
                σwgt[ilon,ilat,idt] = calc / psfc[ilon,ilat,idt]
            else
                pwgt[ilon,ilat,idt] = NaN32
                σwgt[ilon,ilat,idt] = NaN32
            end

        end

    end

    close(wds)
    close(pds)
    close(tds)
    close(sds)

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","2D","p_wwgt-$timestr.nc")
    else
        fnc = datadir("wrf3","2D","p_wwgt-$timestr-$smthstr.nc")
    end
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["time"]      = ndt

    nctime = defVar(ds,"time",Int32,("time",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncpwgt = defVar(ds,"p_wwgt",Float32,("longitude","latitude","time",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "Pa",
    ))

    ncσwgt = defVar(ds,"σ_wwgt",Float32,("longitude","latitude","time",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "0-1",
    ))

    ncpsfc = defVar(ds,"PSFC",Float32,("longitude","latitude","time",),attrib=Dict(
        "long_name" => "surface_pressure",
        "full_name" => "Surface Pressure",
        "units"     => "Pa",
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpwgt[:,:,:] = pwgt
    ncσwgt[:,:,:] = σwgt
    ncpsfc[:,:,:] = psfc

    close(ds)

end

function wrfwwgtpre_compiled(;
    start :: Date,
    stop  :: Date
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    pbse = ds["pressure_base"][:,:,:]
    close(ds)

    nlon = size(lon,1)
    nlat = size(lat,2)
    nlvl = 50

    dtvec = start : Day(1) : stop
    ndt = length(dtvec)
 
    wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr.nc"))
    pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
    tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr.nc"))
    sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr.nc"))
    rds = NCDataset(datadir("wrf3","2D","RAINNC-daily-$timestr.nc"))

    warr = dropdims(mean(wds["W"][:,:,:,1:(end-1)],dims=4),dims=4)
    parr = dropdims(mean(pds["P"][:,:,:,1:(end-1)],dims=4),dims=4)
    tarr = dropdims(mean(tds["T"][:,:,:,1:(end-1)],dims=4),dims=4)
    psfc = dropdims(mean(sds["PSFC"][:,:,1:(end-1)],dims=3),dims=3)
    rain = dropdims(mean(rds["RAINNC"][:,:,1:(end-1)],dims=3),dims=3)

    close(wds)
    close(pds)
    close(tds)
    close(sds)
    close(rds)

    for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
        parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
        tarr[ilon,ilat,ilvl] += 290
        tarr[ilon,ilat,ilvl]  = tarr[ilon,ilat,ilvl] * (100000 / parr[ilon,ilat,ilvl]) ^ (287/1004)
    end

    tmp_wvec = zeros(Float32,52)
    tmp_pvec = zeros(Float32,52)
    tmp_ρvec = zeros(Float32,52)

    pwgt = zeros(Float32,nlon,nlat)
    σwgt = zeros(Float32,nlon,nlat)     

    for ilat = 1 : nlat, ilon = 1 : nlon
        
        tmp_pvec[1] = psfc[ilon,ilat]

        for ilvl = 1 : nlvl
            tmp_ρvec[ilvl+1]  =  parr[ilon,ilat,ilvl] / Rd / tarr[ilon,ilat,ilvl]
            tmp_wvec[ilvl+1]  = (warr[ilon,ilat,ilvl] + warr[ilon,ilat,ilvl+1]) / 2
            tmp_wvec[ilvl+1] *= (tmp_ρvec[ilvl+1] * -9.81)
            tmp_pvec[ilvl+1]  =  parr[ilon,ilat,ilvl]
        end

        calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
        if (calc > 0) & (calc < psfc[ilon,ilat])
            pwgt[ilon,ilat] = calc
            σwgt[ilon,ilat] = calc / psfc[ilon,ilat]
        else
            pwgt[ilon,ilat] = NaN32
            σwgt[ilon,ilat] = NaN32
        end

    end

    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","p_wwgt-compiledwrf-$timestr.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat

    ncpwgt = defVar(ds,"p_wwgt",Float32,("longitude","latitude"),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "Pa",
    ))

    ncσwgt = defVar(ds,"σ_wwgt",Float32,("longitude","latitude"),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "0-1",
    ))

    ncrain = defVar(ds,"RAINNC",Float32,("longitude","latitude"),attrib=Dict(
        "long_name" => "total_precipitation",
        "full_name" => "Total Precipitation",
        "units"     => "mm",
    ))

    ncpsfc = defVar(ds,"PSFC",Float32,("longitude","latitude"),attrib=Dict(
        "long_name" => "surface_pressure",
        "full_name" => "Surface Pressure",
        "units"     => "Pa",
    ))

    ncpwgt[:,:] = pwgt
    ncσwgt[:,:] = σwgt
    ncrain[:,:] = rain
    ncpsfc[:,:] = psfc

    close(ds)

end

function wrfwwgtpre_monthly(;
    start :: Date,
    stop  :: Date
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    pbse = ds["pressure_base"][:,:,:]
    close(ds)

    nlon = size(lon,1)
    nlat = size(lat,2)
    nlvl = 50

    dtvec = start : Month(1) : stop
    ndt = length(dtvec)

    warr = zeros(Float32,nlon,nlat,nlvl+1)
    parr = zeros(Float32,nlon,nlat,nlvl)
    tarr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)

    tmp_wvec = zeros(Float32,52)
    tmp_pvec = zeros(Float32,52)
    tmp_ρvec = zeros(Float32,52)

    pwgt = zeros(Float32,nlon,nlat,ndt)
    σwgt = zeros(Float32,nlon,nlat,ndt)
    psfc = zeros(Float32,nlon,nlat,ndt)
    rain = zeros(Float32,nlon,nlat,ndt)
 
    wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr.nc"))
    pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr.nc"))
    tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr.nc"))
    sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr.nc"))
    rds = NCDataset(datadir("wrf3","2D","RAINNC-daily-$timestr.nc"))
    
    tt = rds["time"][:]

    for it = 1 : ndt

        iyr = year(dtvec[it]); imo = month(dtvec[it]); ndy = daysinmonth(iyr,imo)
        ibeg = findfirst(tt .>= Date(iyr,imo,1))
        iend = findlast(tt .<= Date(iyr,imo,ndy))
        if iend == length(tt); iend -= 1 end

        warr = dropdims(mean(wds["W"][:,:,:,ibeg:iend],dims=4),dims=4)
        parr = dropdims(mean(pds["P"][:,:,:,ibeg:iend],dims=4),dims=4)
        tarr = dropdims(mean(tds["T"][:,:,:,ibeg:iend],dims=4),dims=4)
        psfc[:,:,it] = dropdims(mean(sds["PSFC"][:,:,ibeg:iend],dims=3),dims=3)
        rain[:,:,it] = dropdims(mean(rds["RAINNC"][:,:,ibeg:iend],dims=3),dims=3)

        for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
            parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            tarr[ilon,ilat,ilvl] += 290
            tarr[ilon,ilat,ilvl]  = tarr[ilon,ilat,ilvl] * (100000 / parr[ilon,ilat,ilvl]) ^ (287/1004)
        end        

        for ilat = 1 : nlat, ilon = 1 : nlon
            
            tmp_pvec[1] = psfc[ilon,ilat,it]

            for ilvl = 1 : nlvl
                tmp_ρvec[ilvl+1]  =  parr[ilon,ilat,ilvl] / Rd / tarr[ilon,ilat,ilvl]
                tmp_wvec[ilvl+1]  = (warr[ilon,ilat,ilvl] + warr[ilon,ilat,ilvl+1]) / 2
                tmp_wvec[ilvl+1] *= (tmp_ρvec[ilvl+1] * -9.81)
                tmp_pvec[ilvl+1]  =  parr[ilon,ilat,ilvl]
            end

            calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
            if (calc > 0) & (calc < psfc[ilon,ilat,it])
                pwgt[ilon,ilat,it] = calc
                σwgt[ilon,ilat,it] = calc / psfc[ilon,ilat,it]
            else
                pwgt[ilon,ilat,it] = NaN32
                σwgt[ilon,ilat,it] = NaN32
            end

        end

    end

    close(wds)
    close(pds)
    close(tds)
    close(sds)
    close(rds)

    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","p_wwgt-compiledwrf-monthly-$timestr.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["month"]     = ndt

    ncdate = defVar(ds,"date",Float32,("month",),attrib=Dict(
        "long_name" => "months since $(start) 00:00:00.0",
    ))

    ncpwgt = defVar(ds,"p_wwgt",Float32,("longitude","latitude","month"),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "Pa",
    ))

    ncσwgt = defVar(ds,"σ_wwgt",Float32,("longitude","latitude","month"),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "0-1",
    ))

    ncrain = defVar(ds,"RAINNC",Float32,("longitude","latitude","month"),attrib=Dict(
        "long_name" => "total_precipitation",
        "full_name" => "Total Precipitation",
        "units"     => "mm",
    ))

    ncpsfc = defVar(ds,"PSFC",Float32,("longitude","latitude","month"),attrib=Dict(
        "long_name" => "surface_pressure",
        "full_name" => "Surface Pressure",
        "units"     => "Pa",
    ))

    ncdate[:] = collect(0:(ndt-1)) .+ 0.5
    ncpwgt[:,:,:] = pwgt
    ncσwgt[:,:,:] = σwgt
    ncrain[:,:,:] = rain
    ncpsfc[:,:,:] = psfc

    close(ds)

end
