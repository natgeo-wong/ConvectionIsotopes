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
        dqdpds = NCDataset(datadir("wrf3","processed","$(geo.ID)-dqdp-daily-$timestr.nc"))
    else
        wds = NCDataset(datadir("wrf3","3D","W-daily-$timestr-$smthstr.nc"))
        pds = NCDataset(datadir("wrf3","3D","P-daily-$timestr-$smthstr.nc"))
        tds = NCDataset(datadir("wrf3","3D","T-daily-$timestr-$smthstr.nc"))
        sds = NCDataset(datadir("wrf3","2D","PSFC-daily-$timestr-$smthstr.nc"))
        dqdpds = NCDataset(datadir("wrf3","processed","$(geo.ID)-dqdp-daily-$timestr-$smthstr.nc"))
    end

    @info "$(now()) - Loading data in the GeoRegion ..."; flush(stderr)
    warr = wds["W"].var[lon1:lon2,lat1:lat2,:,:]
    parr = pds["P"].var[lon1:lon2,lat1:lat2,:,:] .+ pbse
    tarr = tds["T"].var[lon1:lon2,lat1:lat2,:,:] .+ 300
    psfc = sds["PSFC"].var[lon1:lon2,lat1:lat2,:]

    dqdp = dqdpds["dqdp"][:,:]

    close(wds)
    close(pds)
    close(tds)
    close(sds)
    close(dqdpds)

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

        tmp_dqdp = vcat(0,dqdp[:,it],0)
        tmp_wvec = dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2)) .* tmp_dqdp
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

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt2-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt2-daily-$timestr-$smthstr.nc")
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


    ncw = defVar(ds,"OMEGA2",Float32,("levels","date",),attrib=Dict(
        "long_name" => "pressure_velocity",
        "full_name" => "Pressure Velocity",
        "units"     => "kg kg**-1 s**-1",
    ))

    nctime.var[:] = collect(0 : (ndt-1))
    ncpwgt[:] = pwgt
    ncσwgt[:] = σwgt
    ncp[:,:] = pvec
    ncw[:,:] = wvec

    close(ds)

end

function wrfwwgtpre2(
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
    @views @. tarr += 300
    @views @. tarr *= (100000 / parr) ^ (287/1004)

    for igeo = 1 : ngeo

        lon1 = lon1v[igeo]; lon2 = lon2v[igeo]; nilon = lon2 - lon1 + 1; lonv = lon1:lon2
        lat1 = lat1v[igeo]; lat2 = lat2v[igeo]; nilat = lat2 - lat1 + 1; latv = lat1:lat2

        tmp_wmat = zeros(Float32,nilon,nilat,nlvl+2)
        tmp_pmat = zeros(Float32,nilon,nilat,nlvl+2)
        tmp_ρmat = zeros(Float32,nilon,nilat,nlvl)

        if iszero(days)
            dqdpds = NCDataset(datadir("wrf3","processed","$(geov[igeo].ID)-dqdp-daily-$timestr.nc"))
        else
            dqdpds = NCDataset(datadir("wrf3","processed","$(geov[igeo].ID)-dqdp-daily-$timestr-$smthstr.nc"))
        end

        dqdp = dqdpds["dqdp"][:,:]
        close(dqdpds)

        for it in 1 : ndt

            @info "$(now()) - Calculation for Day $it of $ndt for $(geov[igeo].ID) ..."; flush(stderr)

            @views @. tmp_pmat[:,:,1] = psfc[lonv,latv,it]
            for ilvl = 1 : nlvl
                @views @. tmp_wmat[:,:,ilvl+1] = (warr[lonv,latv,ilvl,it] + 
                                                  warr[lonv,latv,ilvl+1,it]) / 2
                @views @. tmp_ρmat[:,:,ilvl] = parr[lonv,latv,ilvl,it] / Rd / 
                                                 tarr[lonv,latv,ilvl,it]
            end
            @views @. tmp_wmat[:,:,2:(nlvl+1)] *= (tmp_ρmat * -9.81)
            @views @. tmp_pmat[:,:,2:(nlvl+1)]  = parr[lonv,latv,:,it]

            tmp_dqdp = vcat(0,dqdp[:,it],0)
            tmp_wvec = dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2)) .* tmp_dqdp
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
            fnc = datadir("wrf3","processed","$(geov[igeo].ID)-p_wwgt2-daily-$timestr.nc")
        else
            fnc = datadir("wrf3","processed","$(geov[igeo].ID)-p_wwgt2-daily-$timestr-$smthstr.nc")
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


        ncw = defVar(ds,"OMEGA2",Float32,("levels","date",),attrib=Dict(
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