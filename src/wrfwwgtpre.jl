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

    warr = zeros(Float32,nlon,nlat,nlvl+1)
    parr = zeros(Float32,nlon,nlat,nlvl)
    tarr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)

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

    for ii in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for Day $ii of $ndt"
        flush(stderr)

        NCDatasets.load!(wds["W"].var,warr,lon1:lon2,lat1:lat2,:,ii)
        NCDatasets.load!(pds["P"].var,parr,lon1:lon2,lat1:lat2,:,ii)
        NCDatasets.load!(tds["T"].var,tarr,lon1:lon2,lat1:lat2,:,ii)
        NCDatasets.load!(sds["PSFC"].var,psfc,lon1:lon2,lat1:lat2,ii)

        for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
            parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            tarr[ilon,ilat,ilvl] += 290
            tarr[ilon,ilat,ilvl]  = tarr[ilon,ilat,ilvl] * (100000 / parr[ilon,ilat,ilvl]) ^ (287/1004)
        end

        for ilat = 1 : nlat, ilon = 1 : nlon
            
            tmp_pmat[ilon,ilat,1] = psfc[ilon,ilat]

            for ilvl = 1 : nlvl
                tmp_wmat[ilon,ilat,ilvl+1]  = (warr[ilon,ilat,ilvl] + warr[ilon,ilat,ilvl+1]) / 2
                tmp_ρmat[ilon,ilat,ilvl+1]  =  parr[ilon,ilat,ilvl] / Rd / tarr[ilon,ilat,ilvl]
                tmp_wmat[ilon,ilat,ilvl+1] *= (tmp_ρmat[ilon,ilat,ilvl+1] * -9.81)
                tmp_pmat[ilon,ilat,ilvl+1]  = parr[ilon,ilat,ilvl]
            end

        end

        tmp_wvec = dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2))
        tmp_pvec = dropdims(mean(tmp_pmat,dims=(1,2)),dims=(1,2))

        calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
        if (calc > 0) & (calc < mean(psfc))
            pwgt[ii] = calc
            σwgt[ii] = calc / mean(psfc)
        else
            pwgt[ii] = NaN32
            σwgt[ii] = NaN32
        end
        wvec[:,ii] = reverse(tmp_wvec)
        pvec[:,ii] = reverse(tmp_pvec)

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

function wrfwwgtpre_monthly()

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    close(ds)

    nlon = size(lon,1)
    nlat = size(lat,2)
    nlvl = 50

    warr = zeros(Float32,nlon,nlat,nlvl+1)
    parr = zeros(Float32,nlon,nlat,nlvl)
    tarr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)

    tmp_wvec = zeros(Float32,52)
    tmp_pvec = zeros(Float32,52)
    tmp_ρvec = zeros(Float32,52)

    pwgt = zeros(Float32,nlon,nlat,12)
    σwgt = zeros(Float32,nlon,nlat,12)
    psfc = zeros(Float32,nlon,nlat,12)
    rain = zeros(Float32,nlon,nlat,12)

    pds  = NCDataset(datadir("wrf3","3D","PB-daily.nc"))
    pbse = pds["PB"][:,:,:,1]
    tt   = pds["time"][:]
    close(pds)

    wds = NCDataset(datadir("wrf3","3D","W-daily.nc"))
    pds = NCDataset(datadir("wrf3","3D","P-daily.nc"))
    tds = NCDataset(datadir("wrf3","3D","T-daily.nc"))
    sds = NCDataset(datadir("wrf3","2D","PSFC-daily.nc"))
    rds = NCDataset(datadir("wrf3","2D","RAINNC-daily.nc"))

    for it = 1 : 12

        ibeg = findfirst(month.(tt) .== it)
        iend = findlast(month.(tt) .== it)

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

        mkpath(datadir("wrf3","processed"))
        fnc = datadir("wrf3","processed","p_wwgt-wrf-monthly.nc")
        if isfile(fnc); rm(fnc,force=true) end

        ds = NCDataset(fnc,"c")
        ds.dim["longitude"] = nlon
        ds.dim["latitude"]  = nlat
        ds.dim["month"]     = 12

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

        ncpwgt[:,:,:] = pwgt
        ncσwgt[:,:,:] = σwgt
        ncrain[:,:,:] = rain
        ncpsfc[:,:,:] = psfc

        close(ds)

    end

    close(wds)
    close(pds)
    close(tds)
    close(sds)
    close(rds)

end
