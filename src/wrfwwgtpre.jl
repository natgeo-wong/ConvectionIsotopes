using Dates
using GeoRegions
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

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    close(ds)

    ggrd = RegionGrid(geo,lon,lat)
    apnt = findall(ggrd.mask .== 1)
    npnt = length(apnt)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findfirst(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findfirst(ggrd.mask .== 1)[2]

    for ipnt = 2 : npnt
        ilon = apnt[ipnt][1]; ilat = apnt[ipnt][2]
        if ilon < lon1; lon1 = ilon; end
        if ilon > lon2; lon2 = ilon; end
        if ilat < lat1; lat1 = ilat; end
        if ilat > lat2; lat2 = ilat; end
    end

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    warr = zeros(Float32,nlon,nlat,nlvl+1)
    parr = zeros(Float32,nlon,nlat,nlvl)
    tarr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)

    tmp_wmat = zeros(Float32,nlon,nlat,52); tmp_wvec = zeros(Float32,52)
    tmp_pmat = zeros(Float32,nlon,nlat,52); tmp_pvec = zeros(Float32,52)
    tmp_ρmat = zeros(Float32,nlon,nlat,52)

    pwgt = zeros(Float32,ndt)
    σwgt = zeros(Float32,ndt)

    pds  = NCDataset(datadir("wrf","3D","PB-daily.nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    if iszero(days)
        wds = NCDataset(datadir("wrf","3D","W-daily.nc"))
        pds = NCDataset(datadir("wrf","3D","P-daily.nc"))
        tds = NCDataset(datadir("wrf","3D","T-daily.nc"))
        sds = NCDataset(datadir("wrf","2D","PSFC-daily.nc"))
    else
        smthstr = "smooth_$(@sprintf("%02d",days))days"
        wds = NCDataset(datadir("wrf","3D","W-daily-$smthstr.nc"))
        pds = NCDataset(datadir("wrf","3D","P-daily-$smthstr.nc"))
        tds = NCDataset(datadir("wrf","3D","T-daily-$smthstr.nc"))
        sds = NCDataset(datadir("wrf","2D","PSFC-daily-$smthstr.nc"))
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

        tmp_wvec = reverse(dropdims(mean(tmp_wmat,dims=(1,2)),dims=(1,2)))
        tmp_pvec = reverse(dropdims(mean(tmp_pmat,dims=(1,2)),dims=(1,2)))

        pwgt[ii] = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
        
        σwgt[ii] = pwgt[ii] / mean(psfc) 

    end

    close(wds)
    close(pds)
    close(tds)
    close(sds)

    mkpath(datadir("wrf","processed"))
    if iszero(days)
        fnc = datadir("wrf","processed","$(geo.ID)-p_wwgt-daily.nc")
    else
        smthstr = "smooth_$(@sprintf("%02d",days))days"
        fnc = datadir("wrf","processed","$(geo.ID)-p_wwgt-daily-$smthstr.nc")
    end

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]      = ndt

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
    ncpwgt[:] = pwgt
    ncσwgt[:] = σwgt

    close(ds)

end

function wrfwwgtpre()

    Rd = 287.053
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
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

    pwgt = zeros(Float32,nlon,nlat)
    σwgt = zeros(Float32,nlon,nlat)

    pds  = NCDataset(datadir("wrf","3D","PB-daily.nc"))
    pbse = pds["PB"][:,:,:,1]
    close(pds)

    wds = NCDataset(datadir("wrf","3D","W-daily.nc"))
    pds = NCDataset(datadir("wrf","3D","P-daily.nc"))
    tds = NCDataset(datadir("wrf","3D","T-daily.nc"))
    sds = NCDataset(datadir("wrf","2D","PSFC-daily.nc"))

    warr = dropdims(mean(wds["W"][:,:,:,:],dims=4),dims=4)
    parr = dropdims(mean(pds["P"][:,:,:,:],dims=4),dims=4)
    tarr = dropdims(mean(tds["T"][:,:,:,:],dims=4),dims=4)
    psfc = dropdims(mean(sds["PSFC"][:,:,:],dims=3),dims=3)

    close(wds)
    close(pds)
    close(tds)
    close(sds)

    for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
        parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
        tarr[ilon,ilat,ilvl] += 290
        tarr[ilon,ilat,ilvl]  = tarr[ilon,ilat,ilvl] * (100000 / parr[ilon,ilat,ilvl]) ^ (287/1004)
    end        

    for ilat = 1 : nlat, ilon = 1 : nlon
        
        tmp_pvec[1] = psfc[ilon,ilat]

        for ilvl = 1 : nlvl
            tmp_ρvec[ilvl+1]  =  parr[ilon,ilat,ilvl] / Rd / tarr[ilon,ilat,ilvl]
            tmp_wvec[ilvl+1]  = (warr[ilon,ilat,ilvl] + warr[ilon,ilat,ilvl+1]) / 2
            tmp_wvec[ilvl+1] *= (tmp_ρvec[ilvl+1] * -9.81)
            tmp_pvec[ilvl+1]  =  parr[ilon,ilat,ilvl]
        end

        calc = trapz(tmp_pvec,tmp_wvec.*tmp_pvec) / trapz(tmp_pvec,tmp_wvec)
        pwgt[ilon,ilat] = calc
        σwgt[ilon,ilat] = calc / psfc[ilon,ilat]

    end

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","p_wwgt-wrf.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat

    ncpwgt = defVar(ds,"p_wwgt",Float32,("longitude","latitude",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name" => "Vertical Wind Weighted Column Pressure",
        "units"     => "Pa",
    ))

    ncσwgt = defVar(ds,"σ_wwgt",Float32,("longitude","latitude",),attrib=Dict(
        "long_name" => "column_mean_lagrangian_tendency_of_sigma",
        "full_name" => "Vertical Wind Weighted Column Sigma",
        "units"     => "0-1",
    ))

    ncpwgt[:,:] = pwgt
    ncσwgt[:,:] = σwgt

    close(ds)

end