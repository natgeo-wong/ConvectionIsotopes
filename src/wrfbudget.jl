using Distances
using GeoRegions
using NCDatasets
using StatsBase
using Trapz

include(srcdir("backend.jl"))

function wrfqbudget(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
    lon  = ds["longitude"][:]
    lat  = ds["latitude"][:]
    close(ds)

    dtvec = start : Day(1) : stop; ndt = length(dtvec)

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

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtm = sum(wgts)
    wgt1 = sum(wgts[1,:])
    wgt2 = sum(wgts[end,:])
    wgt3 = sum(wgts[:,1])
    wgt4 = sum(wgts[:,end])
    wgtv = weights(wgts)

    if iso != ""; iso = "$(iso)_" end

    tmp1      = zeros(Float32,nlon,nlat,8)
    tmp2      = zeros(Float32,nlon,nlat)
    tmpqflx_1 = zeros(Float32,nlat,8)
    tmpqflx_2 = zeros(Float32,nlat,8)
    tmpqflx_3 = zeros(Float32,nlon,8)
    tmpqflx_4 = zeros(Float32,nlon,8)
    tmpqflx_5 = zeros(Float32,nlat)
    tmpqflx_6 = zeros(Float32,nlat)
    tmpqflx_7 = zeros(Float32,nlon)
    tmpqflx_8 = zeros(Float32,nlon)
    
    prcp = zeros(8,ndt)
    evap = zeros(8,ndt)
    tcwv = zeros(8,ndt)
    qflx = zeros(8,ndt)

    arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
    arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
    arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
    arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        ds1 = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]).nc"))
        fncii = datadir("wrf","raw","2D","$(dtvec[idt]+Day(1))-e.nc")
        if isfile(fncii)
            @info "$(now()) - ConvectionIsotopes - Tail end"
            ds2 = NCDataset(fncii)
        else
            ds2 = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]+Day(1)).nc"))
        end

        NCDatasets.load!(ds1["$(iso)RAINNC"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)RAINNC"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 8, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        prcp[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

        NCDatasets.load!(ds1["$(iso)QFX"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)QFX"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 8, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        evap[:,idt] = (vcat(tmp3[2:end],tmp4) .+ tmp3) / 2

        NCDatasets.load!(ds1["$(iso)VAPORWP"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)VAPORWP"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 8, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        tcwv[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

        NCDatasets.load!(ds1["$(iso)IWTX"].var,tmpqflx_1,lon1,lat1:lat2,:)
        NCDatasets.load!(ds1["$(iso)IWTX"].var,tmpqflx_2,lon2,lat1:lat2,:)
        NCDatasets.load!(ds1["$(iso)IWTY"].var,tmpqflx_3,lon1:lon2,lat1,:)
        NCDatasets.load!(ds1["$(iso)IWTY"].var,tmpqflx_4,lon1:lon2,lat2,:)
        for ii = 1 : 8, ilat = 1 : nlat
            tmpqflx_1[ilat,ii] *= wgts[1,ilat]
            tmpqflx_2[ilat,ii] *= wgts[end,ilat]
        end
        for ii = 1 : 8, ilon = 1 : nlon
            tmpqflx_3[ilon,ii] *= wgts[ilon,1]
            tmpqflx_4[ilon,ii] *= wgts[ilon,end]
        end
        
        tmp3 = dropdims(sum(tmpqflx_2,dims=1),dims=1) / wgt2 * arc2 .+ 
               dropdims(sum(tmpqflx_4,dims=1),dims=1) / wgt4 * arc4 .-
               dropdims(sum(tmpqflx_1,dims=1),dims=1) / wgt1 * arc1 .-
               dropdims(sum(tmpqflx_3,dims=1),dims=1) / wgt3 * arc3

        NCDatasets.load!(ds2["$(iso)IWTX"].var,tmpqflx_5,lon1,lat1:lat2,1)
        NCDatasets.load!(ds2["$(iso)IWTX"].var,tmpqflx_6,lon2,lat1:lat2,1)
        NCDatasets.load!(ds2["$(iso)IWTY"].var,tmpqflx_7,lon1:lon2,lat1,1)
        NCDatasets.load!(ds2["$(iso)IWTY"].var,tmpqflx_8,lon1:lon2,lat2,1)
        for ilat = 1 : nlat
            tmpqflx_5[ilat] *= wgts[1,ilat]
            tmpqflx_6[ilat] *= wgts[end,ilat]
        end
        for ilon = 1 : nlon
            tmpqflx_7[ilon] *= wgts[ilon,1]
            tmpqflx_8[ilon] *= wgts[ilon,end]
        end
        
        tmp4 = sum(tmpqflx_6) / wgt2 * arc2 .+ sum(tmpqflx_8) / wgt4 * arc4 .-
               sum(tmpqflx_5) / wgt1 * arc1 .- sum(tmpqflx_7) / wgt3 * arc3

        qflx[:,idt] = (vcat(tmp3[2:end],tmp4) .+ tmp3) / 2

        close(ds1)
        close(ds2)

    end

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","$(geo.ID)-$(iso)QBUDGET-daily.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt * 8

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncprcp = defVar(ds,"$(iso)P",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2",
        "long_name" => "Accumulated 3-hour Precipitation"
    ))

    ncevap = defVar(ds,"$(iso)E",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2 s**-1",
        "long_name" => "Evaporation Rate"
    ))

    nctcwv = defVar(ds,"$(iso)ΔWVP",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2",
        "long_name" => "Change in Water Vapor Path"
    ))

    ncqflx = defVar(ds,"$(iso)∇",Float32,("date",),attrib=Dict(
        "units"     => "kg m**-2 s**-1",
        "long_name" => "Divergence"
    ))

    nctime.var[:] = (collect(0 : (ndt*8 -1)) .+ 0.5) * 3
    ncprcp[:] = prcp[:]
    ncevap[:] = evap[:]
    nctcwv[:] = tcwv[:]
    ncqflx[:] = qflx[:] * 4 / ((arc2+arc4)*(arc1+arc3))

    close(ds)

end

function wrfqdiv(
    wvar  :: AbstractString,
    geo   :: GeoRegion;
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf","raw","3D","$(start).nc"))
    lon  = ds["XLONG"][:,:,1]; lat  = ds["XLAT"][:,:,1]
    nlvl = size(ds[wvar])[3]
    attrib = Dict(ds[wvar].attrib)
    attrib["units"] = "kg m**-2 (if HDO or O18, relative to SMOW)"
    close(ds)

    dtvec = start : Day(1) : stop; ndt = length(dtvec)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1] - 1; lon2 = findlast(ggrd.mask .== 1)[1] + 1
    lat1 = findfirst(ggrd.mask .== 1)[2] - 1; lat2 = findlast(ggrd.mask .== 1)[2] + 1

    nlon = lon2 - lon1
    nlat = lat2 - lat1

    uarr = zeros(Float32,nlon  ,nlat+1,nlvl)
    varr = zeros(Float32,nlon+1,nlat  ,nlvl)
    qarr = zeros(Float32,nlon+1,nlat+1,nlvl)
    parr = zeros(Float32,nlon+1,nlat+1,nlvl)
    psfc = zeros(Float32,nlon+1,nlat+1)
    qflx = zeros(8,ndt)

    pds = NCDataset(datadir("wrf","3D","PB-daily.nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        ds2 = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]).nc"))
        ds3 = NCDataset(datadir("wrf","raw","3D","$(dtvec[idt]).nc"))

        for ii = 1 : 8

            NCDatasets.load!(ds3["U"].var,uarr,(lon1+1):lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds3["V"].var,varr,lon1:lon2,(lat1+1):lat2,:,ii)
            NCDatasets.load!(ds3[wvar].var,qarr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds3["P"].var,parr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds2["PSFC"].var,psfc,lon1:lon2,lat1:lat2,ii)

            for ilvl = 1 : nlvl, ilat = 1 : (nlat+1), ilon = 1 : (nlon+1)
                parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            end

            for ilat = 2 : nlat
                qavg = dropdims(mean(qarr[1:2,ilat,:],dims=1),dims=1)
                qlat = uarr[1,ilat,:] .* qavg
                plat = parr[1,ilat,:]
                qflx[ii,idt] -= trapz(
                    reverse(vcat(psfc[1,ilat],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                ) / (nlat-1)
                qavg = dropdims(mean(qarr[nlon.+(0:1),ilat,:],dims=1),dims=1)
                qlat = uarr[end,ilat,:] .* qavg
                plat = parr[end,ilat,:]
                qflx[ii,idt] += trapz(
                    reverse(vcat(psfc[end,ilat],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                ) / (nlat-1)
            end

            for ilon = 2 : nlon
                qavg = dropdims(mean(qarr[ilon,1:2,:],dims=1),dims=1)
                qlat = varr[ilon,1,:] .* qavg
                plat = parr[ilon,1,:]
                qflx[ii,idt] -= trapz(
                    reverse(vcat(psfc[ilon,1],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                ) / (nlon-1)
                qavg = dropdims(mean(qarr[ilon,nlat.+(0:1),:],dims=1),dims=1)
                qlat = varr[ilon,end,:] .* qavg
                plat = parr[ilon,end,:]
                qflx[ii,idt] += trapz(
                    reverse(vcat(psfc[ilon,end],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                ) / (nlon-1)
            end

        end

        close(ds2)
        close(ds3)

    end

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","$(geo.ID)-∇_$(wvar)-daily.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt * 8

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncvar = defVar(ds,"∇_$(wvar)",Float32,("date",),attrib=attrib)

    lonA = (lon[lon1,lat1]+lon[lon1,lat1+1]+lon[lon1+1,lat1]+lon[lon1+1,lat1+1]) / 4
    lonB = (lon[lon1,lat2]+lon[lon1,lat2-1]+lon[lon1+1,lat2]+lon[lon1+1,lat2-1]) / 4
    lonC = (lon[lon2,lat2]+lon[lon2,lat2-1]+lon[lon2-1,lat2]+lon[lon2-1,lat2-1]) / 4
    lonD = (lon[lon2,lat1]+lon[lon2,lat1+1]+lon[lon2-1,lat1]+lon[lon2-1,lat1+1]) / 4

    latA = (lat[lon1,lat1]+lat[lon1,lat1+1]+lat[lon1+1,lat1]+lat[lon1+1,lat1+1]) / 4
    latB = (lat[lon1,lat2]+lat[lon1,lat2-1]+lat[lon1+1,lat2]+lat[lon1+1,lat2-1]) / 4
    latC = (lat[lon2,lat2]+lat[lon2,lat2-1]+lat[lon2-1,lat2]+lat[lon2-1,lat2-1]) / 4
    latD = (lat[lon2,lat1]+lat[lon2,lat1+1]+lat[lon2-1,lat1]+lat[lon2-1,lat1+1]) / 4

    nctime.var[:] = (collect(0 : (ndt*8 -1))) * 3
    ncvar[:] = qflx[:] / 9.81 * 3600 * 3 * (
        haversine((lonA,latA),(lonB,latB)) + haversine((lonB,latB),(lonC,latC)) +
        haversine((lonC,latC),(lonD,latD)) + haversine((lonD,latD),(lonA,latA))
    ) / (
        (haversine((lonA,latA),(lonB,latB)) + haversine((lonC,latC),(lonD,latD))) *
        (haversine((lonB,latB),(lonC,latC)) + haversine((lonD,latD),(lonA,latA)))
    ) * 4

    close(ds)

end

function wrfqdivvsiwt(;
    iso   :: AbstractString,
    geo   :: GeoRegion,
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
    lon  = ds["longitude"][:,:,1]
    lat  = ds["latitude"][:,:,1]
    close(ds)

    dtvec = start : Day(1) : stop; ndt = length(dtvec)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50

    if iso != ""; iso = "$(iso)_" end

    uarr = zeros(Float32,nlon+1,nlat,nlvl)
    varr = zeros(Float32,nlon,nlat+1,nlvl)
    qarr = zeros(Float32,nlon,nlat,nlvl)
    parr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)
    qflx_u = zeros(nlon,nlat,8,ndt)
    qflx_v = zeros(nlon,nlat,8,ndt)
    qflxwu = zeros(Float32,nlon,nlat,8,ndt)
    qflxwv = zeros(Float32,nlon,nlat,8,ndt)

    pds = NCDataset(datadir("wrf","3D","PB-daily.nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        iiqflxu = @view qflxwu[:,:,:,idt]
        iiqflxv = @view qflxwv[:,:,:,idt]

        ds2 = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]).nc"))
        ds3 = NCDataset(datadir("wrf","raw","3D","$(dtvec[idt]).nc"))

        NCDatasets.load!(ds2["$(iso)IWTX"].var,iiqflxu,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)IWTY"].var,iiqflxv,lon1:lon2,lat1:lat2,:)

        for ii = 1 : 8

            NCDatasets.load!(ds3["U"].var,uarr,lon1:(lon2+1),lat1:lat2,:,ii)
            NCDatasets.load!(ds3["V"].var,varr,lon1:lon2,lat1:(lat2+1),:,ii)
            NCDatasets.load!(ds3["$(iso)QVAPOR"].var,qarr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds3["P"].var,parr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds2["PSFC"].var,psfc,lon1:lon2,lat1:lat2,ii)

            for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
                parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            end

            for ilat = 1 : nlat, ilon = 1 : nlon
                uavg = dropdims(mean(uarr[ilon.+(0:1),ilat,:],dims=1),dims=1)
                qlat = qarr[ilon,ilat,:] .* uavg
                plat = parr[ilon,ilat,:]
                qflx_u[ilon,ilat,ii,idt] = trapz(
                    reverse(vcat(psfc[ilon,ilat],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                )
                vavg = dropdims(mean(varr[ilon,ilat.+(0:1),:],dims=1),dims=1)
                qlat = qarr[ilon,ilat,:] .* vavg
                plat = parr[ilon,ilat,:]
                qflx_v[ilon,ilat,ii,idt] = trapz(
                    reverse(vcat(psfc[ilon,ilat],plat,0)),
                    reverse(vcat(qlat[1],qlat,0))
                )
            end

        end

        close(ds2)
        close(ds3)

    end

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","$(geo.ID)-$(iso)IWT_wrfvscalc-daily.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["date"]      = ndt * 8

    nclon = defVar(ds,"longitude",Float32,("longitude","latitude"),attrib=Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("longitude","latitude"),attrib=Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncvar_iwtxwrf = defVar(
        ds,"$(iso)IWTX_wrf",Float32,("longitude","latitude","date",),
        attrib=Dict(
            "units"     => "kg m**-1 s**-1",
            "long_name" => "WRF IVT in the X-Direction (relative to SMOW)"
    ))

    ncvar_iwtywrf = defVar(
        ds,"$(iso)IWTY_wrf",Float32,("longitude","latitude","date",),
        attrib=Dict(
            "units"     => "kg m**-1 s**-1",
            "long_name" => "WRF IVT in the Y-Direction (relative to SMOW)"
    ))

    ncvar_iwtxclc = defVar(
        ds,"$(iso)IWTX_clc",Float32,("longitude","latitude","date",),
        attrib=Dict(
            "units"     => "kg m**-1 s**-1",
            "long_name" => "Calculated IVT in the X-Direction (relative to SMOW)"
    ))

    ncvar_iwtyclc = defVar(
        ds,"$(iso)IWTY_clc",Float32,("longitude","latitude","date",),
        attrib=Dict(
            "units"     => "kg m**-1 s**-1",
            "long_name" => "Calculated IVT in the Y-Direction (relative to SMOW)"
    ))

    nclon[:] = lon[lon1:lon2,lat1:lat2]
    nclat[:] = lat[lon1:lon2,lat1:lat2]
    nctime.var[:] = (collect(0 : (ndt*8 -1))) * 3
    ncvar_iwtxwrf[:] = reshape(qflxwu,nlon,nlat,:)
    ncvar_iwtywrf[:] = reshape(qflxwv,nlon,nlat,:)
    ncvar_iwtxclc[:] = reshape(qflx_u,nlon,nlat,:)
    ncvar_iwtyclc[:] = reshape(qflx_v,nlon,nlat,:)

    close(ds)

end

function wrfqdivdecompose(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
    lon  = ds["longitude"][:]
    lat  = ds["latitude"][:]
    close(ds)

    dtvec = start : Day(1) : stop; ndt = length(dtvec)

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

    nlon = lon2 - lon1
    nlat = lat2 - lat1
    nlvl = 50

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtv  = weights(wgts)
    lon1w = weights(wgts[1,:]); lon2w = weights(wgts[end,:])
    lat1w = weights(wgts[:,1]); lat2w = weights(wgts[:,end])

    if iso != ""; iso = "$(iso)_" end

    utmp1 = zeros(Float32,nlon+2,nlat+1,nlvl); utmp2 = zeros(Float32,nlon+2,nlat+1,nlvl)
    vtmp1 = zeros(Float32,nlon+1,nlat+2,nlvl); vtmp2 = zeros(Float32,nlon+1,nlat+2,nlvl)

    u1 = zeros(Float32,nlon+1,nlat+1,nlvl); u2 = zeros(Float32,nlon+1,nlat+1,nlvl)
    v1 = zeros(Float32,nlon+1,nlat+1,nlvl); v2 = zeros(Float32,nlon+1,nlat+1,nlvl)
    q1 = zeros(Float32,nlon+1,nlat+1,nlvl); q2 = zeros(Float32,nlon+1,nlat+1,nlvl)
    p1 = zeros(Float32,nlon+1,nlat+1,nlvl); p2 = zeros(Float32,nlon+1,nlat+1,nlvl)

    us1 = zeros(Float32,nlon+1,nlat+1); us2 = zeros(Float32,nlon+1,nlat+1)
    vs1 = zeros(Float32,nlon+1,nlat+1); vs2 = zeros(Float32,nlon+1,nlat+1)
    ps1 = zeros(Float32,nlon+1,nlat+1); ps2 = zeros(Float32,nlon+1,nlat+1)

    μq = zeros(nlvl+2); Δqu = zeros(nlvl+2); Δqv = zeros(nlvl+2)
    μu = zeros(nlvl+2); Δu  = zeros(nlvl+2)
    μv = zeros(nlvl+2); Δv  = zeros(nlvl+2)
    
    qadv = zeros(8,ndt)
    qdiv = zeros(8,ndt)

    arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
    arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
    arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
    arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        ds1_3D = NCDataset(datadir("wrf","raw","3D","$(dtvec[idt]).nc"))
        ds1_2D = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]).nc"))
        fncii = datadir("wrf","raw","3D","$(dtvec[idt]+Day(1))-e.nc")
        if isfile(fncii)
            @info "$(now()) - ConvectionIsotopes - Tail end"
            ds2_3D = NCDataset(fncii)
            ds2_2D = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]+Day(1))-e.nc"))
        else
            ds2_3D = NCDataset(datadir("wrf","raw","3D","$(dtvec[idt]+Day(1)).nc"))
            ds2_2D = NCDataset(datadir("wrf","raw","2D","$(dtvec[idt]+Day(1)).nc"))
        end

        for it = 1 : 8

            NCDatasets.load!(ds1_3D["$(iso)QVAPOR"].var,q1,lon1:lon2,lat1:lat2,:,it)
            NCDatasets.load!(ds1_3D["P"].var,p1,lon1:lon2,lat1:lat2,:,it)
            NCDatasets.load!(ds1_3D["U"].var,utmp1,lon1:lon2+1,lat1:lat2,:,it)
            NCDatasets.load!(ds1_3D["V"].var,vtmp1,lon1:lon2,lat1:lat2+1,:,it)

            NCDatasets.load!(ds1_2D["PSFC"].var,ps1,lon1:lon2,lat1:lat2,it)
            NCDatasets.load!(ds1_2D["U10"].var,us1,lon1:lon2,lat1:lat2,it)
            NCDatasets.load!(ds1_2D["V10"].var,vs1,lon1:lon2,lat1:lat2,it)

            if it < 8
                NCDatasets.load!(ds1_3D["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1_3D["P"].var,p2,lon1:lon2,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1_3D["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1_3D["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,it+1)

                NCDatasets.load!(ds1_2D["PSFC"].var,ps2,lon1:lon2,lat1:lat2,it+1)
                NCDatasets.load!(ds1_2D["U10"].var,us2,lon1:lon2,lat1:lat2,it+1)
                NCDatasets.load!(ds1_2D["V10"].var,vs2,lon1:lon2,lat1:lat2,it+1)
            else
                NCDatasets.load!(ds2_3D["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,1)
                NCDatasets.load!(ds2_3D["P"].var,p2,lon1:lon2,lat1:lat2,:,1)
                NCDatasets.load!(ds2_3D["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,1)
                NCDatasets.load!(ds2_3D["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,1)
                
                NCDatasets.load!(ds2_2D["PSFC"].var,ps2,lon1:lon2,lat1:lat2,1)
                NCDatasets.load!(ds2_2D["U10"].var,us2,lon1:lon2,lat1:lat2,1)
                NCDatasets.load!(ds2_2D["V10"].var,vs2,lon1:lon2,lat1:lat2,1)
            end

            for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
                u1[ilon,ilat,ilvl] = (utmp1[ilon,ilat,ilvl] + utmp1[ilon+1,ilat,ilvl]) / 2
                u2[ilon,ilat,ilvl] = (utmp2[ilon,ilat,ilvl] + utmp2[ilon+1,ilat,ilvl]) / 2
                v1[ilon,ilat,ilvl] = (vtmp1[ilon,ilat,ilvl] + vtmp1[ilon,ilat+1,ilvl]) / 2
                v2[ilon,ilat,ilvl] = (vtmp2[ilon,ilat,ilvl] + vtmp2[ilon,ilat+1,ilvl]) / 2
            end

            for ilvl = 1 : nlvl

                μq[ilvl+1] = (mean(q1[:,:,ilvl],wgtv) + mean(q2[:,:,ilvl],wgtv)) / 2
                μu[ilvl+1] = (mean(u1[:,:,ilvl],wgtv) + mean(u2[:,:,ilvl],wgtv)) / 2
                μv[ilvl+1] = (mean(v1[:,:,ilvl],wgtv) + mean(v2[:,:,ilvl],wgtv)) / 2
                μp[ilvl+1] = (mean(p1[:,:,ilvl],wgtv) + mean(p2[:,:,ilvl],wgtv)) / 2
                
                Δqu[ilvl+1] = ((mean(q1[end,:,ilvl],lon2w) + 
                                mean(q2[end,:,ilvl],lon2w)) * arc2 -
                               (mean(q2[1,:,ilvl],lon1w) + 
                                mean(q2[1,:,ilvl],lon1w))   * arc1) / 2
                Δqv[ilvl+1] = ((mean(q1[:,end,ilvl],lat2w) + 
                                mean(q2[:,end,ilvl],lat2w)) * arc4 -
                               (mean(q2[:,1,ilvl],lat1w) + 
                                mean(q2[:,1,ilvl],lat1w))   * arc3) / 2
                Δu[ilvl+1] = ((mean(u1[end,:,ilvl],lon2w) + 
                               mean(u2[end,:,ilvl],lon2w)) * arc2 -
                              (mean(u2[1,:,ilvl],lon1w) + 
                               mean(u2[1,:,ilvl],lon1w))   * arc1) / 2
                Δv[ilvl+1] = ((mean(v1[:,end,ilvl],lat2w) + 
                               mean(v2[:,end,ilvl],lat2w)) * arc4 -
                              (mean(v2[:,1,ilvl],lat1w) + 
                               mean(v2[:,1,ilvl],lat1w))   * arc3) / 2
                
            end

            μp[1] = (mean(ps1,wgtv) + mean(ps2,wgtv)) / 2
            μu[1] = (mean(us1,wgtv) + mean(us2,wgtv)) / 2
            μv[1] = (mean(vs1,wgtv) + mean(vs2,wgtv)) / 2
            μq[1] = μq[2]

            Δqu[1] = Δqu[2]
            Δqv[1] = Δqv[2]

            qadv[it,idt] = trapz(reverse(μp),reverse(Δqu.*μu)) .+ 
                           trapz(reverse(μp),reverse(Δqv.*μv))
            qdiv[it,idt] = trapz(reverse(μp),reverse(Δu .*μq)) .+ 
                           trapz(reverse(μp),reverse(Δv .*μq))
        
        end

        close(ds1)
        close(ds2)

    end

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","$(geo.ID)-$(iso)∇decompose-daily.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt * 8

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncqdiv = defVar(ds,"$(iso)P",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2 s**-1",
        "long_name" => "Divergence component of ∇"
    ))

    ncqadv = defVar(ds,"$(iso)E",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2 s**-1",
        "long_name" => "Advection component of ∇"
    ))

    nctime.var[:] = (collect(0 : (ndt*8 -1)) .+ 0.5) * 3
    ncqdiv[:] = qdiv[:] * 4 / ((arc2+arc4)*(arc1+arc3))
    ncqadv[:] = qadv[:] * 4 / ((arc2+arc4)*(arc1+arc3))

    close(ds)

end