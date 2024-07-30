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
    
    ds   = NCDataset(datadir("wrf2","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    close(ds)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    ndt  = length(dtvec)

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

    tmp1      = zeros(Float32,nlon,nlat,24)
    tmp2      = zeros(Float32,nlon,nlat)
    tmpqflx_1 = zeros(Float32,nlat,24)
    tmpqflx_2 = zeros(Float32,nlat,24)
    tmpqflx_3 = zeros(Float32,nlon,24)
    tmpqflx_4 = zeros(Float32,nlon,24)
    tmpqflx_5 = zeros(Float32,nlat)
    tmpqflx_6 = zeros(Float32,nlat)
    tmpqflx_7 = zeros(Float32,nlon)
    tmpqflx_8 = zeros(Float32,nlon)
    
    prcp = zeros(24,ndt)
    evap = zeros(24,ndt)
    tcwv = zeros(24,ndt)
    qflx = zeros(24,ndt)

    arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
    arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
    arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
    arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        ds1 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]).nc"))
        fncii = datadir("wrf2","raw","$(dtvec[idt]+Day(1))-e.nc")
        if isfile(fncii)
            @info "$(now()) - ConvectionIsotopes - Tail end"
            ds2 = NCDataset(fncii)
        else
            ds2 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]+Day(1)).nc"))
        end

        NCDatasets.load!(ds1["$(iso)RAINNC"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)RAINNC"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 24, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        prcp[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

        NCDatasets.load!(ds1["$(iso)QFX"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)QFX"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 24, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        evap[:,idt] = (vcat(tmp3[2:end],tmp4) .+ tmp3) / 2

        NCDatasets.load!(ds1["$(iso)VAPORWP"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)VAPORWP"].var,tmp2,lon1:lon2,lat1:lat2,1)
        for ii = 1 : 24, ilat = 1 : nlat, ilon = 1 : nlon
            tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
        end
        tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
        tmp4 = mean(tmp2,wgtv)
        tcwv[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

        NCDatasets.load!(ds1["$(iso)IWTX"].var,tmpqflx_1,lon1,lat1:lat2,:)
        NCDatasets.load!(ds1["$(iso)IWTX"].var,tmpqflx_2,lon2,lat1:lat2,:)
        NCDatasets.load!(ds1["$(iso)IWTY"].var,tmpqflx_3,lon1:lon2,lat1,:)
        NCDatasets.load!(ds1["$(iso)IWTY"].var,tmpqflx_4,lon1:lon2,lat2,:)
        for ii = 1 : 24, ilat = 1 : nlat
            tmpqflx_1[ilat,ii] *= wgts[1,ilat]
            tmpqflx_2[ilat,ii] *= wgts[end,ilat]
        end
        for ii = 1 : 24, ilon = 1 : nlon
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

    mkpath(datadir("wrf2","processed"))
    fnc = datadir("wrf2","processed","$(geo.ID)-$(iso)QBUDGET.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt * 24

    nctime = defVar(ds,"time",Float32,("date",),attrib=Dict(
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

    nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
    ncprcp[:] = prcp[:]
    ncevap[:] = evap[:]
    nctcwv[:] = tcwv[:]
    ncqflx[:] = qflx[:] * 4 / ((arc2+arc4)*(arc1+arc3))

    close(ds)

end

function wrfqdiv(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date
)
        
    ds   = NCDataset(datadir("wrf2","grid.nc"))
    lon  = ds["longitude"][:,:]
    lat  = ds["latitude"][:,:]
    close(ds)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    if iso != ""; iso = "$(iso)_" end

    utmp1 = zeros(Float32,nlon+1,nlat  ,nlvl); utmp2 = zeros(Float32,nlon+1,nlat  ,nlvl)
    vtmp1 = zeros(Float32,nlon  ,nlat+1,nlvl); vtmp2 = zeros(Float32,nlon  ,nlat+1,nlvl)

    u1 = zeros(Float32,nlon,nlat,nlvl); u2 = zeros(Float32,nlon,nlat,nlvl)
    v1 = zeros(Float32,nlon,nlat,nlvl); v2 = zeros(Float32,nlon,nlat,nlvl)
    q1 = zeros(Float32,nlon,nlat,nlvl); q2 = zeros(Float32,nlon,nlat,nlvl)
    p1 = zeros(Float32,nlon,nlat,nlvl); p2 = zeros(Float32,nlon,nlat,nlvl)

    us1 = zeros(Float32,nlon,nlat); us2 = zeros(Float32,nlon,nlat)
    vs1 = zeros(Float32,nlon,nlat); vs2 = zeros(Float32,nlon,nlat)
    ps1 = zeros(Float32,nlon,nlat); ps2 = zeros(Float32,nlon,nlat)

    u = zeros(nlon,nlat,nlvl+2)
    v = zeros(nlon,nlat,nlvl+2)
    q = zeros(nlon,nlat,nlvl+2)
    p = zeros(nlon,nlat,nlvl+2)

    ∇ = zeros(24,ndt)

    arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
    arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
    arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
    arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

    pds = NCDataset(datadir("wrf2","raw","$(dtvec[1]).nc"))
    pb = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)
    
    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        ds1 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]).nc"))
        fncii = datadir("wrf2","raw","$(dtvec[idt]+Day(1))-e.nc")
        if isfile(fncii)
            @info "$(now()) - ConvectionIsotopes - Tail end"
            ds2 = NCDataset(fncii)
        else
            ds2 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]+Day(1)).nc"))
        end

        for it = 1 : 24

            NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q1,lon1:lon2,lat1:lat2,:,it)
            NCDatasets.load!(ds1["P"].var,p1,lon1:lon2,lat1:lat2,:,it)
            NCDatasets.load!(ds1["U"].var,utmp1,lon1:lon2+1,lat1:lat2,:,it)
            NCDatasets.load!(ds1["V"].var,vtmp1,lon1:lon2,lat1:lat2+1,:,it)

            NCDatasets.load!(ds1["PSFC"].var,ps1,lon1:lon2,lat1:lat2,it)
            NCDatasets.load!(ds1["U10"].var,us1,lon1:lon2,lat1:lat2,it)
            NCDatasets.load!(ds1["V10"].var,vs1,lon1:lon2,lat1:lat2,it)

            if it < 24
                NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1["P"].var,p2,lon1:lon2,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,it+1)
                NCDatasets.load!(ds1["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,it+1)

                NCDatasets.load!(ds1["PSFC"].var,ps2,lon1:lon2,lat1:lat2,it+1)
                NCDatasets.load!(ds1["U10"].var,us2,lon1:lon2,lat1:lat2,it+1)
                NCDatasets.load!(ds1["V10"].var,vs2,lon1:lon2,lat1:lat2,it+1)
            else
                NCDatasets.load!(ds2["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,1)
                NCDatasets.load!(ds2["P"].var,p2,lon1:lon2,lat1:lat2,:,1)
                NCDatasets.load!(ds2["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,1)
                NCDatasets.load!(ds2["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,1)
                
                NCDatasets.load!(ds2["PSFC"].var,ps2,lon1:lon2,lat1:lat2,1)
                NCDatasets.load!(ds2["U10"].var,us2,lon1:lon2,lat1:lat2,1)
                NCDatasets.load!(ds2["V10"].var,vs2,lon1:lon2,lat1:lat2,1)
            end

            for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
                u1[ilon,ilat,ilvl] = (utmp1[ilon,ilat,ilvl] + utmp1[ilon+1,ilat,ilvl]) / 2
                u2[ilon,ilat,ilvl] = (utmp2[ilon,ilat,ilvl] + utmp2[ilon+1,ilat,ilvl]) / 2
                v1[ilon,ilat,ilvl] = (vtmp1[ilon,ilat,ilvl] + vtmp1[ilon,ilat+1,ilvl]) / 2
                v2[ilon,ilat,ilvl] = (vtmp2[ilon,ilat,ilvl] + vtmp2[ilon,ilat+1,ilvl]) / 2
            end

            for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
                u[ilon,ilat,ilvl+1] = (u1[ilon,ilat,ilvl] + u2[ilon,ilat,ilvl]) / 2
                v[ilon,ilat,ilvl+1] = (v1[ilon,ilat,ilvl] + v2[ilon,ilat,ilvl]) / 2
                q[ilon,ilat,ilvl+1] = (q1[ilon,ilat,ilvl] + q2[ilon,ilat,ilvl]) / 2
                p[ilon,ilat,ilvl+1] = (p1[ilon,ilat,ilvl] + p2[ilon,ilat,ilvl]) / 2 + 
                                       pb[ilon,ilat,ilvl]
            end

            for ilat = 1 : nlat, ilon = 1 : nlon
                u[ilon,ilat,1] = (us1[ilon,ilat,1] + us2[ilon,ilat,1]) / 2
                v[ilon,ilat,1] = (vs1[ilon,ilat,1] + vs2[ilon,ilat,1]) / 2
                p[ilon,ilat,1] = (ps1[ilon,ilat,1] + ps2[ilon,ilat,1]) / 2
                q[ilon,ilat,1] = q[ilon,ilat,2]
            end

            for ilat = 2 : (nlat-1)
                ∇[it,idt] -= trapz(
                    reverse(p[1,ilat,:]),
                    reverse(q[1,ilat,:] .* u[1,ilat,:])
                ) / (nlat-1) * arc1
                ∇[it,idt] += trapz(
                    reverse(p[end,ilat,:]),
                    reverse(q[end,ilat,:] .* u[end,ilat,:])
                ) / (nlat-1) * arc2
            end

            for ilon = 2 : (nlon-1)
                ∇[it,idt] -= trapz(
                    reverse(p[ilon,1,:]),
                    reverse(q[ilon,1,:] .* v[ilon,1,:])
                ) / (nlon-1) * arc3
                ∇[it,idt] += trapz(
                    reverse(p[ilon,end,:]),
                    reverse(q[ilon,end,:] .* v[ilon,end,:])
                ) / (nlon-1) * arc4
            end

            for ilat in [1, nlat]
                ∇[it,idt] -= trapz(
                    reverse(p[1,ilat,:]),
                    reverse(q[1,ilat,:] .* u[1,ilat,:])
                ) / (nlat-1) * arc1 / 2
                ∇[it,idt] += trapz(
                    reverse(p[end,ilat,:]),
                    reverse(q[end,ilat,:] .* u[end,ilat,:])
                ) / (nlat-1) * arc2 / 2
            end

            for ilon in [1, nlon]
                ∇[it,idt] -= trapz(
                    reverse(p[ilon,1,:]),
                    reverse(q[ilon,1,:] .* v[ilon,1,:])
                ) / (nlon-1) * arc3 / 2
                ∇[it,idt] += trapz(
                    reverse(p[ilon,end,:]),
                    reverse(q[ilon,end,:] .* v[ilon,end,:])
                ) / (nlon-1) * arc4 / 2
            end
        
        end

        close(ds1_2D); close(ds1_3D)
        close(ds2_2D); close(ds2_3D)

    end

    mkpath(datadir("wrf2","processed"))
    fnc = datadir("wrf2","processed","$(geo.ID)-$(iso)∇.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt * 24

    nctime = defVar(ds,"time",Float32,("date",),attrib=Dict(
        "units"     => "hours since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncqdiv = defVar(ds,"$(iso)∇",Float32,("date",),attrib=Dict(
        "units"     => "kg m**-2 s**-1",
        "long_name" => "Divergence"
    ))

    nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
    ncqdiv[:] = ∇[:] * 4 / ((arc2+arc4)*(arc1+arc3)) / 9.81

    close(ds)

end

function wrfqdivdecompose(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date,
    overwrite :: Bool = true
)

    if iso != ""; iso = "$(iso)_" end
    fnc = datadir("wrf2","processed","$(geo.ID)-$(iso)∇decompose.nc")
    if !isfile(fnc) || overwrite

        rm(fnc,force=true)
    
        ds   = NCDataset(datadir("wrf2","grid.nc"))
        lon  = ds["longitude"][:,:]
        lat  = ds["latitude"][:,:]
        close(ds)

        ggrd = RegionGrid(geo,lon,lat)
        lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
        lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

        dtvec = start : Day(1) : stop

        nlon = lon2 - lon1 + 1
        nlat = lat2 - lat1 + 1
        nlvl = 50
        ndt  = length(dtvec)

        wgts = ones(nlon,nlat)
        wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
        wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
        wgtv  = weights(wgts)
        lon1w = weights(wgts[1,:]); lon2w = weights(wgts[end,:])
        lat1w = weights(wgts[:,1]); lat2w = weights(wgts[:,end])

        utmp1 = zeros(Float32,nlon+1,nlat  ,nlvl); utmp2 = zeros(Float32,nlon+1,nlat  ,nlvl)
        vtmp1 = zeros(Float32,nlon  ,nlat+1,nlvl); vtmp2 = zeros(Float32,nlon  ,nlat+1,nlvl)

        u1 = zeros(Float32,nlon,nlat,nlvl); u2 = zeros(Float32,nlon,nlat,nlvl)
        v1 = zeros(Float32,nlon,nlat,nlvl); v2 = zeros(Float32,nlon,nlat,nlvl)
        q1 = zeros(Float32,nlon,nlat,nlvl); q2 = zeros(Float32,nlon,nlat,nlvl)
        p1 = zeros(Float32,nlon,nlat,nlvl); p2 = zeros(Float32,nlon,nlat,nlvl)

        us1 = zeros(Float32,nlon,nlat); us2 = zeros(Float32,nlon,nlat)
        vs1 = zeros(Float32,nlon,nlat); vs2 = zeros(Float32,nlon,nlat)
        ps1 = zeros(Float32,nlon,nlat); ps2 = zeros(Float32,nlon,nlat)

        μq = zeros(nlvl+2); Δqu = zeros(nlvl+2); Δqv = zeros(nlvl+2)
        μu = zeros(nlvl+2); Δu  = zeros(nlvl+2)
        μv = zeros(nlvl+2); Δv  = zeros(nlvl+2)
        μp = zeros(nlvl+2); Δv  = zeros(nlvl+2)
        
        qadv = zeros(24,ndt)
        qdiv = zeros(24,ndt)

        arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
        arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
        arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
        arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

        pds = NCDataset(datadir("wrf2","raw","$(dtvec[1]).nc"))
        pb = pds["PB"][lon1:lon2,lat1:lat2,:,1]
        close(pds)

        for idt in 1 : ndt

            @info "$(now()) - ConvectionIsotopes - Extracting $(iso)QVAPOR data during $(dtvec[idt]) for $(geo.name)"
            flush(stderr)

            ds1 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]).nc"))
            fncii = datadir("wrf2","raw","$(dtvec[idt]+Day(1))-e.nc")
            if isfile(fncii)
                @info "$(now()) - ConvectionIsotopes - Tail end"
                ds2 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]+Day(1))-e.nc"))
            else
                ds2 = NCDataset(datadir("wrf2","raw","$(dtvec[idt]+Day(1)).nc"))
            end

            for it = 1 : 24

                NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q1,lon1:lon2,lat1:lat2,:,it)
                NCDatasets.load!(ds1["P"].var,p1,lon1:lon2,lat1:lat2,:,it)
                NCDatasets.load!(ds1["U"].var,utmp1,lon1:lon2+1,lat1:lat2,:,it)
                NCDatasets.load!(ds1["V"].var,vtmp1,lon1:lon2,lat1:lat2+1,:,it)
                NCDatasets.load!(ds1["PSFC"].var,ps1,lon1:lon2,lat1:lat2,it)
                NCDatasets.load!(ds1["U10"].var,us1,lon1:lon2,lat1:lat2,it)
                NCDatasets.load!(ds1["V10"].var,vs1,lon1:lon2,lat1:lat2,it)

                if it < 24
                    NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,it+1)
                    NCDatasets.load!(ds1["P"].var,p2,lon1:lon2,lat1:lat2,:,it+1)
                    NCDatasets.load!(ds1["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,it+1)
                    NCDatasets.load!(ds1["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,it+1)

                    NCDatasets.load!(ds1["PSFC"].var,ps2,lon1:lon2,lat1:lat2,it+1)
                    NCDatasets.load!(ds1["U10"].var,us2,lon1:lon2,lat1:lat2,it+1)
                    NCDatasets.load!(ds1["V10"].var,vs2,lon1:lon2,lat1:lat2,it+1)
                else
                    NCDatasets.load!(ds2["$(iso)QVAPOR"].var,q2,lon1:lon2,lat1:lat2,:,1)
                    NCDatasets.load!(ds2["P"].var,p2,lon1:lon2,lat1:lat2,:,1)
                    NCDatasets.load!(ds2["U"].var,utmp2,lon1:lon2+1,lat1:lat2,:,1)
                    NCDatasets.load!(ds2["V"].var,vtmp2,lon1:lon2,lat1:lat2+1,:,1)
                    
                    NCDatasets.load!(ds2["PSFC"].var,ps2,lon1:lon2,lat1:lat2,1)
                    NCDatasets.load!(ds2["U10"].var,us2,lon1:lon2,lat1:lat2,1)
                    NCDatasets.load!(ds2["V10"].var,vs2,lon1:lon2,lat1:lat2,1)
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
                    μp[ilvl+1] = (mean(p1[:,:,ilvl],wgtv) + mean(p2[:,:,ilvl],wgtv)) / 2 .+ 
                                mean(pb[:,:,ilvl],wgtv)
                    
                    Δqu[ilvl+1] = ((mean(q1[end,:,ilvl],lon2w) + 
                                    mean(q2[end,:,ilvl],lon2w)) * arc2 -
                                (mean(q1[1,:,ilvl],lon1w) + 
                                    mean(q2[1,:,ilvl],lon1w))   * arc1) / 2
                    Δqv[ilvl+1] = ((mean(q1[:,end,ilvl],lat2w) + 
                                    mean(q2[:,end,ilvl],lat2w)) * arc4 -
                                (mean(q1[:,1,ilvl],lat1w) + 
                                    mean(q2[:,1,ilvl],lat1w))   * arc3) / 2
                    Δu[ilvl+1] = ((mean(u1[end,:,ilvl],lon2w) + 
                                mean(u2[end,:,ilvl],lon2w)) * arc2 -
                                (mean(u1[1,:,ilvl],lon1w) + 
                                mean(u2[1,:,ilvl],lon1w))   * arc1) / 2
                    Δv[ilvl+1] = ((mean(v1[:,end,ilvl],lat2w) + 
                                mean(v2[:,end,ilvl],lat2w)) * arc4 -
                                (mean(v1[:,1,ilvl],lat1w) + 
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

        mkpath(datadir("wrf2","processed"))

        ds = NCDataset(fnc,"c")
        ds.dim["date"] = ndt * 24

        nctime = defVar(ds,"time",Float32,("date",),attrib=Dict(
            "units"     => "hours since $(start) 00:00:00.0",
            "long_name" => "time",
            "calendar"  => "gregorian"
        ))

        ncqdiv = defVar(ds,"$(iso)DIV",Float32,("date",),attrib=Dict(
            "units" => "kg m**-2 s**-1",
            "long_name" => "Divergence component of ∇"
        ))

        ncqadv = defVar(ds,"$(iso)ADV",Float32,("date",),attrib=Dict(
            "units" => "kg m**-2 s**-1",
            "long_name" => "Advection component of ∇"
        ))

        nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
        ncqdiv[:] = qdiv[:] * 4 / ((arc2+arc4)*(arc1+arc3)) / 9.81
        ncqadv[:] = qadv[:] * 4 / ((arc2+arc4)*(arc1+arc3)) / 9.81

        close(ds)

    else
        @warn "$(now()) - ConvectionIsotopes - Skipping $(iso)QVAPOR calculation for $(geo.name)" 
    end

end

function wrfqdivvsiwt(;
    iso   :: AbstractString,
    geo   :: GeoRegion,
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf2","grid.nc"))
    lon  = ds["longitude"][:,:,1]
    lat  = ds["latitude"][:,:,1]
    close(ds)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

    dtvec = start : Day(1) : stop

    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlvl = 50
    ndt  = length(dtvec)

    if iso != ""; iso = "$(iso)_" end

    uarr = zeros(Float32,nlon+1,nlat,nlvl)
    varr = zeros(Float32,nlon,nlat+1,nlvl)
    qarr = zeros(Float32,nlon,nlat,nlvl)
    parr = zeros(Float32,nlon,nlat,nlvl)
    psfc = zeros(Float32,nlon,nlat)
    qflx_u = zeros(nlon,nlat,24,ndt)
    qflx_v = zeros(nlon,nlat,24,ndt)
    qflxwu = zeros(Float32,nlon,nlat,24,ndt)
    qflxwv = zeros(Float32,nlon,nlat,24,ndt)

    pds = NCDataset(datadir("wrf2","raw","$(dtvec[1]).nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        iiqflxu = @view qflxwu[:,:,:,idt]
        iiqflxv = @view qflxwv[:,:,:,idt]

        ds = NCDataset(datadir("wrf2","raw","$(dtvec[idt]).nc"))

        NCDatasets.load!(ds["$(iso)IWTX"].var,iiqflxu,lon1:lon2,lat1:lat2,:)

        for ii = 1 : 24

            NCDatasets.load!(ds["U"].var,uarr,lon1:(lon2+1),lat1:lat2,:,ii)
            NCDatasets.load!(ds["V"].var,varr,lon1:lon2,lat1:(lat2+1),:,ii)
            NCDatasets.load!(ds["$(iso)QVAPOR"].var,qarr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds["P"].var,parr,lon1:lon2,lat1:lat2,:,ii)
            NCDatasets.load!(ds["PSFC"].var,psfc,lon1:lon2,lat1:lat2,ii)

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

        close(ds)

    end

    mkpath(datadir("wrf2","processed"))
    fnc = datadir("wrf2","processed","$(geo.ID)-$(iso)IWT_wrfvscalc.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["date"]      = ndt * 24

    nclon = defVar(ds,"longitude",Float32,("longitude","latitude"),attrib=Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("longitude","latitude"),attrib=Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctime = defVar(ds,"time",Float32,("date",),attrib=Dict(
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
    nctime.var[:] = collect(0 : (ndt*24 -1))
    ncvar_iwtxwrf[:] = reshape(qflxwu,nlon,nlat,:)
    ncvar_iwtywrf[:] = reshape(qflxwv,nlon,nlat,:)
    ncvar_iwtxclc[:] = reshape(qflx_u,nlon,nlat,:)
    ncvar_iwtyclc[:] = reshape(qflx_v,nlon,nlat,:)

    close(ds)

end
