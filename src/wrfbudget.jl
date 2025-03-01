using Distances
using GeoRegions
using RegionGrids
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

    wgts = ones(nlon,nlat)
    wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
    wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
    wgtm = sum(wgts)
    wgt1 = sum(wgts[1,:])
    wgt2 = sum(wgts[end,:])
    wgt3 = sum(wgts[:,1])
    wgt4 = sum(wgts[:,end])
    wgtv = weights(wgts)

    if iso == "H2O"; iso = "" end
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
    
    prcp = zeros(24,ndt) * NaN
    evap = zeros(24,ndt) * NaN
    tcwv = zeros(24,ndt) * NaN
    qflx = zeros(24,ndt) * NaN

    arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
    arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
    arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
    arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        fnc1 = datadir("wrf3","raw","$(dtvec[idt]).nc")
        if isfile(fnc1)

            ds1 = NCDataset(fnc1)

            if ds1.dim["Time"] == 24
                
                fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1))-e.nc")
                if !isfile(fnc2)
                    fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1)).nc")
                else
                    @info "$(now()) - ConvectionIsotopes - Tail end"
                end
                ds2 = NCDataset(fnc2)

                NCDatasets.load!(ds1["$(iso)RAINNC"].var,tmp1,lon1:lon2,lat1:lat2,:)
                NCDatasets.load!(ds2["$(iso)RAINNC"].var,tmp2,lon1:lon2,lat1:lat2,1)
                for ii = 1 : 24, ilat = 1 : nlat, ilon = 1 : nlon
                    tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
                end
                tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
                tmp4 = mean(tmp2,wgtv)
                prcp[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

                NCDatasets.load!(ds1["$(iso)SFCEVP"].var,tmp1,lon1:lon2,lat1:lat2,:)
                NCDatasets.load!(ds2["$(iso)SFCEVP"].var,tmp2,lon1:lon2,lat1:lat2,1)
                for ii = 1 : 24, ilat = 1 : nlat, ilon = 1 : nlon
                    tmp1[ilon,ilat,ii] *= wgts[ilon,ilat]
                end
                tmp3 = dropdims(sum(tmp1,dims=(1,2)),dims=(1,2)) / wgtm
                tmp4 = mean(tmp2,wgtv)
                evap[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

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

        end

    end

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)QBUDGET-$(dtbegstr)_$(dtbegend).nc")
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
        "long_name" => "Accumulated hourly Precipitation"
    ))

    ncevap = defVar(ds,"$(iso)E",Float32,("date",),attrib=Dict(
        "units" => "kg m**-2",
        "long_name" => "Accumulated hourly Evaporation"
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
    gvec  :: Vector{GeoRegion};
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date,
    overwrite :: Bool = true
)

    if iso == "H2O"; iso = "" end
    if iso != ""; iso = "$(iso)_" end
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]; nlon,nlat = size(lon)
    lat  = ds["latitude"][:,:]
    pb   = ds["pressure_base"][:,:,:]
    close(ds)
    nlvl = 50

    ngeo = length(gvec)
    lon1vec = zeros(Int64,ngeo)
    lat1vec = zeros(Int64,ngeo)
    lon2vec = zeros(Int64,ngeo)
    lat2vec = zeros(Int64,ngeo)
    for igeo = 1 : ngeo
        ggrd = RegionGrid(gvec[igeo],Point2.(lon,lat))
        lon1vec[igeo] = minimum(ggrd.ilon); lon2vec[igeo] = maximum(ggrd.ilon)
        lat1vec[igeo] = minimum(ggrd.ilat); lat2vec[igeo] = maximum(ggrd.ilat)
    end

    dtvec = start : Day(1) : stop
    ndt  = length(dtvec)

    utmp1 = zeros(Float32,nlon+1,nlat  ,nlvl); utmp2 = zeros(Float32,nlon+1,nlat  ,nlvl)
    vtmp1 = zeros(Float32,nlon  ,nlat+1,nlvl); vtmp2 = zeros(Float32,nlon  ,nlat+1,nlvl)

    u1 = zeros(Float32,nlon,nlat,nlvl); u2 = zeros(Float32,nlon,nlat,nlvl)
    v1 = zeros(Float32,nlon,nlat,nlvl); v2 = zeros(Float32,nlon,nlat,nlvl)
    q1 = zeros(Float32,nlon,nlat,nlvl); q2 = zeros(Float32,nlon,nlat,nlvl)
    p1 = zeros(Float32,nlon,nlat,nlvl); p2 = zeros(Float32,nlon,nlat,nlvl)

    us1 = zeros(Float32,nlon,nlat); us2 = zeros(Float32,nlon,nlat)
    vs1 = zeros(Float32,nlon,nlat); vs2 = zeros(Float32,nlon,nlat)
    ps1 = zeros(Float32,nlon,nlat); ps2 = zeros(Float32,nlon,nlat)

    u = zeros(Float32,nlon,nlat,nlvl+2)
    v = zeros(Float32,nlon,nlat,nlvl+2)
    q = zeros(Float32,nlon,nlat,nlvl+2)
    p = zeros(Float32,nlon,nlat,nlvl+2)

    ∇ = zeros(Float32,24,ndt,ngeo)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting $(iso)QVAPOR data during $(dtvec[idt])"
        flush(stderr)

        fnc1 = datadir("wrf3","raw","$(dtvec[idt]).nc")
        if isfile(fnc1)

            ds1 = NCDataset(fnc1)

            if ds1.dim["Time"] == 24
                
                fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1))-e.nc")
                if !isfile(fnc2)
                    fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1)).nc")
                else
                    @info "$(now()) - ConvectionIsotopes - Tail end"
                end
                ds2 = NCDataset(fnc2)

                for it = 1 : 24

                    NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q1,:,:,:,it)
                    NCDatasets.load!(ds1["P"].var,p1,:,:,:,it)
                    NCDatasets.load!(ds1["U"].var,utmp1,:,:,:,it)
                    NCDatasets.load!(ds1["V"].var,vtmp1,:,:,:,it)
                    NCDatasets.load!(ds1["PSFC"].var,ps1,:,:,it)
                    NCDatasets.load!(ds1["U10"].var,us1,:,:,it)
                    NCDatasets.load!(ds1["V10"].var,vs1,:,:,it)

                    if it < 24
                        NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q2,:,:,:,it+1)
                        NCDatasets.load!(ds1["P"].var,p2,:,:,:,it+1)
                        NCDatasets.load!(ds1["U"].var,utmp2,:,:,:,it+1)
                        NCDatasets.load!(ds1["V"].var,vtmp2,:,:,:,it+1)
                        NCDatasets.load!(ds1["PSFC"].var,ps2,:,:,it+1)
                        NCDatasets.load!(ds1["U10"].var,us2,:,:,it+1)
                        NCDatasets.load!(ds1["V10"].var,vs2,:,:,it+1)
                    else
                        NCDatasets.load!(ds2["$(iso)QVAPOR"].var,q2,:,:,:,1)
                        NCDatasets.load!(ds2["P"].var,p2,:,:,:,1)
                        NCDatasets.load!(ds2["U"].var,utmp2,:,:,:,1)
                        NCDatasets.load!(ds2["V"].var,vtmp2,:,:,:,1)
                        NCDatasets.load!(ds2["PSFC"].var,ps2,:,:,1)
                        NCDatasets.load!(ds2["U10"].var,us2,:,:,1)
                        NCDatasets.load!(ds2["V10"].var,vs2,:,:,1)
                    end

                    @views @. u1 = (utmp1[1:(end-1),:,:] + utmp1[2:end,:,:]) / 2
                    @views @. u2 = (utmp2[1:(end-1),:,:] + utmp2[2:end,:,:]) / 2
                    @views @. v1 = (vtmp1[:,1:(end-1),:] + vtmp1[:,2:end,:]) / 2
                    @views @. v2 = (vtmp2[:,1:(end-1),:] + vtmp2[:,2:end,:]) / 2

                    @views @. u[:,:,2:(nlvl+1)] = (u1 + u2) / 2
                    @views @. v[:,:,2:(nlvl+1)] = (v1 + v2) / 2
                    @views @. q[:,:,2:(nlvl+1)] = (q1 + q2) / 2
                    @views @. p[:,:,2:(nlvl+1)] = (p1 + p2) / 2 + pb

                    @views @. u[:,:,1] = (us1 + us2) / 2
                    @views @. v[:,:,1] = (vs1 + vs2) / 2
                    @views @. p[:,:,1] = (ps1 + ps2) / 2
                    @views @. q[:,:,1] = q[:,:,2]

                    for igeo in 1 : ngeo

                        lon1 = lon1vec[igeo]; lon2 = lon2vec[igeo]; nilon = lon2 - lon1 + 1
                        lat1 = lat1vec[igeo]; lat2 = lat2vec[igeo]; nilat = lat2 - lat1 + 1

                        arc1 = haversine(
                            (lon[lon1,lat1],lat[lon1,lat1]),
                            (lon[lon1,lat2],lat[lon1,lat2])
                        )
                        arc2 = haversine(
                            (lon[lon2,lat1],lat[lon2,lat1]),
                            (lon[lon2,lat2],lat[lon2,lat2])
                        )
                        arc3 = haversine(
                            (lon[lon1,lat1],lat[lon1,lat1]),
                            (lon[lon2,lat1],lat[lon2,lat1])
                        )
                        arc4 = haversine(
                            (lon[lon1,lat2],lat[lon1,lat2]),
                            (lon[lon2,lat2],lat[lon2,lat2])
                        )

                        for ilat = (lat1+1) : (lat2-1)
                            ∇[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]),
                                reverse(q[lon1,ilat,:] .* u[lon1,ilat,:])
                            ) / (nilat-1) * arc1
                            ∇[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]),
                                reverse(q[lon2,ilat,:] .* u[lon2,ilat,:])
                            ) / (nilat-1) * arc2
                        end
    
                        for ilon = (lon1+1) : (lon2-1)
                            ∇[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]),
                                reverse(q[ilon,lat1,:] .* v[ilon,lat1,:])
                            ) / (nilon-1) * arc3
                            ∇[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]),
                                reverse(q[ilon,lat2,:] .* v[ilon,lat2,:])
                            ) / (nilon-1) * arc4
                        end
    
                        for ilat in [lat1, lat2]
                            ∇[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]),
                                reverse(q[lon1,ilat,:] .* u[lon1,ilat,:])
                            ) / (nilat-1) * arc1 / 2
                            ∇[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]),
                                reverse(q[lon2,ilat,:] .* u[lon2,ilat,:])
                            ) / (nilat-1) * arc2 / 2
                        end
    
                        for ilon in [lon1, lon2]
                            ∇[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]),
                                reverse(q[ilon,lat1,:] .* v[ilon,lat1,:])
                            ) / (nilon-1) * arc3 / 2
                            ∇[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]),
                                reverse(q[ilon,lat2,:] .* v[ilon,lat2,:])
                            ) / (nilon-1) * arc4 / 2
                        end

                        ∇[it,idt,igeo] *= (4 / ((arc1+arc2)*(arc3+arc4)) / 9.81)

                    end
                
                end
                close(ds2)

            end

            close(ds1)

        end

    end


    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    for igeo = 1 : ngeo
        fnc = datadir("wrf3","processed","$(gvec[igeo].ID)-$(iso)∇-$(dtbegstr)_$(dtbegend).nc")
        if !isfile(fnc) || overwrite
            rm(fnc,force=true)
        end
        mkpath(datadir("wrf3","processed"))

        ds = NCDataset(fnc,"c")
        ds.dim["date"] = ndt * 24

        nctime = defVar(ds,"time",Float32,("date",),attrib=Dict(
            "units"     => "hours since $(start) 00:00:00.0",
            "long_name" => "time",
            "calendar"  => "gregorian"
        ))

        ncqdiv = defVar(ds,"$(iso)∇",Float32,("date",),attrib=Dict(
            "units" => "kg m**-2 s**-1",
            "long_name" => "Divergence component of ∇"
        ))

        nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
        ncqdiv[:] = ∇[:,:,igeo][:]

        close(ds)
    end

end

function wrfqdivdecompose(
    gvec  :: Vector{GeoRegion};
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date,
    overwrite :: Bool = true
)

    if iso == "H2O"; iso = "" end
    if iso != ""; iso = "$(iso)_" end
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    lon  = ds["longitude"][:,:]; nlon,nlat = size(lon)
    lat  = ds["latitude"][:,:]
    pb   = ds["pressure_base"][:,:,:]
    close(ds)
    nlvl = 50

    ngeo = length(gvec)
    lon1vec = zeros(Int64,ngeo)
    lat1vec = zeros(Int64,ngeo)
    lon2vec = zeros(Int64,ngeo)
    lat2vec = zeros(Int64,ngeo)
    for igeo = 1 : ngeo
        ggrd = RegionGrid(gvec[igeo],Point2.(lon,lat))
        lon1vec[igeo] = minimum(ggrd.ilon); lon2vec[igeo] = maximum(ggrd.ilon)
        lat1vec[igeo] = minimum(ggrd.ilat); lat2vec[igeo] = maximum(ggrd.ilat)
    end

    dtvec = start : Day(1) : stop
    ndt  = length(dtvec)

    utmp1 = zeros(Float32,nlon+1,nlat  ,nlvl); utmp2 = zeros(Float32,nlon+1,nlat  ,nlvl)
    vtmp1 = zeros(Float32,nlon  ,nlat+1,nlvl); vtmp2 = zeros(Float32,nlon  ,nlat+1,nlvl)

    u1 = zeros(Float32,nlon,nlat,nlvl); u2 = zeros(Float32,nlon,nlat,nlvl)
    v1 = zeros(Float32,nlon,nlat,nlvl); v2 = zeros(Float32,nlon,nlat,nlvl)
    q1 = zeros(Float32,nlon,nlat,nlvl); q2 = zeros(Float32,nlon,nlat,nlvl)
    p1 = zeros(Float32,nlon,nlat,nlvl); p2 = zeros(Float32,nlon,nlat,nlvl)

    us1 = zeros(Float32,nlon,nlat); us2 = zeros(Float32,nlon,nlat)
    vs1 = zeros(Float32,nlon,nlat); vs2 = zeros(Float32,nlon,nlat)
    ps1 = zeros(Float32,nlon,nlat); ps2 = zeros(Float32,nlon,nlat)

    u = zeros(Float32,nlon,nlat,nlvl+2)
    v = zeros(Float32,nlon,nlat,nlvl+2)
    q = zeros(Float32,nlon,nlat,nlvl+2)
    p = zeros(Float32,nlon,nlat,nlvl+2)

    μq = zeros(Float32,nlvl+2)
    μu = zeros(Float32,nlvl+2)
    μv = zeros(Float32,nlvl+2)
    μp = zeros(Float32,nlvl+2)
    
    qadv = zeros(Float32,24,ndt,ngeo) * NaN
    qdiv = zeros(Float32,24,ndt,ngeo) * NaN

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting $(iso)QVAPOR data during $(dtvec[idt])"
        flush(stderr)

        fnc1 = datadir("wrf3","raw","$(dtvec[idt]).nc")
        if isfile(fnc1)

            ds1 = NCDataset(fnc1)

            if ds1.dim["Time"] == 24
                
                fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1))-e.nc")
                if !isfile(fnc2)
                    fnc2 = datadir("wrf3","raw","$(dtvec[idt]+Day(1)).nc")
                else
                    @info "$(now()) - ConvectionIsotopes - Tail end"
                end
                ds2 = NCDataset(fnc2)

                for it = 1 : 24

                    NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q1,:,:,:,it)
                    NCDatasets.load!(ds1["P"].var,p1,:,:,:,it)
                    NCDatasets.load!(ds1["U"].var,utmp1,:,:,:,it)
                    NCDatasets.load!(ds1["V"].var,vtmp1,:,:,:,it)
                    NCDatasets.load!(ds1["PSFC"].var,ps1,:,:,it)
                    NCDatasets.load!(ds1["U10"].var,us1,:,:,it)
                    NCDatasets.load!(ds1["V10"].var,vs1,:,:,it)

                    if it < 24
                        NCDatasets.load!(ds1["$(iso)QVAPOR"].var,q2,:,:,:,it+1)
                        NCDatasets.load!(ds1["P"].var,p2,:,:,:,it+1)
                        NCDatasets.load!(ds1["U"].var,utmp2,:,:,:,it+1)
                        NCDatasets.load!(ds1["V"].var,vtmp2,:,:,:,it+1)
                        NCDatasets.load!(ds1["PSFC"].var,ps2,:,:,it+1)
                        NCDatasets.load!(ds1["U10"].var,us2,:,:,it+1)
                        NCDatasets.load!(ds1["V10"].var,vs2,:,:,it+1)
                    else
                        NCDatasets.load!(ds2["$(iso)QVAPOR"].var,q2,:,:,:,1)
                        NCDatasets.load!(ds2["P"].var,p2,:,:,:,1)
                        NCDatasets.load!(ds2["U"].var,utmp2,:,:,:,1)
                        NCDatasets.load!(ds2["V"].var,vtmp2,:,:,:,1)
                        NCDatasets.load!(ds2["PSFC"].var,ps2,:,:,1)
                        NCDatasets.load!(ds2["U10"].var,us2,:,:,1)
                        NCDatasets.load!(ds2["V10"].var,vs2,:,:,1)
                    end

                    @views @. u1 = (utmp1[1:(end-1),:,:] + utmp1[2:end,:,:]) / 2
                    @views @. u2 = (utmp2[1:(end-1),:,:] + utmp2[2:end,:,:]) / 2
                    @views @. v1 = (vtmp1[:,1:(end-1),:] + vtmp1[:,2:end,:]) / 2
                    @views @. v2 = (vtmp2[:,1:(end-1),:] + vtmp2[:,2:end,:]) / 2

                    @views @. u[:,:,2:(nlvl+1)] = (u1 + u2) / 2
                    @views @. v[:,:,2:(nlvl+1)] = (v1 + v2) / 2
                    @views @. q[:,:,2:(nlvl+1)] = (q1 + q2) / 2
                    @views @. p[:,:,2:(nlvl+1)] = (p1 + p2) / 2 + pb

                    @views @. u[:,:,1] = (us1 + us2) / 2
                    @views @. v[:,:,1] = (vs1 + vs2) / 2
                    @views @. p[:,:,1] = (ps1 + ps2) / 2
                    @views @. q[:,:,1] = q[:,:,2]

                    for igeo in 1 : ngeo

                        lon1 = lon1vec[igeo]; lon2 = lon2vec[igeo]; lonr = lon1 : lon2
                        lat1 = lat1vec[igeo]; lat2 = lat2vec[igeo]; latr = lat1 : lat2

                        nglon = lon2 - lon1 + 1
                        nglat = lat2 - lat1 + 1

                        wgts = ones(nglon,nglat)
                        wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
                        wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
                        wgtv  = weights(wgts)

                        arc1 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon1,lat2],lat[lon1,lat2]))
                        arc2 = haversine((lon[lon2,lat1],lat[lon2,lat1]),(lon[lon2,lat2],lat[lon2,lat2]))
                        arc3 = haversine((lon[lon1,lat1],lat[lon1,lat1]),(lon[lon2,lat1],lat[lon2,lat1]))
                        arc4 = haversine((lon[lon1,lat2],lat[lon1,lat2]),(lon[lon2,lat2],lat[lon2,lat2]))

                        for ilvl = 1 : (nlvl+1)

                            μq[ilvl] = mean(view(q,lonr,latr,ilvl),wgtv)
                            μu[ilvl] = mean(view(u,lonr,latr,ilvl),wgtv)
                            μv[ilvl] = mean(view(v,lonr,latr,ilvl),wgtv)
                            μp[ilvl] = mean(view(p,lonr,latr,ilvl),wgtv)
                            
                        end

                        qadv[it,idt,igeo] = 0
                        qdiv[it,idt,igeo] = 0

                        for ilat = (lat1+1) : (lat2-1)
                            qadv[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]), reverse(q[lon1,ilat,:] .* μu)
                            ) / (nglat-1) * arc1
                            qadv[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]), reverse(q[lon2,ilat,:] .* μu)
                            ) / (nglat-1) * arc2
                            qdiv[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]), reverse(u[lon1,ilat,:] .* μq)
                            ) / (nglat-1) * arc1
                            qdiv[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]), reverse(u[lon2,ilat,:] .* μq)
                            ) / (nglat-1) * arc2
                        end
    
                        for ilon = (lon1+1) : (lon2-1)
                            qadv[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]), reverse(q[ilon,lat1,:] .* μv)
                            ) / (nglon-1) * arc3
                            qadv[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]), reverse(q[ilon,lat2,:] .* μv)
                            ) / (nglon-1) * arc4
                            qdiv[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]), reverse(v[ilon,lat1,:] .* μq)
                            ) / (nglon-1) * arc3
                            qdiv[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]), reverse(v[ilon,lat2,:] .* μq)
                            ) / (nglon-1) * arc4
                        end
    
                        for ilat in [lat1, lat2]
                            qadv[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]), reverse(q[lon1,ilat,:] .* μu)
                            ) / (nglat-1) * arc1 / 2
                            qadv[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]), reverse(q[lon2,ilat,:] .* μu)
                            ) / (nglat-1) * arc2 / 2
                            qdiv[it,idt,igeo] -= trapz(
                                reverse(p[lon1,ilat,:]), reverse(u[lon1,ilat,:] .* μq)
                            ) / (nglat-1) * arc1 / 2
                            qdiv[it,idt,igeo] += trapz(
                                reverse(p[lon2,ilat,:]), reverse(u[lon2,ilat,:] .* μq)
                            ) / (nglat-1) * arc2 / 2
                        end
    
                        for ilon in [lon1, lon2]
                            qadv[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]), reverse(q[ilon,lat1,:] .* μv)
                            ) / (nglon-1) * arc3 / 2
                            qadv[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]), reverse(q[ilon,lat2,:] .* μv)
                            ) / (nglon-1) * arc4 / 2
                            qdiv[it,idt,igeo] -= trapz(
                                reverse(p[ilon,lat1,:]), reverse(v[ilon,lat1,:] .* μq)
                            ) / (nglon-1) * arc3 / 2
                            qdiv[it,idt,igeo] += trapz(
                                reverse(p[ilon,lat2,:]), reverse(v[ilon,lat2,:] .* μq)
                            ) / (nglon-1) * arc4 / 2
                        end

                        qadv[it,idt,igeo] *= (4 / ((arc1+arc2)*(arc3+arc4)) / 9.81)
                        qdiv[it,idt,igeo] *= (4 / ((arc1+arc2)*(arc3+arc4)) / 9.81)

                    end
                
                end
                close(ds2)

            end

            close(ds1)

        end

    end


    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    for igeo = 1 : ngeo
        fnc = datadir("wrf3","processed","$(gvec[igeo].ID)-$(iso)∇decompose-$(dtbegstr)_$(dtbegend).nc")
        if !isfile(fnc) || overwrite
            rm(fnc,force=true)
        end
        mkpath(datadir("wrf3","processed"))

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
        ncqdiv[:] = qdiv[:,:,igeo][:]
        ncqadv[:] = qadv[:,:,igeo][:]

        close(ds)
    end

end

function wrfqdivdecomposecompile(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
    start :: Date,
    stop  :: Date,
    overwrite :: Bool = true
)

    if iso == "H2O"; iso = "" end
    if iso != ""; iso = "$(iso)_" end

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    dtvec = start : Day(1) : stop; ndt = length(dtvec) * 24

    DIV = []
    ADV = []
    for idt = start : Month(1) : stop
        iyr = year(idt)
        imo = month(idt)
        ndy = daysinmonth(iyr,imo)
        iidtbegstr = Dates.format(Date(iyr,imo,1),dateformat"yyyymmdd")
        iidtbegend = Dates.format(Date(iyr,imo,ndy),dateformat"yyyymmdd")
        fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)∇decompose-$(iidtbegstr)_$(iidtbegend).nc")
        ds = NCDataset(fnc)
        DIV = vcat(DIV,ds["$(iso)DIV"])
        ADV = vcat(ADV,ds["$(iso)ADV"])
        close(ds)

    end

    fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)∇decompose-$(dtbegstr)_$(dtbegend).nc")
    if !isfile(fnc) || overwrite
        rm(fnc,force=true)
    end
    mkpath(datadir("wrf3","processed"))

    ds = NCDataset(fnc,"c")
    ds.dim["date"] = ndt

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

    nctime.var[:] = collect(0 : (ndt-1)) .+ 0.5
    ncqdiv[:] = DIV
    ncqadv[:] = ADV

    close(ds)

end

function wrfqdivvsiwt(;
    iso   :: AbstractString,
    geo   :: GeoRegion,
    start :: Date,
    stop  :: Date
)
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
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

    pds = NCDataset(datadir("wrf3","raw","$(dtvec[1]).nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    for idt in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
        flush(stderr)

        iiqflxu = @view qflxwu[:,:,:,idt]
        iiqflxv = @view qflxwv[:,:,:,idt]

        ds = NCDataset(datadir("wrf3","raw","$(dtvec[idt]).nc"))

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

    mkpath(datadir("wrf3","processed"))
    fnc = datadir("wrf3","processed","$(geo.ID)-$(iso)IWT_wrfvscalc.nc")
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
