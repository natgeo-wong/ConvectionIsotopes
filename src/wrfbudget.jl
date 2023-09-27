using Distances
using GeoRegions
using NCDatasets
using Statistics
using Trapz

include(srcdir("backend.jl"))

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

            for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
                parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
            end

            for ilat = 2 : nlat
                qavg = dropdims(mean(qarr[1:2,ilat,:],dims=1),dims=1)
                qlat = uarr[1,ilat,:] .* qavg
                plat = parr[1,ilat,:]
                qflx[ii,idt] -= trapz(vcat(psfc[1,ilat],plat,0),vcat(qlat[1],qlat,0))
                qavg = dropdims(mean(qarr[nlon.+(0:1),ilat,:],dims=1),dims=1)
                qlat = uarr[end,ilat,:] .* qavg
                plat = parr[end,ilat,:]
                qflx[ii,idt] += trapz(vcat(psfc[end,ilat],plat,0),vcat(qlat[1],qlat,0))
            end

            for ilon = 2 : nlon
                qavg = dropdims(mean(qarr[ilon,1:2,:],dims=1),dims=1)
                qlat = varr[ilon,1,:] .* qavg
                plat = parr[ilon,1,:]
                qflx[ii,idt] -= trapz(vcat(psfc[ilon,1],plat,0),vcat(qlat[1],qlat,0))
                qavg = dropdims(mean(qarr[ilon,nlat.+(0:1),:],dims=1),dims=1)
                qlat = varr[ilon,end,:] .* qavg
                plat = parr[ilon,end,:]
                qflx[ii,idt] += trapz(vcat(psfc[ilon,end],plat,0),vcat(qlat[1],qlat,0))
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
    ncvar[:] = qflx[:] / 9.81 / 1000 * 3600 * 3 * (
        haversine((lonA,latA),(lonB,latB)) + haversine((lonB,latB),(lonC,latC)) +
        haversine((lonC,latC),(lonD,latD)) + haversine((lonD,latD),(lonA,latA))
    ) / (
        (haversine((lonA,latA),(lonB,latB)) + haversine((lonC,latC),(lonD,latD))) *
        (haversine((lonB,latB),(lonC,latC)) + haversine((lonD,latD),(lonA,latA)))
    ) * 4

    close(ds)

end

function wrfqbudget(
    geo   :: GeoRegion;
    iso   :: AbstractString = "",
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

    if iso != ""; iso = "$(iso)_" end

    tmp1 = zeros(Float32,nlon,nlat,8)
    tmp2 = zeros(Float32,nlon,nlat)
    prcp = zeros(8,ndt)
    evap = zeros(8,ndt)
    tcwv = zeros(8,ndt)

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
        tmp3 = dropdims(mean(tmp1,dims=(1,2)),dims=(1,2))
        tmp4 = mean(tmp2)
        prcp[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

        NCDatasets.load!(ds1["$(iso)QFX"].var,tmp1,lon1:lon2,lat1:lat2,:)
        evap[:,idt] = dropdims(mean(tmp1,dims=(1,2)),dims=(1,2))

        NCDatasets.load!(ds1["$(iso)VAPORWP"].var,tmp1,lon1:lon2,lat1:lat2,:)
        NCDatasets.load!(ds2["$(iso)VAPORWP"].var,tmp2,lon1:lon2,lat1:lat2,1)
        tmp3 = dropdims(mean(tmp1,dims=(1,2)),dims=(1,2))
        tmp4 = mean(tmp2)
        tcwv[:,idt] = vcat(tmp3[2:end],tmp4) .- tmp3

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

    nctime.var[:] = (collect(0 : (ndt*8 -1))) * 3
    ncprcp.var[:] = prcp[:]
    ncevap.var[:] = evap[:]
    nctcwv.var[:] = tcwv[:]

    close(ds)

end