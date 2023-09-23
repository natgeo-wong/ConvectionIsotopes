using GeoRegions
using NCDatasets
using Statistics
using Trapz

include(srcdir("backend.jl"))

function wrfbudget(
    wvar :: AbstractString,
    geo  :: GeoRegion,
)
    
    ds   = NCDataset(datadir("wrf","3D","$(wvar)-daily.nc"))
    lon  = ds["longitude"][:]
    lat  = ds["latitude"][:]
    nlvl = ds.dim["levels"]
    ndt  = ds.dim["date"]
    close(ds)

    ggrd = RegionGrid(geo,lon,lat)
    lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findlast(ggrd.mask .== 1)[1]
    lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findlast(ggrd.mask .== 1)[2]

    nlon = lon2 - lon1
    nlat = lat2 - lat1

    uarr = zeros(Float32,nlon+2,nlat+1,nlvl)
    varr = zeros(Float32,nlon+1,nlat+2,nlvl)
    qarr = zeros(Float32,nlon+1,nlat+1,nlvl)
    parr = zeros(Float32,nlon+1,nlat+1,nlvl)
    psfc = zeros(Float32,nlon+1,nlat+1,nlvl)
    qflx = zeros(ndt)


    pds = NCDataset(datadir("wrf","3D","PB-daily.nc"))
    pbse = pds["PB"][lon1:lon2,lat1:lat2,:,1]
    close(pds)

    uds = NCDataset(datadir("wrf","3D","U-daily.nc"))
    vds = NCDataset(datadir("wrf","3D","V-daily.nc"))
    qds = NCDataset(datadir("wrf","3D","$(wvar)-daily.nc"))
    pds = NCDataset(datadir("wrf","3D","P-daily.nc"))
    sds = NCDataset(datadir("wrf","3D","PSFC-daily.nc"))

    for ii in 1 : ndt

        @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[ii])"
        flush(stderr)

        NCDatasets.load!(uds["U"].var,uarr,lon1:(lon2+1),lat1:lat2,:,ii)
        NCDatasets.load!(vds["V"].var,varr,lon1:lon2,lat1:(lat2+1),:,ii)
        NCDatasets.load!(qds[wvar].var,qarr,lon1:lon2,lat1:lat2,:,ii)
        NCDatasets.load!(pds["P"].var,parr,lon1:lon2,lat1:lat2,:,ii)
        NCDatasets.load!(sds["PSFC"].var,psfc,lon1:lon2,lat1:lat2,:,ii)

        for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
            parr[ilon,ilat,ilvl] += pbse[ilon,ilat,ilvl]
        end

        for ilat = 2 : nlat
            uavg = dropdims(mean(uarr[1:2,ilat,:],dims=1),dims=1)
            qlat = qarr[1,ilat] .* uavg
            plat = parr[1,ilat]
            qflx[ii] += trapz(vcat(psfc,plat,0),vcat(0,qlat,0))
            uavg = dropdims(mean(uarr[nlon.+(1:2),ilat,:],dims=1),dims=1)
            qlat = qarr[end,ilat] .* uavg
            plat = parr[end,ilat]
            qflx[ii] -= trapz(vcat(psfc,plat,0),vcat(0,qlat,0))
        end

        for ilat = [1, nlat+1]
            uavg = dropdims(mean(uarr[1:2,ilat,:],dims=1),dims=1)
            qlat = qarr[1,ilat] .* uavg
            plat = parr[1,ilat]
            qflx[ii] += trapz(vcat(psfc,plat,0),vcat(0,qlat,0)) * 0.5
            uavg = dropdims(mean(uarr[nlon.+(1:2),ilat,:],dims=1),dims=1)
            qlat = qarr[end,ilat] .* uavg
            plat = parr[end,ilat]
            qflx[ii] -= trapz(vcat(psfc,plat,0),vcat(0,qlat,0)) * 0.5
        end

        for ilon = 2 : nlon
            vavg = dropdims(mean(uarr[ilon,1:2,:],dims=1),dims=1)
            qlat = qarr[ilon,1] .* vavg
            plat = parr[ilon,1]
            qflx[ii] += trapz(vcat(psfc,plat,0),vcat(0,qlat,0))
            vavg = dropdims(mean(uarr[ilon,nlat.+(1:2),:],dims=1),dims=1)
            qlat = qarr[ilon,end] .* vavg
            plat = parr[ilon,end]
            qflx[ii] -= trapz(vcat(psfc,plat,0),vcat(0,qlat,0))
        end

        for ilon = [1, nlon+1]
            vavg = dropdims(mean(uarr[ilon,1:2,:],dims=1),dims=1)
            qlat = qarr[ilon,1] .* vavg
            plat = parr[ilon,1]
            qflx[ii] += trapz(vcat(psfc,plat,0),vcat(0,qlat,0)) * 0.5
            vavg = dropdims(mean(uarr[ilon,nlat.+(1:2),:],dims=1),dims=1)
            qlat = qarr[ilon,end] .* vavg
            plat = parr[ilon,end]
            qflx[ii] -= trapz(vcat(psfc,plat,0),vcat(0,qlat,0)) * 0.5
        end

    end

    close(uds)
    close(vds)
    close(qds)
    close(pds)
    close(sds)

    mkpath(datadir("wrf","processed"))
    fnc = datadir("wrf","processed","FLUX_$(wvar)-daily.nc")
    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]      = ndt

    nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncvar = defVar(ds,"FLUX_$(wvar)",Float32,("date",),attrib=attrib)

    nctime.var[:] = collect(0 : (ndt-1))
    ncvar[:] = qflx

    close(ds)

end