using GeoRegions
using RegionGrids
using NCDatasets
using StatsBase

include(srcdir("backend.jl"))

# function wrfdqdp(
#     geo   :: GeoRegion;
#     iso   :: AbstractString = "",
#     start :: Date,
#     stop  :: Date,
# 	days  :: Int = 0
# )
    
#     ds   = NCDataset(datadir("wrf","grid.nc"))
#     lon  = ds["longitude"][:,:]
#     lat  = ds["latitude"][:,:]
#     close(ds)

#     ggrd = RegionGrid(geo,lon,lat)
#     apnt = findall(ggrd.mask .== 1)
#     npnt = length(apnt)
#     lon1 = findfirst(ggrd.mask .== 1)[1]; lon2 = findfirst(ggrd.mask .== 1)[1]
#     lat1 = findfirst(ggrd.mask .== 1)[2]; lat2 = findfirst(ggrd.mask .== 1)[2]

#     for ipnt = 2 : npnt
#         ilon = apnt[ipnt][1]; ilat = apnt[ipnt][2]
#         if ilon < lon1; lon1 = ilon; end
#         if ilon > lon2; lon2 = ilon; end
#         if ilat < lat1; lat1 = ilat; end
#         if ilat > lat2; lat2 = ilat; end
#     end

#     dtvec = start : Day(1) : stop

#     nlon = lon2 - lon1 + 1
#     nlat = lat2 - lat1 + 1
#     nlvl = 50
#     ndt  = length(dtvec)

#     wgts = ones(nlon,nlat)
#     wgts[1,:] *= 0.5; wgts[end,:] *= 0.5
#     wgts[:,1] *= 0.5; wgts[:,end] *= 0.5
#     wgtm = sum(wgts)

#     if iso != ""; iso = "$(iso)_" end
# 	if !iszero(days); days = "-smooth_$(@sprintf("%02d",days))days"; else; days = "" end

#     tmpq = zeros(Float32,nlon,nlat,nlvl)
#     tmpp = zeros(Float32,nlon,nlat,nlvl)
    
#     pvec = zeros(nlvl,ndt)
#     dqdp = zeros(nlvl,ndt)

#     dsp = NCDataset(datadir("wrf","3D","PB-daily.nc"))
#     pbs = dsp["PB"][lon1:lon2,lat1:lat2,:,1]
#     close(dsp)
#     pbs = dropdims(sum(pbs .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm

# 	dsq = NCDataset(datadir("wrf","3D","$(iso)QVAPOR-daily$(days).nc"))
# 	dsp = NCDataset(datadir("wrf","3D","P-daily$(days).nc"))

#     for idt in 1 : ndt

#         @info "$(now()) - ConvectionIsotopes - Extracting data for $(dtvec[idt])"
#         flush(stderr)

#         NCDatasets.load!(dsq["$(iso)QVAPOR"].var,tmpq,lon1:lon2,lat1:lat2,:,idt)
#         NCDatasets.load!(dsp["P"].var,tmpp,lon1:lon2,lat1:lat2,:,idt)

# 		q            = dropdims(sum(tmpq .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm
# 		pvec[:,idt] .= dropdims(sum(tmpp .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm .+ pbs

# 		for ilvl = 2 : (nlvl-1)
# 			dqdp[ilvl,idt] = (q[ilvl+1]-q[ilvl-1]) ./ (pvec[ilvl+1,idt]-pvec[ilvl-1,idt])
# 		end
# 		dqdp[1,idt] = (q[2]-q[1]) / (pvec[2,idt]-pvec[1,idt])
# 		dqdp[end,idt] = (q[end]-q[end-1]) / (pvec[end,idt]-pvec[end-1,idt])

#     end

# 	close(dsq)
# 	close(dsp)

#     mkpath(datadir("wrf","processed"))
#     fnc = datadir("wrf","processed","$(geo.ID)-$(iso)dqdp$(days).nc")
#     if isfile(fnc); rm(fnc,force=true) end

#     ds = NCDataset(fnc,"c")
#     ds.dim["level"] = nlvl
#     ds.dim["date"]  = ndt

#     nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
#         "units"     => "hours since $(start) 00:00:00.0",
#         "long_name" => "time",
#         "calendar"  => "gregorian"
#     ))

#     ncpres = defVar(ds,"P",Float64,("level","date",),attrib=Dict(
#         "units" => "Pa",
#         "long_name" => "Pressure"
#     ))

#     ncdqdp = defVar(ds,"$(iso)dqdp",Float64,("level","date",),attrib=Dict(
#         "units" => "kg kg**-1 Pa**-1",
#         "long_name" => "Gradient of $(iso)VAPOR against pressure"
#     ))

#     nctime.var[:] = collect(0 : (ndt-1))
#     ncpres[:,:] = pvec
#     ncdqdp[:,:] = dqdp

#     close(ds)

# end

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
    pbs = ds["pressure_base"][:,:]
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

    tmpq = zeros(Float32,nlon,nlat,nlvl)
    tmph = zeros(Float32,nlon,nlat,nlvl)
    tmpp = zeros(Float32,nlon,nlat,nlvl)
    
    pvec  = zeros(Float32,nlvl,ndt)
    hq    = zeros(Float32,nlvl,ndt)
    dhqdp = zeros(Float32,nlvl,ndt)

    pbs = dropdims(sum(pbs[lon1:lon2,lat1:lat2] .* wgts,dims=(1,2)),dims=(1,2)) ./ wgtm

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