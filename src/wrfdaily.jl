using DrWatson
using Dates
using Logging
using NCDatasets
using Statistics

function wrf3Ddaily(
    wvar :: AbstractString;
	start :: Date,
	stop  :: Date,
)

    ds  = NCDataset(datadir("wrf","raw","3D","$(start).nc"))
	nlon,nlat,nlvl = size(ds[wvar])[[1,2,3]]
	if wvar == "U"
		lon = ds["XLONG_U"][:,:,1]
		lat = ds["XLAT_U"][:,:,1]
	elseif wvar == "V"
		lon = ds["XLONG_V"][:,:,1]
		lat = ds["XLAT_V"][:,:,1]
	else
		lon = ds["XLONG"][:,:,1]
		lat = ds["XLAT"][:,:,1]
	end
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	dtvec = start : Day(1) : stop; ndt = length(dtvec)

	oarr = zeros(Float32,nlon,nlat,nlvl,8)
	narr = zeros(Float32,nlon,nlat,nlvl,ndt)

	for ii in 1 : ndt

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for $(dtvec[ii])"
		flush(stderr)

		ids = NCDataset(datadir("wrf","raw","3D","$(dtvec[ii]).nc"))
		NCDatasets.load!(ids[wvar].var,oarr,:,:,:,:)
		close(ids)

		for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
			narr[ilon,ilat,ilvl,ii] = Float32(mean(view(oarr,ilon,ilat,ilvl,:)))
		end

	end

	fnc = datadir("wrf","3D","$wvar-daily.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["levels"]    = nlvl
	ds.dim["date"]      = ndt

	nclon = defVar(ds,"longitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
		"units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
	))

	ncvar = defVar(ds,wvar,Float32,("longitude","latitude","levels","date"),attrib=attrib)

	nclon[:] = lon
	nclat[:] = lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:] = narr

	close(ds)

end

function wrf2Ddaily(
    wvar :: AbstractString;
	start :: Date,
	stop  :: Date,
	isaccum :: Bool = false
)

    ds  = NCDataset(datadir("wrf","raw","2D","$(start).nc"))
	nlon,nlat = size(ds[wvar])[[1,2]]
	lon = ds["XLONG"][:,:,1]
	lat = ds["XLAT"][:,:,1]
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	dtvec = start : Day(1) : stop; ndt = length(dtvec)

	oarr = zeros(Float32,nlon,nlat,8)
	narr = zeros(Float32,nlon,nlat,ndt)
	if isaccum
		aarr = zeros(Float32,nlon,nlat)
	end

	for ii in 1 : ndt

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for $(dtvec[ii])"
		flush(stderr)

		ids = NCDataset(datadir("wrf","raw","2D","$(dtvec[ii]).nc"))
		NCDatasets.load!(ids[wvar].var,oarr,:,:,:)
		close(ids)
		if isaccum
			fncii = datadir("wrf","raw","2D","$(dtvec[ii]+Day(1))-e.nc")
			if isfile(fncii)
				@info "$(now()) - ConvectionIsotopes - Tail end"
				ids = NCDataset(fncii)
			else
				ids = NCDataset(datadir("wrf","raw","2D","$(dtvec[ii]+Day(1)).nc"))
			end
			NCDatasets.load!(ids[wvar].var,aarr,:,:,1)
			close(ids)
		end


		if isaccum
			for ilat = 1 : nlat, ilon = 1 : nlon
				narr[ilon,ilat,ii] = aarr[ilon,ilat] - oarr[ilon,ilat,1]
			end
		else
			for ilat = 1 : nlat, ilon = 1 : nlon
				narr[ilon,ilat,ii] = Float32(mean(view(oarr,ilon,ilat,:)))
			end
		end

	end

	fnc = datadir("wrf","2D","$wvar-daily.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["date"]      = ndt

	nclon = defVar(ds,"longitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
		"units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
	))

	ncvar = defVar(ds,wvar,Float32,("longitude","latitude","date"),attrib=attrib)

	nclon[:] = lon
	nclat[:] = lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:] = narr

	close(ds)

end