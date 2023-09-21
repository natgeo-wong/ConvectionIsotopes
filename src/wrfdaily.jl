using Statistics

function wrf3Ddaily(
    wvar :: AbstracString;
	start :: DateTime,
	stop  :: DateTime
)

    ds  = NCDataset(datadir("wrf","raw","3D","$(start).nc"))
	nlon,nlat,nlvl = size(ds[wvar])[[1,2,3]]
	lon = ds["XLONG"][:,:,1]
	lat = ds["XLAT"][:,:,1]
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	dtvec = start : Day(1) : stop; ndt = length(dtvec)

	oarr = zeros(Float32,nlon,nlat,nlvl,8)
	narr = zeros(Float32,nlon,nlat,nlvl,ndt)

	for ii in 1 : ndt

		ids = NCDataset(datadir(fol,"$(dtvec[ii]).nc"))
		NCDatasets.load!(oarr,ids[wvar].var,:,:,:,:)
		close(ids)

		for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
			narr[ilon,ilat,ilvl,ii] = Float32(mean(view(oarr[ilon,ilat,ilvl,:])))
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
    wvar :: AbstracString;
	start :: DateTime,
	stop  :: DateTime
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

	for ii in 1 : ndt

		ids = NCDataset(datadir(fol,"$(dtvec[ii]).nc"))
		NCDatasets.load!(oarr,ids[wvar].var,:,:,:)
		close(ids)

		for ilat = 1 : nlat, ilon = 1 : nlon
			narr[ilon,ilat,ii] = Float32(mean(view(oarr[ilon,ilat,:])))
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