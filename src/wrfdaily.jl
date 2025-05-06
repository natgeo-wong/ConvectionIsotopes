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

    ds  = NCDataset(datadir("wrf3","grid.nc"))
	if wvar == "U"
		lon = ds["longitude_u"][:,:]
		lat = ds["latitude_u"][:,:]
	elseif wvar == "V"
		lon = ds["longitude_v"][:,:]
		lat = ds["latitude_v"][:,:]
	else
		lon = ds["longitude"][:,:]
		lat = ds["latitude"][:,:]
	end
	nlon,nlat = size(lon)[[1,2]]
	close(ds)

	ds  = NCDataset(datadir("wrf3","raw","$start.nc"))
	nlvl = size(ds[wvar])[3]
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	dtvec = start : Day(1) : stop; ndt = length(dtvec)

	oarr = zeros(Float32,nlon,nlat,nlvl,24)
	narr = zeros(Float32,nlon,nlat,nlvl,ndt)

	for ii in 1 : ndt

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for $(dtvec[ii])"
		flush(stderr)

		fnc = datadir("wrf3","raw","$(dtvec[ii]).nc")
		if isfile(fnc)

			ids = NCDataset(fnc)
			if ids.dim["Time"] == 24
				NCDatasets.load!(ids[wvar].var,oarr,:,:,:,:)
				close(ids)

				for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
					narr[ilon,ilat,ilvl,ii] = Float32(mean(view(oarr,ilon,ilat,ilvl,:)))
				end
			else
				@warn "$(now()) - ConvectionIsotopes - Unable to extract $wvar data for $(dtvec[ii]), 24 hours not given, set to NaN"
				close(ids)
				narr[:,:,:,ii] .= NaN32
			end

		else

			@warn "$(now()) - ConvectionIsotopes - Unable to extract $wvar data for $(dtvec[ii]), the file does not exist, set to NaN"
			narr[:,:,:,ii] .= NaN32

		end

	end

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	fnc = datadir("wrf3","3D","$wvar-daily-$(dtbegstr)_$(dtbegend).nc")
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

	nclon[:,:] = lon
	nclat[:,:] = lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:,:,:,:] = narr

	close(ds)

end

function wrf2Ddaily(
    wvar :: AbstractString;
	start :: Date,
	stop  :: Date,
	isaccum   :: Bool = false,
)

	ds  = NCDataset(datadir("wrf3","grid.nc"))
	lon = ds["longitude"][:,:]
	lat = ds["latitude"][:,:]
	nlon,nlat = size(lon)[[1,2]]
	close(ds)

	ds  = NCDataset(datadir("wrf3","raw","$start.nc"))
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	dtvec = start : Day(1) : stop; ndt = length(dtvec)

	oarr = zeros(Float32,nlon,nlat,24)
	narr = zeros(Float32,nlon,nlat,ndt)
	if isaccum
		aarr = zeros(Float32,nlon,nlat)
	end

	for ii in 1 : ndt

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for $(dtvec[ii])"
		flush(stderr)

		fnc = datadir("wrf3","raw","$(dtvec[ii]).nc")
		if isfile(fnc)

			ids = NCDataset(fnc)
			if ids.dim["Time"] == 24
				NCDatasets.load!(ids[wvar].var,oarr,:,:,:)
				close(ids)
				if isaccum
					fncii = datadir("wrf3","raw","$(dtvec[ii]+Day(1))-e.nc")
					if isfile(fncii)
						@info "$(now()) - ConvectionIsotopes - Tail end"
						ids = NCDataset(fncii)
					else
						ids = NCDataset(datadir("wrf3","raw","$(dtvec[ii]+Day(1)).nc"))
					end
					NCDatasets.load!(ids[wvar].var,aarr,:,:,1)
					close(ids)
				end


				if isaccum
					for ilat = 1 : nlat, ilon = 1 : nlon
						if aarr[ilon,ilat] >= oarr[ilon,ilat,1]
							narr[ilon,ilat,ii] = aarr[ilon,ilat] - oarr[ilon,ilat,1]
						else
							narr[ilon,ilat,ii] = oarr[ilon,ilat,end] - oarr[ilon,ilat,1]
						end
					end
				else
					for ilat = 1 : nlat, ilon = 1 : nlon
						narr[ilon,ilat,ii] = Float32(mean(view(oarr,ilon,ilat,:)))
					end
				end
			else
				@warn "$(now()) - ConvectionIsotopes - Unable to extract $wvar data for $(dtvec[ii]), 24 hours not given, set to NaN"
				close(ids)
				narr[:,:,ii] .= NaN32
			end

		else

			@warn "$(now()) - ConvectionIsotopes - Unable to extract $wvar data for $(dtvec[ii]), the file does not exist"
			narr[:,:,ii] .= NaN32

		end

	end

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	fnc = datadir("wrf3","2D","$wvar-daily-$(dtbegstr)_$(dtbegend).nc")
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

	nclon[:,:] = lon
	nclat[:,:] = lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:,:,:] = narr

	close(ds)

end