using DrWatson
using Dates
using Logging
using NCDatasets
using Printf
using Statistics

function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

function wrf3Dsmooth(
    wvar  :: AbstractString;
	start :: Date,
	stop  :: Date,
	days  :: Int = 1
)

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	fnc = datadir("wrf3","3D","$wvar-daily-$(dtbegstr)_$(dtbegend).nc")
	
    ds   = NCDataset(fnc)
	lon  = ds["longitude"][:,:]; nlon = ds.dim["longitude"]
    lat  = ds["latitude"][:,:];  nlat = ds.dim["latitude"]
	nlvl = ds.dim["levels"]
    ndt  = ds.dim["date"]; start = ds["time"][1]
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	buffer,weights = calculatebufferweights(days)

    oarr = zeros(Float32,nlon,nlat,nlvl,1+buffer*2)
	narr = fill(NaN32,(nlon,nlat,nlvl,ndt))
	smth = zeros(Float32,1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for Day $ii of $ndt"
		flush(stderr)

		ids = NCDataset(datadir("wrf3","3D","$(wvar)-daily.nc"))
		NCDatasets.load!(ids[wvar].var,oarr,:,:,:,ii.+(-buffer:buffer))
		close(ids)

		for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon
			for ismth = 1 .+ (0 : buffer*2)
				smth[ismth] = oarr[ilon,ilat,ilvl,ismth] * weights[ismth]
			end
			narr[ilon,ilat,ilvl,ii] = Float32(sum(smth))
		end

	end

	fnc = datadir("wrf3","3D","$wvar-daily-$(dtbegstr)_$(dtbegend)-smooth_$(@sprintf("%02d",days))days.nc")
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

function wrf2Dsmooth(
    wvar :: AbstractString;
	start :: Date,
	stop  :: Date,
	days :: Int = 1
)

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	fnc = datadir("wrf3","2D","$wvar-daily-$(dtbegstr)_$(dtbegend).nc")

    ds   = NCDataset(fnc)
	lon  = ds["longitude"][:]; nlon = ds.dim["longitude"]
    lat  = ds["latitude"][:];  nlat = ds.dim["latitude"]
    ndt  = ds.dim["date"]; start = ds["time"][1]
	attrib = Dict(ds[wvar].attrib)
	close(ds)

	buffer,weights = calculatebufferweights(days)

    oarr = zeros(Float32,nlon,nlat,1+buffer*2)
	narr = fill(NaN32,(nlon,nlat,ndt))
	smth = zeros(Float32,1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		@info "$(now()) - ConvectionIsotopes - Extracting $wvar data for Day $ii of $ndt"
		flush(stderr)

		ids = NCDataset(datadir("wrf3","2D","$(wvar)-daily.nc"))
		NCDatasets.load!(ids[wvar].var,oarr,:,:,ii.+(-buffer:buffer))
		close(ids)

		for ilat = 1 : nlat, ilon = 1 : nlon
			for ismth = 1 .+ (0 : buffer*2)
				smth[ismth] = oarr[ilon,ilat,ismth] * weights[ismth]
			end
			narr[ilon,ilat,ii] = Float32(sum(smth))
		end

	end

	fnc = datadir("wrf3","2D","$wvar-daily-$(dtbegstr)_$(dtbegend)-smooth_$(@sprintf("%02d",days))days.nc")
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