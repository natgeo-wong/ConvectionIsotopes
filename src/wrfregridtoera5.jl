using Dates
using ERA5Reanalysis
using NCDatasets
using Statistics

function wrfnewregridera5(
	wvar  :: AbstractString;
	start :: Date,
	stop  :: Date,
	days  :: Int = 1
)

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

	if iszero(days)
        wds = NCDataset(datadir("wrf3","2D","$wvar-daily-$timestr.nc"))
    else
        wds = NCDataset(datadir("wrf3","2D","$wvar-daily-$timestr-$smthstr.nc"))
    end
	attrib = Dict(wds[wvar].attrib)
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    wlon = ds["longitude"][:,:]; nx,ny = size(wlon)
    wlat = ds["latitude"][:,:]
    close(ds)

    dtvec = start : Day(1) : stop
    ndt  = length(dtvec)

	e5ds = ERA5Dummy(path=datadir())
    geo  = ERA5Region(GeoRegion("OTREC",path=srcdir()))
	lsd  = getLandSea(e5ds,geo)
	nlon = length(lsd.lon)
	nlat = length(lsd.lat)

	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wlon[ilon,ilat].-lsd.lon.+360))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlat[ilon,ilat].-lsd.lat))
	end

	tdata = zeros(Float32,nx,ny)
	ndata = zeros(nlon,nlat,ndt)
	for idt = 1 : ndt
		@info "$(now()) - ConvectionIsotopes - Regridding $wvar data for Day $idt of $ndt"
		flush(stderr)
		for ilat = 1 : nlat, ilon = 1 : nlon
			ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
			NCDatasets.load!(wds[wvar].var,tdata,:,:,idt)
			idata = @view tdata[ind]; ndata[ilon,ilat,idt] = mean(idata[.!isnan.(idata)])
		end
	end

	close(wds)

	if !isdir(datadir("wrf3","regridded")); mkpath(datadir("wrf3","regridded")) end
	fnc = datadir("wrf3","regridded","$wvar-$timestr-$smthstr.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["date"]      = ndt

	nclon = defVar(ds,"longitude",Float32,("longitude",),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("latitude",),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
		"units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
	))

	ncvar = defVar(ds,wvar,Float32,("longitude","latitude","date"),attrib=attrib)

	nclon[:] = lsd.lon
	nclat[:] = lsd.lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:,:,:] = ndata

	close(ds)

end

function wrfoldregridera5(
	wvar  :: AbstractString;
	days  :: Int = 1
)

    smthstr = "smooth_$(@sprintf("%02d",days))days"

	if iszero(days)
        wds = NCDataset(datadir("wrf","2D","$wvar-daily.nc"))
    else
        wds = NCDataset(datadir("wrf","2D","$wvar-daily-$smthstr.nc"))
    end
	attrib = Dict(wds[wvar].attrib)
	ndt  = wds.dim["time"]; start = wds["time"][1]
    
    ds   = NCDataset(datadir("wrf","grid.nc"))
    wlon = ds["longitude"][:,:]; nx,ny = size(wlon)
    wlat = ds["latitude"][:,:]
    close(ds)

	e5ds = ERA5Dummy(path=datadir())
    geo  = ERA5Region(GeoRegion("OTREC",path=srcdir()))
	lsd  = getLandSea(e5ds,geo)
	nlon = length(lsd.lon)
	nlat = length(lsd.lat)

	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wlon[ilon,ilat].-lsd.lon.+360))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlat[ilon,ilat].-lsd.lat))
	end

	tdata = zeros(Float32,nx,ny)
	ndata = zeros(nlon,nlat,ndt)
	for idt = 1 : ndt
		@info "$(now()) - ConvectionIsotopes - Regridding $wvar data for Day $idt of $ndt"
		flush(stderr)
		for ilat = 1 : nlat, ilon = 1 : nlon
			ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
			NCDatasets.load!(wds[wvar].var,tdata,:,:,idt)
			idata = @view tdata[ind]; ndata[ilon,ilat,idt] = mean(idata[.!isnan.(idata)])
		end
	end

	close(wds)

	if !isdir(datadir("wrf","regridded")); mkpath(datadir("wrf","regridded")) end
	fnc = datadir("wrf","regridded","$wvar-$smthstr.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["date"]      = ndt

	nclon = defVar(ds,"longitude",Float32,("longitude",),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("latitude",),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nctime = defVar(ds,"time",Int32,("date",),attrib=Dict(
		"units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
	))

	ncvar = defVar(ds,wvar,Float32,("longitude","latitude","date"),attrib=attrib)

	nclon[:] = lsd.lon
	nclat[:] = lsd.lat
	nctime.var[:] = collect(0 : (ndt-1))
	ncvar[:,:,:] = ndata

	close(ds)

end