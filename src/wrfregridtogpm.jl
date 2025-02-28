using Base.Threads
using Dates
using NASAPrecipitation
using NCDatasets
using Statistics

function wrfnewregridgpm2D(
	wvar  :: AbstractString;
	start :: Date,
	stop  :: Date,
)

	dtbegstr = Dates.format(start,dateformat"yyyymmdd")
	dtbegend = Dates.format(stop,dateformat"yyyymmdd")
	timestr = "$(dtbegstr)_$(dtbegend)"
    
    ds   = NCDataset(datadir("wrf3","grid.nc"))
    wlon = ds["longitude"][:,:]; nx,ny = size(wlon)
    wlat = ds["latitude"][:,:]
    close(ds)

	ds  = NCDataset(datadir("wrf3","raw","$start.nc"))
	attrib = Dict(ds[wvar].attrib)
	close(ds)

    dtvec = start : Day(1) : stop
    ndt  = length(dtvec)

	npd  = IMERGDummy(path=datadir())
    geo  = GeoRegion("OTREC_wrf_d02",path=srcdir())
	lsd  = getLandSea(npd,geo)
	nlon = length(lsd.lon)
	nlat = length(lsd.lat)

	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wlon[ilon,ilat].-lsd.lon))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlat[ilon,ilat].-lsd.lat))
	end

	tmp1 = zeros(Float32,nx,ny,24)
	tmp2 = zeros(Float32,nx,ny)
	tmp3 = zeros(Float32,nx,ny,24)
	ndata = zeros(nlon,nlat,24,ndt) * NaN
	for idt = 1 : ndt

		@info "$(now()) - ConvectionIsotopes - Regridding $wvar data for Day $idt of $ndt"
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

				NCDatasets.load!(ds1[wvar].var,tmp1,:,:,:)
                NCDatasets.load!(ds2[wvar].var,tmp2,:,:,1)
                tmp3 .= cat(tmp1[:,:,2:end],tmp2,dims=3) .- tmp1

				Threads.@threads for idx in 1:(24 * nlat * nlon)
					# Compute indices
					ihr = div(idx - 1, (nlat * nlon)) + 1
					ilat = div(mod(idx - 1, (nlat * nlon)), nlon) + 1
					ilon = mod(idx - 1, nlon) + 1

					ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
					idata = @view tmp3[:,:,ihr]
					idata = @view idata[ind]
					ndata[ilon,ilat,ihr,idt] = mean(idata[.!isnan.(idata)])
				end
				close(ds2)

			end

			close(ds1)

		end

	end

	if !isdir(datadir("wrf3","regridded")); mkpath(datadir("wrf3","regridded")) end
	fnc = datadir("wrf3","regridded","gpm-$(wvar)-$timestr.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["date"]      = ndt * 24

	nclon = defVar(ds,"longitude",Float32,("longitude",),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("latitude",),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nctime = defVar(ds,"time",Float64,("date",),attrib=Dict(
		"units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
	))

	ncvar = defVar(ds,wvar,Float32,("longitude","latitude","date"),attrib=attrib)

	nclon[:] = lsd.lon
	nclat[:] = lsd.lat
	nctime.var[:] = collect(0 : (ndt*24 -1)) .+ 0.5
	ncvar[:,:,:] = ndata

	close(ds)

end