using DrWatson
@quickactivate "ConvectionIsotopes"

using NCDatasets

function wrfgrid()
    
    ds   = NCDataset(datadir("wrf","raw","2D","2019-08-01.nc"))
    lon   = ds["XLONG"][:,:,1]; nlon = size(lon,1)
    lat   = ds["XLAT"][:,:,1];  nlat = size(lat,2)
    lon_u = ds["XLONG_U"][:,:,1]; lat_u = ds["XLAT_U"][:,:,1]
    lon_v = ds["XLONG_V"][:,:,1]; lat_v = ds["XLAT_V"][:,:,1]
    close(ds)

    fnc = datadir("wrf","grid.nc")
	if isfile(fnc); rm(fnc,force=true) end

	ds = NCDataset(fnc,"c")

	ds.dim["longitude"] = nlon
	ds.dim["latitude"]  = nlat
	ds.dim["longitude_u"] = nlon + 1
	ds.dim["latitude_v"]  = nlat + 1

	nclon = defVar(ds,"longitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("longitude","latitude"),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "latitude",
	))

	nclonu = defVar(ds,"longitude_u",Float32,("longitude_u","latitude"),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "U-grid longitude",
	))

	nclatu = defVar(ds,"latitude_u",Float32,("longitude_u","latitude"),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "U-grid latitude",
	))

	nclonv = defVar(ds,"longitude_v",Float32,("longitude","latitude_v"),attrib=Dict(
		"units"     => "degrees_east",
        "long_name" => "V-grid longitude",
	))

	nclatv = defVar(ds,"latitude_v",Float32,("longitude","latitude_v"),attrib=Dict(
		"units"     => "degrees_north",
        "long_name" => "V-grid latitude",
	))

	nclon[:] = lon; nclonu[:] = lon_u; nclonv[:] = lon_v
	nclat[:] = lat; nclatu[:] = lat_u; nclatv[:] = lat_v

	close(ds)

end

wrfgrid()