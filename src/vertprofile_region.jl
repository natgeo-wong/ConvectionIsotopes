using DrWatson
using ERA5Reanalysis
using Statistics

function vertprofile(
    e5ds :: ERA5Monthly,
    ereg :: ERA5Region,
)

    lsd = ERA5Reanalysis.getLandSea(ereg,eroot=datadir("emask"))
    nlon = length(lsd.lon)
	nlat = length(lsd.lat)
	
	var_wp = SingleVariable("p_wwgt")
	var_w  = Vector{PressureVariable}(undef,32)

	wp  = zeros(nlon,nlat)
	w   = zeros(nlon,nlat,32)
	cnt = 0
	
	pp  = era5Pressures(); pp = pp[pp.>=10]
	for ip = 1 : 32
		var_w[ip] = PressureVariable("w",hPa=pp[ip])
	end

	for dt = Date(2013) : Year(1) : Date(2021)

		ids = read(e5ds,var_wp,ereg,dt)
		wp .= nomissing(ids[var_wp.varID][:,:,1],NaN)
		close(ids)

		for ip = 1 : 32

			ids = read(e5ds,var_w[ip],ereg,dt)
			w[:,:,ip] .= ids[var_w[ip].varID][:,:,1]
			close(ids)

		end

		for ilon = 1 : nlon, ilat = 1 : nlat

			if !isnan(wp[ilon,ilat])
				cnt += 1
			end
			
		end
		
	end

    wprf_tmp = zeros(32,cnt)

    cnt = 0
    for dt = Date(2013) : Year(1) : Date(2021)

		ids = read(e5ds,var_wp,ereg,dt)
		wp .= nomissing(ids[var_wp.varID][:,:,1],NaN)
		close(ids)

		for ip = 1 : 32

			ids = read(e5ds,var_w[ip],ereg,dt)
			w[:,:,ip] .= ids[var_w[ip].varID][:,:,1]
			close(ids)

		end

		for ilon = 1 : nlon, ilat = 1 : nlat

			if !isnan(wp[ilon,ilat])
				cnt += 1
				wprf_tmp[:,cnt] .= w[ilon,ilat,:]
			end
			
		end
		
	end

    wprf_μ = dropdims(mean(wprf_tmp,dims=2),dims=2)
    wprf_σ = dropdims(std( wprf_tmp,dims=2),dims=2)

	return wprf_μ, wprf_σ

end