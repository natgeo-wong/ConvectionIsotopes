### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 5e0b13eb-573e-4460-a4bf-fcb0fdf44838
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ d44838c8-72cc-47ed-ba5b-1c11b4064472
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ 88457e84-fb15-11ec-1f04-4148aa4e224a
md"
# 04d. Exploring WRF Weighted Pressure
"

# ╔═╡ 38760419-698c-437d-936b-8475359dc245
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ cfd0b1c6-d821-450b-8960-19a5ad82fe1f
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ c360328f-58df-44f1-b81a-54034a4370c5
geo = GeoRegion("OTREC")

# ╔═╡ 5bd26a8e-de5e-4693-a9a0-8fd229452ba8
egeo = ERA5Region(geo)

# ╔═╡ a0b66946-d425-4b96-902d-9b416d6d3c21
lsd = getLandSea(e5ds,egeo)

# ╔═╡ ff4e498d-1884-4c80-9fcd-2dd956f55f86
begin
	wds = NCDataset(datadir("wrf","2D","AUG-p_wwgt"))
	wln = wds["longitude"][:] .+ 360; nwlon = size(wln,1)
	wlt = wds["latitude"][:]; 		  nwlat = size(wlt,2)
	close(wds)
	md"Finding Grid Size and Points for the WRF simulations ..."
end

# ╔═╡ e45ca11d-270b-4aa5-a988-a96e45fdde6f
begin
	ipnt_lon = zeros(Int,nwlon,nwlat)
	ipnt_lat = zeros(Int,nwlon,nwlat)
	ind      = zeros(Bool,nwlon,nwlat)
	for ilat = 1 : nwlat, ilon = 1 : nwlon
	
		ipnt_lon[ilon,ilat] = argmin(abs.(wln[ilon,ilat].-lsd.lon))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlt[ilon,ilat].-lsd.lat))
	
	end
	md"Finding closest ERA5 points to each of the WRF points ..."
end

# ╔═╡ 8502570a-5634-4a6c-8fc8-53b850b95da4
begin
	wprcp = zeros(length(lsd.lon),length(lsd.lat),5)
	wwgtp = zeros(length(lsd.lon),length(lsd.lat),5)
	md"Preallocation of arrays for precipitation and weighted pressure ..."
end

# ╔═╡ b19951bf-5317-4dbe-91d4-b0c724a608fe
for imo = 1 : 5
	mostr = uppercase(monthabbr(7+imo))
	ds  = NCDataset(datadir("wrf","2D","$mostr-p_wwgt"))
	ipre = ds["p_wwgt"][:] / 100
	close(ds)
	ds  = NCDataset(datadir("wrf","2D","$mostr-RAINNC"))
	iprc = ds["RAINNC"][:]
	close(ds)
	for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
		ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
		iipre = @view ipre[ind]
		wwgtp[ilon,ilat,imo] = mean(iipre[.!isnan.(iipre)])
		wprcp[ilon,ilat,imo] = mean(@view iprc[ind])
	end
	for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
		if wprcp[ilon,ilat,imo] < 5
			wwgtp[ilon,ilat,imo] = NaN
		end
	end
end

# ╔═╡ 6abac2b3-f712-4c56-a946-a206f01d00f8
e5mo = ERA5Monthly(start=Date(2019),stop=Date(2019),path=datadir())

# ╔═╡ d592e8bc-9471-49c5-9308-ffd26b70d1a9
evar = SingleVariable("p_wwgt")

# ╔═╡ e2b1e538-c6ff-470d-95f9-73e48fe44cdb
etpr = SingleVariable("tp")

# ╔═╡ 1c852ca2-9726-4674-a292-33d9e781c73d
begin
	eds = read(e5mo,evar,egeo,Date(2019))
	ewp = nomissing(eds[evar.ID][:,:,8:12],NaN) ./ 100
	close(eds)
	eds = read(e5mo,etpr,egeo,Date(2019))
	etp = nomissing(eds[etpr.ID][:,:,8:12],NaN) * 1000
	close(eds)
end

# ╔═╡ 21c0da18-5d64-42ed-850f-072082b0df9a
npd = IMERGMonthly(start=Date(2019),stop=Date(2019),path=datadir())

# ╔═╡ d9e3e054-30eb-48a5-b69b-53b5c0f35a89
npd_lsd = getLandSea(npd,geo)

# ╔═╡ 3e2cfbf4-d52c-48e8-b36d-89a2194c8177
begin
	npdds = read(npd,geo,npd.start)
	imerg = npdds["precipitation"][:] * 86400
	close(npdds)
end

# ╔═╡ c0ea0006-6b80-49b8-8497-41692097e4ae
begin
	pplt.close(); fig,axs = pplt.subplots(
		nrows=3,ncols=5,axwidth=1.25
	)
	
	c = axs[11].pcolormesh(
		lsd.lon,lsd.lat,wprcp[:,:,1]',levels=5:5:45,cmap="blues",extend="both"
	)

	axs[12].pcolormesh(
		lsd.lon,lsd.lat,wprcp[:,:,2]',
		levels=5:5:45,cmap="blues",extend="both",
	)

	axs[13].pcolormesh(
		lsd.lon,lsd.lat,wprcp[:,:,3]',levels=5:5:45,cmap="blues",extend="both",
	)

	axs[14].pcolormesh(
		lsd.lon,lsd.lat,wprcp[:,:,4]',levels=5:5:45,cmap="blues",extend="both",
	)

	axs[15].pcolormesh(
		lsd.lon,lsd.lat,wprcp[:,:,5]',levels=5:5:45,cmap="blues",extend="both",
	)

	axs[6].pcolormesh(
		lsd.lon,lsd.lat,etp[:,:,1]',levels=5:5:45,cmap="blues",extend="both",
	)
	
	axs[7].pcolormesh(
		lsd.lon,lsd.lat,etp[:,:,2]',levels=5:5:45,cmap="blues",extend="both",
	)
	
	axs[8].pcolormesh(
		lsd.lon,lsd.lat,etp[:,:,3]',levels=5:5:45,cmap="blues",extend="both",
	)
	
	axs[9].pcolormesh(
		lsd.lon,lsd.lat,etp[:,:,4]',levels=5:5:45,cmap="blues",extend="both"
	)
	
	axs[10].pcolormesh(
		lsd.lon,lsd.lat,etp[:,:,5]',levels=5:5:45,cmap="blues",extend="both"
	)
	
	axs[1].pcolormesh(
		npd_lsd.lon.+360,npd_lsd.lat,imerg[:,:,8]',
		levels=5:5:45,cmap="blues",extend="both",
	)
	
	axs[2].pcolormesh(
		npd_lsd.lon.+360,npd_lsd.lat,imerg[:,:,9]',
		levels=5:5:45,cmap="blues",extend="both"
	)
	
	axs[3].pcolormesh(
		npd_lsd.lon.+360,npd_lsd.lat,imerg[:,:,10]',
		levels=5:5:45,cmap="blues",extend="both"
	)
	
	axs[4].pcolormesh(
		npd_lsd.lon.+360,npd_lsd.lat,imerg[:,:,11]',
		levels=5:5:45,cmap="blues",extend="both"
	)

	axs[5].pcolormesh(
		npd_lsd.lon.+360,npd_lsd.lat,imerg[:,:,12]',
		levels=5:5:45,cmap="blues",extend="both",
	)

	for ax in axs
		ax.plot(x,y,lw=1,c="k")
		ax.format(
			xlocator=270:5:285,ylocator=0:5:15,
			xlim=(270,285),ylim=(0,15),
			ylabel=L"Latitude / $\degree$",
			xlabel=L"Longitude / $\degree$",
		)
	end

	axs[1].format(ltitle="(a) AUG 2019")
	axs[2].format(ltitle="(b) SEP 2019")
	axs[3].format(ltitle="(c) OCT 2019")
	axs[4].format(ltitle="(d) NOV 2019")
	axs[5].format(ltitle="(e) DEC 2019",urtitle="IMERG")
	axs[10].format(urtitle="ERA5")
	axs[15].format(urtitle="WRF")

	fig.colorbar(c,length=0.6,label=L"Rainfall Rate / mm day$^{-}$")
	fig.savefig(plotsdir("05a-raindatasets.png"),transparent=false,dpi=400)
	load(plotsdir("05a-raindatasets.png"))
end

# ╔═╡ Cell order:
# ╟─88457e84-fb15-11ec-1f04-4148aa4e224a
# ╟─5e0b13eb-573e-4460-a4bf-fcb0fdf44838
# ╟─d44838c8-72cc-47ed-ba5b-1c11b4064472
# ╟─38760419-698c-437d-936b-8475359dc245
# ╟─cfd0b1c6-d821-450b-8960-19a5ad82fe1f
# ╟─c360328f-58df-44f1-b81a-54034a4370c5
# ╟─5bd26a8e-de5e-4693-a9a0-8fd229452ba8
# ╟─a0b66946-d425-4b96-902d-9b416d6d3c21
# ╟─ff4e498d-1884-4c80-9fcd-2dd956f55f86
# ╟─e45ca11d-270b-4aa5-a988-a96e45fdde6f
# ╟─8502570a-5634-4a6c-8fc8-53b850b95da4
# ╟─b19951bf-5317-4dbe-91d4-b0c724a608fe
# ╟─6abac2b3-f712-4c56-a946-a206f01d00f8
# ╟─d592e8bc-9471-49c5-9308-ffd26b70d1a9
# ╟─e2b1e538-c6ff-470d-95f9-73e48fe44cdb
# ╟─1c852ca2-9726-4674-a292-33d9e781c73d
# ╟─21c0da18-5d64-42ed-850f-072082b0df9a
# ╟─d9e3e054-30eb-48a5-b69b-53b5c0f35a89
# ╟─3e2cfbf4-d52c-48e8-b36d-89a2194c8177
# ╟─c0ea0006-6b80-49b8-8497-41692097e4ae
