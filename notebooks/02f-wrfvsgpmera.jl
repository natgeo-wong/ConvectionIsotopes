### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 1d7446ba-f464-44df-89e2-ae2a5726e849
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ab294fae-101f-4587-a2f4-7d72254dd421
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using NASAPrecipitation
	using ERA5Reanalysis
	using PlutoUI
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 02a. Exploring the GPM IMERG Dataset for the OTREC Domain

In this notebook, we explore the diurnal cycle of the GPM IMERG dataset for the OTREC domain using a gaussian-filtered/smoothed version of the IMERG landsea mask.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Datasets and LandSea Masks ...
"

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 1cee44fe-75b7-42a2-948e-db330cf788e8
npd = IMERGFinalHH(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ 379d74ef-e4e4-4879-8b9a-92cacb26dc3a
e5ds = ERA5Hourly(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ bceeb091-1cae-44cc-9a9a-85f3e7d0865e
evar = SingleVariable("tp")

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC",path=srcdir())

# ╔═╡ 302d3f79-9a9e-4a06-88b4-cc29fbb3e352
egeo = ERA5Region(geo)

# ╔═╡ a01e1721-23cf-4a3f-b5aa-0189b1a113a3
lsd = getLandSea(npd,geo)

# ╔═╡ 6b9b0803-a7a4-415f-ac4c-fe6e05c40c3b
elsd = getLandSea(e5ds,egeo)

# ╔═╡ 7c466ad9-7c99-4af0-bbc2-c81d8575be6e
begin
	prcp_npd = zeros(length(lsd.lon),length(lsd.lat))
	for idt = npd.start : Day(1) : npd.stop
		ds = read(npd,geo,idt)
		prcp_npd[:,:] += sum(ds["precipitation"][:,:,:] * 1800,dims=3)
		close(ds)
	end
	prcp_e5 = zeros(length(elsd.lon),length(elsd.lat))
	for idt = e5ds.start : Month(1) : e5ds.stop
		ds = read(e5ds,evar,egeo,idt)
		prcp_e5[:,:] += sum(ds[evar.ID][:,:,:],dims=3) *1000
		close(ds)
	end
	wds = NCDataset(datadir("wrf3","regridded","era-RAINNC-20190801_20201231.nc"))
	wrfp_era = wds["RAINNC"][:,:,:]
	close(wds)
	wds = NCDataset(datadir("wrf3","regridded","gpm-RAINNC-20190801_20201231.nc"))
	wrfp_gpm = wds["RAINNC"][:,:,:]
	close(wds)
	wrfp_eraμ = zeros(length(elsd.lon),length(elsd.lat))
	for ilat = 1 : length(elsd.lat), ilon = 1 : length(elsd.lon)
		wrfp_eraii = wrfp_era[ilon,ilat,:]
		wrfp_eraμ[ilon,ilat] = sum(wrfp_eraii[.!isnan.(wrfp_eraii)])
	end
	wrfp_gpmμ = zeros(length(lsd.lon),length(lsd.lat))
	for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
		wrfp_gpmii = wrfp_gpm[ilon,ilat,:]
		wrfp_gpmμ[ilon,ilat] = sum(wrfp_gpmii[.!isnan.(wrfp_gpmii)])
	end
end

# ╔═╡ 02f026eb-1a1b-47e9-be60-0b2ef5fbf781
begin
	pplt.close()
	fig,axs = pplt.subplots(nrows=2,ncols=3,axwidth=1.5,hspace=1.5)
	
	axs[1].pcolormesh(lsd.lon,lsd.lat,prcp_npd',levels=(1:15).*1000,extend="both")
	c1 = axs[2].pcolormesh(lsd.lon,lsd.lat,wrfp_gpmμ',levels=(1:15).*1000,extend="both")
	c2 = axs[3].pcolormesh(lsd.lon,lsd.lat,(wrfp_gpmμ.-prcp_npd)',levels=(-5:5).*1000,extend="both",cmap="drywet")
	axs[4].pcolormesh(elsd.lon,elsd.lat,prcp_e5',levels=(1:15).*1000,extend="both")
	axs[5].pcolormesh(elsd.lon,elsd.lat,wrfp_eraμ',levels=(1:15).*1000,extend="both")
	axs[6].pcolormesh(elsd.lon,elsd.lat,(wrfp_eraμ.-prcp_e5)',levels=(-5:5).*1000,extend="both",cmap="drywet")

	for ax in axs
		ax.plot(x,y,c="k",lw=0.5)
		ax.format(xlim=(270,285),ylim=(0,15),ylocator=0:5:15)
	end

	axs[1].format(ultitle="(a) IMERGv7",suptitle="Accumulated Rainfall (Aug 2019 - Dec 2020)")
	axs[2].format(ultitle="(b) WRF (IMERG grid)")
	axs[3].format(ultitle="(c) WRF - IMERGv7")
	axs[4].format(ultitle="(d) ERA5")
	axs[5].format(ultitle="(e) WRF (ERA5 grid)")
	axs[6].format(ultitle="(f) WRF - ERA5")

	fig.colorbar(c1,row=[1],label="mm")
	fig.colorbar(c2,row=[2],label="mm")
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 04bf9f74-475c-4140-be13-d28c9407c8cb
md"
### Station by Station
"

# ╔═╡ 91ae8e84-29c3-4011-8a2e-fd860c498e5b
sgeo = GeoRegion("OTREC_wrf_stn07",path=srcdir())

# ╔═╡ 4a52749a-cfb3-4dd6-8ae9-6d7d4baadb88
begin
	tds = NCDataset(datadir("wrf3","processed","$(sgeo.ID)-gpmrain-20190801_20201231.nc"))
	timetime  = tds["time"][:]
	sprcp_gpm = tds["precipitation"][:] * 3600
	close(tds)
	tds = NCDataset(datadir("wrf3","processed","$(sgeo.ID)-QBUDGET-20190801_20201231.nc"))
	sprcp_wrf = tds["P"][:]
	close(tds)
end

# ╔═╡ b8375412-4c53-4e84-8e40-4fdf15648357
function smooth(vec;days)

	hrs = days * 24
	ndt = length(vec)
	nvec = zeros(length(vec))

	for it = 1 : ndt-hrs
		nvec[it] = mean(vec[it.+(0:(hrs-1))])
	end

	return nvec

end

# ╔═╡ fd1eb701-9641-433f-bd13-8feda2332574
begin
	pplt.close(); f2,a2 = pplt.subplots([1,1,1,2],aspect=3)
	
	a2[1].plot(timetime,smooth(sprcp_gpm,days=30))
	a2[1].plot(timetime,smooth(sprcp_wrf,days=30))

	a2[2].scatter(smooth(sprcp_wrf,days=30),smooth(sprcp_gpm,days=30),s=2)

	a2[1].format(ylabel="Hourly Rainfall",suptitle="30-day moving average ($(sgeo.name))")
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─1cee44fe-75b7-42a2-948e-db330cf788e8
# ╟─379d74ef-e4e4-4879-8b9a-92cacb26dc3a
# ╠═bceeb091-1cae-44cc-9a9a-85f3e7d0865e
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─302d3f79-9a9e-4a06-88b4-cc29fbb3e352
# ╟─a01e1721-23cf-4a3f-b5aa-0189b1a113a3
# ╟─6b9b0803-a7a4-415f-ac4c-fe6e05c40c3b
# ╠═7c466ad9-7c99-4af0-bbc2-c81d8575be6e
# ╠═02f026eb-1a1b-47e9-be60-0b2ef5fbf781
# ╟─04bf9f74-475c-4140-be13-d28c9407c8cb
# ╠═91ae8e84-29c3-4011-8a2e-fd860c498e5b
# ╠═4a52749a-cfb3-4dd6-8ae9-6d7d4baadb88
# ╠═b8375412-4c53-4e84-8e40-4fdf15648357
# ╠═fd1eb701-9641-433f-bd13-8feda2332574
