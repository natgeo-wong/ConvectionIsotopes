### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ e32a00ee-5f32-47a1-a983-91fb77bc5d18
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using ERA5Reanalysis
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 03a. Exploring W-Weighted Mean Pressure
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 75925166-d1e8-450e-bae6-29affd25d635
md"
### A. Defining Datasets, Variables and Regions
"

# ╔═╡ d3f4a6a9-0915-470b-83af-e05e00dff5d2
npd = IMERGMonthly(start=Date(2013),stop=Date(2020),path=datadir())

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 9473b30b-787a-495b-bdfe-b386c9995761
evar = SingleVariable("p_wwgt")

# ╔═╡ 444e1d6a-6edb-4511-a9aa-e251c5b8013b
pvar = SingleVariable("tp")

# ╔═╡ c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
evar_pre = SingleVariable("sp")

# ╔═╡ c7d0ae59-8d78-4946-8348-9fa570197b0b
geo = GeoRegion("OTREC")

# ╔═╡ 3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
ereg = ERA5Region(geo)

# ╔═╡ 692cf678-4e91-4b18-b8e9-9d462cceef04
nlsd = getLandSea(npd,geo,smooth=true,σlon=5,σlat=5)

# ╔═╡ cb7c6118-e25b-462a-84a5-751ea0682b52
elsd = getLandSea(e5ds,ereg,smooth=true,σlon=2,σlat=2)

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ 274db264-f4d0-4989-b547-c50275b46745
ls_threshold = 0.95

# ╔═╡ 59c930cd-5b7f-4047-8660-615148d1bd9f
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
begin
	tmp = zeros(length(elsd.lon),length(elsd.lat),12)
	wp  = zeros(length(elsd.lon),length(elsd.lat),12)
	sp  = zeros(length(elsd.lon),length(elsd.lat),12)
	wσ  = zeros(length(elsd.lon),length(elsd.lat),12)
	cnt = zeros(length(elsd.lon),length(elsd.lat),12)
	for dt in Date(2013) : Year(1) : Date(2020)
		wpds = read(e5ds,evar,ereg,dt,quiet=true)
		spds = read(e5ds,evar_pre,ereg,dt,quiet=true)
		tmp[:,:,:] .= nomissing(wpds[evar.ID][:],0) / 100
		wp[:,:,:]  += tmp
		cnt[:,:,:] += Int.(.!iszero.(tmp))
		sp[:,:,:]  .= nomissing(spds[evar_pre.ID][:]) / 100
		close(wpds)
		close(spds)
		wσ[:,:,:]  += tmp ./ sp
	end
	wp  = wp ./ cnt
	wσ  = wσ ./ cnt
	for ilat = 1 : length(elsd.lat), ilon = 1 : length(elsd.lon)
		if elsd.lsm[ilon,ilat] > ls_threshold
			wp[ilon,ilat,:] .= NaN
			wσ[ilon,ilat,:] .= NaN
		end
	end
	md"Loading and filtering w-weighted pressure and sigma data and counts ..."
end

# ╔═╡ 5fc05bef-09ab-4b71-b3b8-52824a4132e4
begin
	prc = zeros(length(nlsd.lon),length(nlsd.lat),12)
	for dt in Date(2013) : Year(1) : Date(2020)
		ids = read(npd,geo,dt,quiet=true)
		prc[:,:,:] += ids["precipitation"][:] / 8 * 86400
		close(ids)
	end
	md"Loading precipitation data ..."
end

# ╔═╡ 14663694-56f8-4651-ab52-8560374be699
begin
	pplt.close(); asp = (ereg.geo.E-ereg.geo.W)/(ereg.geo.N-ereg.geo.S)
	f1,a1 = pplt.subplots(ncols=4,nrows=2,axwidth=1.5,aspect=asp)
	
	c = a1[1].pcolormesh(
		elsd.lon,elsd.lat,wσ[:,:,1]',levels=0.25:0.05:0.85,
		extend="both",cmap="delta"
	)
	a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].format(ultitle="(a) Jan")

	a1[2].pcolormesh(
		elsd.lon,elsd.lat,wσ[:,:,4]',levels=0.25:0.05:0.85,
		extend="both",cmap="delta"
	)
	a1[2].plot(x,y,c="k",lw=0.5)
	a1[2].format(ultitle="(b) Apr")

	a1[3].pcolormesh(
		elsd.lon,elsd.lat,wσ[:,:,7]',levels=0.25:0.05:0.85,
		extend="both",cmap="delta"
	)
	a1[3].plot(x,y,c="k",lw=0.5)
	a1[3].format(ultitle="(c) Jul")

	a1[4].pcolormesh(
		elsd.lon,elsd.lat,wσ[:,:,10]',levels=0.25:0.05:0.85,
		extend="both",cmap="delta"
	)
	a1[4].plot(x,y,c="k",lw=0.5)
	a1[4].format(ultitle="(d) Oct")
	a1[4].colorbar(c,label=L"$\sigma_w$ / hPa")
	
	c = a1[5].pcolormesh(
		elsd.lon,elsd.lat,wp[:,:,1]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	a1[5].plot(x,y,c="k",lw=0.5)
	a1[5].format(ultitle="(a) Jan")

	a1[6].pcolormesh(
		elsd.lon,elsd.lat,wp[:,:,4]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	a1[6].plot(x,y,c="k",lw=0.5)
	a1[6].format(ultitle="(b) Apr")

	a1[7].pcolormesh(
		elsd.lon,elsd.lat,wp[:,:,7]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	a1[7].plot(x,y,c="k",lw=0.5)
	a1[7].format(ultitle="(c) Jul")

	a1[8].pcolormesh(
		elsd.lon,elsd.lat,wp[:,:,10]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	a1[8].plot(x,y,c="k",lw=0.5)
	a1[8].format(ultitle="(d) Oct")
	a1[8].colorbar(c,label=L"$p_w$ / hPa")

	for ax in a1
		ax.scatter(infody[:,2],infody[:,3],s=10,c="r")
		ax.format(
			xlim=(minimum(elsd.lon)-1,maximum(elsd.lon)+1),xlocator=270:5:285,
			ylim=(minimum(elsd.lat)-1,maximum(elsd.lat)+1),
			xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
			suptitle="W-weighted Mean Pressure (2013-2020)"
		)
	end

	f1.savefig(plotsdir("03a-wwgt_pre-$(ereg.ID).png"),transparent=false,dpi=400)
	load(plotsdir("03a-wwgt_pre-$(ereg.ID).png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─d3f4a6a9-0915-470b-83af-e05e00dff5d2
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─9473b30b-787a-495b-bdfe-b386c9995761
# ╟─444e1d6a-6edb-4511-a9aa-e251c5b8013b
# ╟─c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
# ╟─c7d0ae59-8d78-4946-8348-9fa570197b0b
# ╟─3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
# ╟─692cf678-4e91-4b18-b8e9-9d462cceef04
# ╟─cb7c6118-e25b-462a-84a5-751ea0682b52
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╠═274db264-f4d0-4989-b547-c50275b46745
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
# ╟─5fc05bef-09ab-4b71-b3b8-52824a4132e4
# ╟─14663694-56f8-4651-ab52-8560374be699
