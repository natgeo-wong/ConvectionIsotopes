### A Pluto.jl notebook ###
# v0.20.4

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
	pplt = pyimport("ultraplot")

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

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Monthly(start=Date(1979,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 9473b30b-787a-495b-bdfe-b386c9995761
evar = SingleVariable("p_wwgt",path=srcdir())

# ╔═╡ c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
evar_pre = SingleVariable("sp")

# ╔═╡ c7d0ae59-8d78-4946-8348-9fa570197b0b
egeo = ERA5Region(GeoRegion([0,360,360,0,0],[30,30,-30,-30,30],ID="TRP"))

# ╔═╡ cb7c6118-e25b-462a-84a5-751ea0682b52
elsd = getLandSea(e5ds,egeo)

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
begin
	tmp = zeros(Float32,length(elsd.lon),length(elsd.lat),12)
	wp  = zeros(Float32,length(elsd.lon),length(elsd.lat),12)
	sp  = zeros(Float32,length(elsd.lon),length(elsd.lat),12)
	wσ  = zeros(Float32,length(elsd.lon),length(elsd.lat),12)
	cnt = zeros(Float32,length(elsd.lon),length(elsd.lat),12)
	for dt in e5ds.start : Year(1) : e5ds.stop
		wpds = read(e5ds,evar,egeo,dt,quiet=true)
		spds = read(e5ds,evar_pre,egeo,dt,quiet=true)
		tmp[:,:,:] .= nomissing(wpds[evar.ID][:,:,:],0) / 100
		wp[:,:,:]  += tmp
		cnt[:,:,:] += Int.(.!iszero.(tmp))
		sp[:,:,:]  .= nomissing(spds[evar_pre.ID][:,:,:]) / 100
		close(wpds)
		close(spds)
		wσ[:,:,:]  += tmp ./ sp
	end
	wp  = wp ./ cnt; wp_yr  = dropdims(mean(wp,dims=3),dims=3)
	wσ  = wσ ./ cnt; wσ_yr  = dropdims(mean(wσ,dims=3),dims=3)
	# wp_yr[elsd.z.>500] .= NaN
	# wσ_yr[elsd.z.>500] .= NaN
	md"Loading and filtering w-weighted pressure and sigma data and counts ..."
end

# ╔═╡ 3ee68a0b-82cc-4c38-8545-183e0940bbdc
begin
	ds = NCDataset(datadir("TRPx0.25-p_wwgt-compiled-20190801_20201231.nc"))
	wpii = ds[evar.ID][:,:] ./ 100
	close(ds)
end

# ╔═╡ 4342521a-6a0a-4d30-a5f3-975514d9006f
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=2,aspect=9,axwidth=6)
	
	c1 = axs[1].pcolormesh(elsd.lon,elsd.lat,wpii[:,:,1]',levels=(9:0.5:15).*50,extend="both",cmap="drywet_r")
	c2 = axs[2].pcolormesh(elsd.lon,elsd.lat,wσ_yr',levels=(9:0.5:15)./20,extend="both",cmap="drywet_r")

	for ax in axs
		ax.plot(x,y,c="k",lw=0.5)
		ax.format(
			ylim=(-15,15),ylabel=L"Latitude / $\degree$",ylocator=-15:15:15,
			xlim=(60,330),xlabel=L"Longitude / $\degree$",xlocator=60:30:330,
		)
	end

	fig.colorbar(c1,rows=1,locator=400:100:700,minorlocator=(5:0.5:15).*50,label=L"$p_\omega$")
	fig.colorbar(c2,rows=2,locator=0.4:0.1:0.7,minorlocator=(5:0.5:15)./20,label=L"$\sigma_\omega$")
	fig.savefig(plotsdir("03a-explorewwgtpre.png"),transparent=true,dpi=400)
	load(plotsdir("03a-explorewwgtpre.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─9473b30b-787a-495b-bdfe-b386c9995761
# ╟─c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
# ╠═c7d0ae59-8d78-4946-8348-9fa570197b0b
# ╟─cb7c6118-e25b-462a-84a5-751ea0682b52
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╠═acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
# ╠═3ee68a0b-82cc-4c38-8545-183e0940bbdc
# ╠═4342521a-6a0a-4d30-a5f3-975514d9006f
