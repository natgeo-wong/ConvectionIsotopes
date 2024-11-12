### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 04c. Column-mean $\sigma$ vs Precipitation

Text
"

# ╔═╡ 5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 30d2be4a-9bdd-4939-8785-413cd7a59f78
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 5981e960-e7b4-40a9-9bfe-2ef90220e225
e5ds = ERA5Monthly(
	dtbeg=Date(2013,1,1),dtend=Date(2013,12,31),
	eroot="/n/kuangdss01/lab"
)

# ╔═╡ a9b15e68-ae2d-4698-a650-2f5f24c8888a
evar_wp = SingleVariable("p_wwgt")

# ╔═╡ 0d6e5046-273d-4e8b-aa4e-8ffa8a46fc9e
evar_sp = SingleVariable("sp")

# ╔═╡ b0056575-6850-40e3-88d7-175b11443bc1
evar_tp = SingleVariable("tp")

# ╔═╡ 43cd84a6-53b8-4e5c-b327-3d0c30f90e5c
egeo = GeoRegion("OTREC")

# ╔═╡ 228728ed-bfab-405e-83f8-8236d6952c90
ereg = ERA5Region(egeo)

# ╔═╡ f0751d18-0d6b-4d83-a79b-21d463972edf
md"
### B. Loading a Sample Dataset
"

# ╔═╡ 41157e32-a476-44aa-92a5-b52242c7d3de
begin
	_,lon,lat = read(e5ds,evar_wp,ereg,Date(2014,1,1),lonlat=true)
	nlon = length(lon); nlat = length(lat)
	σ_sfc = zeros(Float32,nlon,nlat,12)
	prcp  = zeros(Float32,nlon,nlat,12)
	for yr in 2013 : 2021
		ds_wp = read(e5ds,evar_wp,ereg,Date(2014,1,1))
		ds_sp = read(e5ds,evar_sp,ereg,Date(2014,1,1))
		ds_tp = read(e5ds,evar_tp,ereg,Date(2014,1,1))
		mo = 9
		p_wwgt = nomissing(ds_wp[evar_wp.varID][:],NaN)
		p_sfc  = nomissing(ds_sp[evar_sp.varID][:],NaN)
		σ_sfc[:,:,:] += p_wwgt ./ p_sfc
		prcp[:,:,:]  += nomissing(ds_tp[evar_tp.varID][:],NaN)
		close(ds_wp)
		close(ds_sp)
		close(ds_tp)
	end
	σ_sfc  = σ_sfc  / 9
	prcp   = prcp   / 9
	md"We load the surface pressure and weighted column-mean pressure"
end

# ╔═╡ 9797dd01-4d3d-48c7-bce4-6b75e1786188
mo = 12

# ╔═╡ c43f9189-ef45-46f9-ab82-10d5e6e2298b
begin
	pbin = 0:0.4:40; np = length(pbin)
	σbin = 0.3:0.006:0.9; nσ = length(σbin)
	pvσm = fit(
		Histogram,
		(prcp[:,:,mo][:]*1000,σ_sfc[:,:,mo][:]),
		(pbin,σbin)
	).weights
	md"Binning precipitation and weighted $\sigma$ data for $(monthname(mo))"
end

# ╔═╡ b401447c-acf6-4eb3-b50b-24e5afaa850d
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=2,sharey=0)

	c = axs[1].contourf(
		lon,lat,prcp[:,:,mo]'*1000,
		levels=(1:0.25:5)*5,cmap="Blues",extend="both"
	)
	axs[1].plot(x,y,c="k")
	axs[1].colorbar(c,loc="r",label="Precipitation / mm")
	
	c = axs[2].contourf(
		lon,lat,σ_sfc[:,:,mo]',
		levels=(3:0.5:9) ./ 10,cmap="drywet_r",extend="both"
	)
	axs[2].plot(x,y,c="k")
	axs[2].colorbar(c,loc="r",label=L"$\sigma_{p_w}$")
	
	c = axs[3].pcolormesh(
		pbin,σbin,pvσm' / sum(pvσm) * np * nσ,
		levels=[0.1,0.2,0.5,1,2,5,10,20,50,100],extend="both"
	)
	axs[3].colorbar(c,loc="r",label="Normalized Density")

	axs[1].format(xlim=(260,305),ylim=(-15,30),suptitle="$(monthname(mo))")
	axs[2].format(xlim=(260,305),ylim=(-15,30))
	axs[3].format(
		xlim=(minimum(pbin),maximum(pbin)),
		ylim=(minimum(σbin),maximum(σbin)),
		grid=true
	)
	
	fig.savefig(plotsdir("04c-p_wwgtVprcp-$(mo).png"),dpi=150,transparent=false)
	load(plotsdir("04c-p_wwgtVprcp-$(mo).png"))
end

# ╔═╡ e84c3b61-e52b-4c52-b7ca-92f281436237
md"Because many of our station datasets occur inland where the topography is non-negligible, as opposed to Torri et al. (2015) where they use the raw values of weighted column-mean pressure, I propose to normalize these values by the surface pressure.  This would allow us to do generalize our comparisons over land and ocean."

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─30d2be4a-9bdd-4939-8785-413cd7a59f78
# ╟─5981e960-e7b4-40a9-9bfe-2ef90220e225
# ╟─a9b15e68-ae2d-4698-a650-2f5f24c8888a
# ╟─0d6e5046-273d-4e8b-aa4e-8ffa8a46fc9e
# ╟─b0056575-6850-40e3-88d7-175b11443bc1
# ╠═43cd84a6-53b8-4e5c-b327-3d0c30f90e5c
# ╟─228728ed-bfab-405e-83f8-8236d6952c90
# ╟─f0751d18-0d6b-4d83-a79b-21d463972edf
# ╟─41157e32-a476-44aa-92a5-b52242c7d3de
# ╠═9797dd01-4d3d-48c7-bce4-6b75e1786188
# ╟─c43f9189-ef45-46f9-ab82-10d5e6e2298b
# ╟─b401447c-acf6-4eb3-b50b-24e5afaa850d
# ╟─e84c3b61-e52b-4c52-b7ca-92f281436237
