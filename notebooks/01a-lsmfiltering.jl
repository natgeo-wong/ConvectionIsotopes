### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a675084e-f638-11eb-2930-b70006b9445c
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 4e08778c-cd46-4a79-9b6a-96ca6447e975
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	using GeoRegions
	using ImageFiltering
	using NCDatasets

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
md"
# 01a. Land-Sea Mask Filtering

In this notebook, we do a gaussian filtering for the land-sea mask in order to determine and better filter out the coastline and near-coastline areas in the region of analysis (the OTREC GeoRegion).  This will help us define further subGeoRegions, namely the Atlantic and Pacific domains, without too much trouble, for future notebooks
"

# ╔═╡ 75ea588a-64df-4b31-9dae-f99f425ee55a
md"
### A. Loading Land-Sea Mask Data
"

# ╔═╡ 0d950c5d-8c29-4eed-9716-e7b979c59ff6
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ 79e35867-824c-4b3a-8c7b-3abb8f9d9761
egeo = ERA5Region("OTREC",resolution=0.25)

# ╔═╡ f3773f6f-55cc-4546-92c7-dd3d8955edab
lsd_raw = getLandSea(e5ds,egeo)

# ╔═╡ c557eeb4-fab4-4751-af4b-bea2115e85a9
lsd_smth = getLandSea(e5ds,egeo,smooth=true,σlon=1,σlat=1)

# ╔═╡ dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]; y = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ a9bae1a3-c994-4b36-b9f8-740aa14a3929
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=1,axwidth=1.5,ncols=2)

	lvls = -10:0; lvls = lvls[lvls.<log10(0.5)]; lvls = vcat(lvls,log10(0.5))
	lsm = log10.(lsd_raw.lsm)
	lsm[lsm.==-Inf] .= -99
	a1[1].pcolormesh(
		lsd_raw.lon,lsd_raw.lat,lsm',
		levels=lvls,extend="both",cmap="Delta"
	)
	a1[2].format(xlim=(270,285),ylim=(0,15))
	
	c = a1[2].pcolormesh(
		lsd_smth.lon,lsd_smth.lat,log10.(lsd_smth.lsm)',
		levels=lvls,extend="both",cmap="Delta"
	)
	a1[2].format(xlim=(270,285),ylim=(0,15))
	a1[2].colorbar(c,loc="r")
	
	f1.savefig(plotsdir("01a-lsmraw.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("01a-lsmraw.png"))
end

# ╔═╡ Cell order:
# ╟─faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
# ╟─a675084e-f638-11eb-2930-b70006b9445c
# ╟─4e08778c-cd46-4a79-9b6a-96ca6447e975
# ╟─75ea588a-64df-4b31-9dae-f99f425ee55a
# ╟─0d950c5d-8c29-4eed-9716-e7b979c59ff6
# ╠═79e35867-824c-4b3a-8c7b-3abb8f9d9761
# ╟─f3773f6f-55cc-4546-92c7-dd3d8955edab
# ╟─c557eeb4-fab4-4751-af4b-bea2115e85a9
# ╟─dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
# ╟─a9bae1a3-c994-4b36-b9f8-740aa14a3929
