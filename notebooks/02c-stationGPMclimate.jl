### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 9db63111-71ec-43b9-9174-d4921191eafd
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ da16ef31-2091-4c36-b3d0-05f4197771f6
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ e0a4eb5e-467e-11ec-21da-af325ae0fd15
md"
# 02c. GPM Climatology for various OTREC Stations

Tex
"

# ╔═╡ a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
md"
### A. Defining the NASAPrecipitation Dataset
"

# ╔═╡ 5dfc4918-a20b-422e-b2bd-fa4805bd3f82
npd = IMERGMonthly(start=Date(2001),stop=Date(2020),path=datadir())

# ╔═╡ 8045657b-83e3-40be-a9ee-e5875d48f73f
begin
	dt = npd.start : Year(1) : npd.stop
	nt = length(dt)
	md"Constructing date vectors ..."
end

# ╔═╡ c05d8e84-c595-42d5-934a-1a2a91f85d97
geo = GeoRegion("OTREC")

# ╔═╡ 665575fe-b149-410c-9bae-53dedf1ec264
lsd = getLandSea(npd,geo)

# ╔═╡ 0ec50c3d-55b1-4015-87ff-30f76412cd96
md"
### B. Loading Station Information
"

# ╔═╡ 704c7bcf-82cb-4662-aecf-540369e8a11b
begin
	infody = stninfody(); nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ 4ddeb466-e4af-4ef2-8328-121289fa9b20
begin
	pdata = zeros(12,nstn)
	md"Preallocating arrays ..."
end

# ╔═╡ 450d68a4-5f3b-4c52-8847-d882787dbc0c
begin
	ilon = zeros(Int,nstn)
	ilat = zeros(Int,nstn)
	for istn = 1 : nstn
		ilon[istn] = argmin(abs.(lsd.lon.-infody[istn,2]))
		ilat[istn] = argmin(abs.(lsd.lat.-infody[istn,3]))
	end
	md"Finding nearest longitude/latitude coordinate points"
end

# ╔═╡ d4945e6a-4aaf-499c-acc9-894b57e04563
infody

# ╔═╡ 8b509f24-da2a-40fd-a69a-d5e1c7356d11
md"
### C. Loading the Precipitation Data
"

# ╔═╡ e095a808-38e9-4015-9224-6295f75470dd
begin
	for tt in dt
		ids = read(npd,geo,tt)
		for istn = 1 : nstn
			pdata[:,istn] .+= ids["precipitation"][ilon[istn],ilat[istn],:] * 3600 / 20
		end
		close(ids)
	end
	md"Loading precipitation data ..."
end

# ╔═╡ 9139653a-d357-456e-bda8-440e2f088403
begin
	pplt.close(); ft,at = pplt.subplots(aspect=2,axwidth=3)

	at[1].plot(0.5:12,pdata[:,1])
	at[1].plot(0.5:12,pdata[:,2:4],c="r")
	at[1].plot(0.5:12,pdata[:,[9,10,12]],c="k")
	at[1].plot(0.5:12,pdata[:,11],c="g")
	at[1].plot(0.5:12,pdata[:,8],c="b")
	at[1].plot(0.5:12,pdata[:,[5,6,7]],c="y")
	at[1].format(ylim=(0,20))
	
	ft.savefig(plotsdir("02c-gpmclimate.png"),transparent=false,dpi=200)
	load(plotsdir("02c-gpmclimate.png"))
end

# ╔═╡ Cell order:
# ╟─e0a4eb5e-467e-11ec-21da-af325ae0fd15
# ╟─9db63111-71ec-43b9-9174-d4921191eafd
# ╟─da16ef31-2091-4c36-b3d0-05f4197771f6
# ╟─a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
# ╟─5dfc4918-a20b-422e-b2bd-fa4805bd3f82
# ╟─8045657b-83e3-40be-a9ee-e5875d48f73f
# ╟─c05d8e84-c595-42d5-934a-1a2a91f85d97
# ╟─665575fe-b149-410c-9bae-53dedf1ec264
# ╟─0ec50c3d-55b1-4015-87ff-30f76412cd96
# ╟─704c7bcf-82cb-4662-aecf-540369e8a11b
# ╟─4ddeb466-e4af-4ef2-8328-121289fa9b20
# ╟─450d68a4-5f3b-4c52-8847-d882787dbc0c
# ╠═d4945e6a-4aaf-499c-acc9-894b57e04563
# ╟─8b509f24-da2a-40fd-a69a-d5e1c7356d11
# ╟─e095a808-38e9-4015-9224-6295f75470dd
# ╠═9139653a-d357-456e-bda8-440e2f088403
