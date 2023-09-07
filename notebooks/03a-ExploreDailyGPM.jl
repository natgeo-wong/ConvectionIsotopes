### A Pluto.jl notebook ###
# v0.19.26

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
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ e0a4eb5e-467e-11ec-21da-af325ae0fd15
md"
# 02a - Exploring the Daily GPM Data over the OTREC Region

In this notebook, we explore daily GPM data over the OTREC Region.  In contrast to notebook `01a` which explores monthly climatology, here we explore the daily data over the rough OTREC region from 2020 to 2021 when there is station data that is of greater than monthly frequency.
"

# ╔═╡ a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
md"
### A. Defining the NASAPrecipitation Dataset
"

# ╔═╡ 5dfc4918-a20b-422e-b2bd-fa4805bd3f82
npd = IMERGFinalHH(start=Date(2020,1,1),stop=Date(2020,11,30),path=datadir())

# ╔═╡ a6f9833e-5349-4aa7-aaa5-e735876b8581
begin
	dt = npd.start : Day(1) : (npd.stop + Day(1))
	nt = length(dt)
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

# ╔═╡ 38941d48-ccf2-464e-a90f-2f0306b5d24f
begin
	tnc = joinpath(
		npd.datapath,geo.regID,"2020","01","imergfinalhh-OTREC-20200101.nc"
	)
	tds  = NCDataset(tnc)
	lon  = tds["longitude"][:]
	lat  = tds["latitude"][:]
	close(tds)
	md"Loading coordinate information ..."
end

# ╔═╡ 4ddeb466-e4af-4ef2-8328-121289fa9b20
begin
	pdata = zeros(48,nt,nstn)
	md"Preallocating arrays ..."
end

# ╔═╡ 450d68a4-5f3b-4c52-8847-d882787dbc0c
begin
	ilon = zeros(Int,nstn)
	ilat = zeros(Int,nstn)
	for istn = 1 : nstn
		ilon[istn] = argmin(abs.(lon.-infody[istn,2]))
		ilat[istn] = argmin(abs.(lat.-infody[istn,3]))
	end
	md"Finding nearest longitude/latitude coordinate points"
end

# ╔═╡ 8b509f24-da2a-40fd-a69a-d5e1c7356d11
md"
### C. Loading the Precipitation Data
"

# ╔═╡ 36b27e6e-a32b-457d-baa2-5c11408d23d3
tshift = 10

# ╔═╡ e095a808-38e9-4015-9224-6295f75470dd
begin
	for it = 1 : nt
		idt = dt[it]
		ifo = joinpath(npd.datapath,geo.regID,Dates.format(idt,dateformat"yyyy/mm"))
		inc = npd.npdID * "-" * geo.regID * "-" * Dates.format(idt,dateformat"yyyymmdd") * ".nc"
		ids = NCDataset(joinpath(ifo,inc))
		for istn = 1 : nstn
			pdata[:,it,istn] .= ids["precipitation"][ilon[istn],ilat[istn],:] * 3600
		end
	end
	npdata = reshape(pdata,:,nstn)
	npdata = npdata[(1+tshift):(end+tshift-48),:]
	npdata = dropdims(mean(reshape(npdata,48,:,nstn),dims=1),dims=1)
	nndata = dropdims(mean(pdata,dims=1),dims=1)[1:(end-1),:]
	md"Loading precipitation data ..."
end

# ╔═╡ 9139653a-d357-456e-bda8-440e2f088403
begin
	pplt.close(); ft,at = pplt.subplots(aspect=1)
	at[1].scatter(nndata,npdata)
	at[1].format(
		xlabel=L"Daily Precipitation Rate, GMT+0 / mm hr$^{-1}$",
		ylabel=L"Daily Precipitation Rate, GMT-5 / mm hr$^{-1}$",
		suptitle="The Effects of Time-shifting"
	)
	ft.savefig(plotsdir("03a-timeshift.png"),transparent=false,dpi=200)
	load(plotsdir("03a-timeshift.png"))
end

# ╔═╡ 03c3f89e-3699-45fa-ae75-f13fa71daddf
md"
### D. Loading spatial precipitation data ...
"

# ╔═╡ d107212e-5a1a-4211-98a4-b647f909b561
begin
	const data = zeros(length(lsd.lon),length(lsd.lat))
	for it = 1 : nt
		idt = dt[it]
		ifo = joinpath(npd.datapath,geo.regID,Dates.format(idt,dateformat"yyyy/mm"))
		inc = npd.npdID * "-" * geo.regID * "-" * Dates.format(idt,dateformat"yyyymmdd") * ".nc"
		ids = NCDataset(joinpath(ifo,inc))
		data[:,:] += mean(ids["precipitation"][:],dims=3) * 86400
		close(ids)
	end
	data /= nt
	md"Loading precipitation data ..."
end

# ╔═╡ aa803db6-b1c4-42ec-a26d-83a9e0998fc7
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=1,axwidth=1.5)
	
	axs[1].pcolormesh(lsd.lon,lsd.lat,data')
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─e0a4eb5e-467e-11ec-21da-af325ae0fd15
# ╟─9db63111-71ec-43b9-9174-d4921191eafd
# ╟─da16ef31-2091-4c36-b3d0-05f4197771f6
# ╟─a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
# ╟─5dfc4918-a20b-422e-b2bd-fa4805bd3f82
# ╟─a6f9833e-5349-4aa7-aaa5-e735876b8581
# ╟─c05d8e84-c595-42d5-934a-1a2a91f85d97
# ╟─665575fe-b149-410c-9bae-53dedf1ec264
# ╟─0ec50c3d-55b1-4015-87ff-30f76412cd96
# ╟─704c7bcf-82cb-4662-aecf-540369e8a11b
# ╟─38941d48-ccf2-464e-a90f-2f0306b5d24f
# ╟─4ddeb466-e4af-4ef2-8328-121289fa9b20
# ╟─450d68a4-5f3b-4c52-8847-d882787dbc0c
# ╟─8b509f24-da2a-40fd-a69a-d5e1c7356d11
# ╠═36b27e6e-a32b-457d-baa2-5c11408d23d3
# ╟─e095a808-38e9-4015-9224-6295f75470dd
# ╟─9139653a-d357-456e-bda8-440e2f088403
# ╟─03c3f89e-3699-45fa-ae75-f13fa71daddf
# ╠═d107212e-5a1a-4211-98a4-b647f909b561
# ╠═aa803db6-b1c4-42ec-a26d-83a9e0998fc7
