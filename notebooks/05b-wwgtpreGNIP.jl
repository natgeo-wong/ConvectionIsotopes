### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ fba33adb-d8a9-495e-b927-b9b62aaa74ef
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9149a006-add6-467e-9ad1-be497a53fff7
begin
	@quickactivate "ColombiaIsotope"
	using DataFrames
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using Statistics
	using XLSX

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ 8e744f98-21b4-11ed-017a-8571e8807e61
md"
# 05b. W-Weighted Mean Pressure for GNIP
"

# ╔═╡ 3b799f5e-0789-49cb-a1c7-895edbf9c7e5
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	clon = coast[:,1]; clat = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 65059102-8d7c-4481-8931-2603480ba231
md"
### A. Defining Datasets and Regions
"

# ╔═╡ 2e783000-8fb9-4e48-9471-16ecf4780596
e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 0ffb4e9f-41e4-4714-8b3d-8795f0df9f5c
evar = SingleVariable("p_wwgt")

# ╔═╡ 42381a85-bf1b-4fcf-b50c-2c590e94f99c
geo = GeoRegion("TRP")

# ╔═╡ e160b2ed-5c60-4156-876c-b2a5bfde3631
ereg = ERA5Region(geo)

# ╔═╡ 95d5d27f-69c2-405b-8376-1d02ae52a5d1
npd = IMERGMonthly(start=Date(2001),stop=Date(2020),path=datadir())

# ╔═╡ 8dca0258-9a69-4770-bf00-e8919c18e0df
md"
### B. Loading Datasets
"

# ╔═╡ 9f498a11-07d4-4c13-b932-4679f1056978
begin
	xf = XLSX.readxlsx(datadir("GNIP-SEA.xlsx"))["Data"]["A1:I973"]
	md"Loading GNIP data ..."
end

# ╔═╡ b671a99c-9824-4fb8-aca4-1ac35fbeb3eb
lsd_npd = getIMERGlsd(geo,path=datadir("imergmask"))

# ╔═╡ d4185fbd-d8e1-4c1c-9119-c22841462970
lsd_era = ERA5Reanalysis.getLandSea(ereg,path=datadir("emask"))

# ╔═╡ ef55c91a-a9cb-4749-89b1-fe9c7f2fd13c
begin
	ds  = NCDataset(datadir("flsm","flsm-TRP.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	lsm = ds["flsm"][:]
	close(ds)
	md"Loading filtered ERA5 land-sea mask ..."
end

# ╔═╡ d2001450-fe8b-4412-b2d3-6f659a737f9d
md"
### C. Precipitation Data Comparison
"

# ╔═╡ b3f3d984-d9b4-4869-946c-1a683f2e051e
begin
	gnip_prcp = xf[2:end,end]
	gnip_prcp[ismissing.(gnip_prcp)] .= NaN
	gnip_prcp = convert(Vector{Float64},gnip_prcp) / 30
	gnip_δ18O = xf[2:end,end-2]
	gnip_δ18O[ismissing.(gnip_δ18O)] .= NaN
	gnip_δ18O = convert(Vector{Float64},gnip_δ18O)
	npd_prcp  = zeros(length(gnip_prcp)) * NaN
	era5_wgtp = zeros(length(gnip_prcp)) * NaN
	md"Filtering out relevant GNIP data and preallocating arrays for precipiation and ERA5 data ..."
end

# ╔═╡ 7025fc11-43b2-46c1-9a30-e4bbda22563b
for idata = 1 : length(npd_prcp)
	idt = Date(xf[idata+1,6])
	if idt >= Date(2001)
		imo = month(idt)
		ilon = argmin(abs.(lsd_era.lon .- xf[idata+1,4]))
		ilat = argmin(abs.(lsd_era.lat .- xf[idata+1,3]))
		erads = read(e5ds,evar,ereg,idt)
		iip = erads["p_wwgt"][ilon,ilat,imo]
		if !ismissing(iip)
			era5_wgtp[idata] = iip / 100
		else; era5_wgtp[idata] = NaN
		end
		close(erads)
	end
end

# ╔═╡ 79bf1a50-ad23-418f-9b82-f717b01d4b23
begin
	rainb = 2.5  : 2.5 : 97.5; nr = length(rainb)
	rainu = 5 : 2.5 : 100
	rainm = 2.5  : 2.5 : 100
	# rainb = 5  : 5 : 95; nr = length(rainb)
	# rainu = 10 : 5 : 100
	# rainm = 5  : 5 : 100
	pwgtb = 100 : 25 : 925; np = length(pwgtb)
	pwgtu = 125 : 25 : 950
	pwgtm = 100 : 25 : 950'
	md"Defining ranges for binning of isotopes and rainfall ..."
end

# ╔═╡ f41e7540-da7d-4806-afaa-83cd7cf7af43
begin
	binmat = zeros(nr,np)
	binmt2 = zeros(nr,np)
	binstd = zeros(nr,np)
	binnum = zeros(nr,np)
	for ip = 1 : np, ir = 1 : nr
		ibin = gnip_δ18O[(era5_wgtp.>pwgtb[ip]) .& (era5_wgtp.<pwgtu[ip]) .& (gnip_prcp.>rainb[ir]) .& (gnip_prcp.<rainu[ir])]
		if !isempty(ibin)
			binmat[ir,ip] = mean(ibin[.!isnan.(ibin)])
			binnum[ir,ip] = sum(.!isnan.(ibin))
		else
			binmat[ir,ip] = NaN
			binnum[ir,ip] = NaN
		end
	end
	for ip = 1 : np, ir = 1 : nr
		ibin = gnip_δ18O[(era5_wgtp.>pwgtb[ip]) .& (era5_wgtp.<pwgtu[ip]) .& (gnip_prcp.>rainb[ir]) .& (gnip_prcp.<rainu[ir])]
		if !isempty(ibin)
			binstd[ir,ip] = sum((ibin[.!isnan.(ibin)] .- binmat[ir,ip]).^2)
		else
			binstd[ir,ip] = NaN
		end
	end
	binmt2 = deepcopy(binmat)
	binmt2[binnum.<2] .= NaN
	binstd[binnum.<2]  .= NaN
	md"Finding mean and standard deviation for station isotopic data ..."
end

# ╔═╡ 83209cf1-3cf0-4dc8-9056-bc0240dc3311
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=2,aspect=1/2,axwidth=1.5)

	# binmat[binnum.<2] .= NaN
	binstd[binnum.<2] .= NaN
	a1[1].pcolormesh(
		rainm,pwgtm,binmat',alpha=0.2,
		cmap="viridis",levels=-15:-5,extend="both"
	)
	
	c1_1 = a1[1].pcolormesh(
		rainm,pwgtm,binmt2',
		cmap="viridis",levels=-15:-5,extend="both"
	)
	a1[1].format(urtitle="(a)")
	
	c1_2 = a1[2].pcolormesh(
		rainm,pwgtm,sqrt.(binstd./binnum)',
		levels=0:0.5:4,extend="both"
	)
	a1[2].format(urtitle="(b)")

	for ax in a1
		ax.format(
			ylim=(950,100),xlim=(0,50),ylabel=L"p$_w$ / hPa",
			xlabel=L"Rain Rate / mm day$^{-1}$",yminorlocator=100:25:950
		)
	end

	lbl = L"$\delta^{18}$O / $\perthousand$"
	a1[1].colorbar(c1_1,loc="b",label= L"$\mu$ " * lbl,locator=-15:5:-5)
	a1[2].colorbar(c1_2,loc="b",label=L"$\sigma$ " * lbl,locator=0:4)
	f1.savefig(plotsdir("05b-wwgtprevrcp_GNIP.png"),transparent=false,dpi=150)
	load(plotsdir("05b-wwgtprevrcp_GNIP.png"))
end

# ╔═╡ Cell order:
# ╟─8e744f98-21b4-11ed-017a-8571e8807e61
# ╟─fba33adb-d8a9-495e-b927-b9b62aaa74ef
# ╟─9149a006-add6-467e-9ad1-be497a53fff7
# ╟─3b799f5e-0789-49cb-a1c7-895edbf9c7e5
# ╟─65059102-8d7c-4481-8931-2603480ba231
# ╟─2e783000-8fb9-4e48-9471-16ecf4780596
# ╟─0ffb4e9f-41e4-4714-8b3d-8795f0df9f5c
# ╟─42381a85-bf1b-4fcf-b50c-2c590e94f99c
# ╟─e160b2ed-5c60-4156-876c-b2a5bfde3631
# ╟─95d5d27f-69c2-405b-8376-1d02ae52a5d1
# ╟─8dca0258-9a69-4770-bf00-e8919c18e0df
# ╟─9f498a11-07d4-4c13-b932-4679f1056978
# ╟─b671a99c-9824-4fb8-aca4-1ac35fbeb3eb
# ╟─d4185fbd-d8e1-4c1c-9119-c22841462970
# ╟─ef55c91a-a9cb-4749-89b1-fe9c7f2fd13c
# ╟─d2001450-fe8b-4412-b2d3-6f659a737f9d
# ╟─b3f3d984-d9b4-4869-946c-1a683f2e051e
# ╟─7025fc11-43b2-46c1-9a30-e4bbda22563b
# ╟─79bf1a50-ad23-418f-9b82-f717b01d4b23
# ╟─f41e7540-da7d-4806-afaa-83cd7cf7af43
# ╟─83209cf1-3cf0-4dc8-9056-bc0240dc3311
