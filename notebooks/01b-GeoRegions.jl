### A Pluto.jl notebook ###
# v0.19.26

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
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using PlutoUI

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01b. Creating GeoRegions

In this notebook, we define GeoRegions based on the filtered Land-Sea mask we created in notebook `01a-lsmfiltering.jl`.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Datasets and LandSea Masks ...
"

# ╔═╡ 189e1048-c92d-457e-a30e-f4e523b80afc
begin
	infoall = stninfoall()
	infody  = stninfody()
	infomo  = stninfomo()
	md"Loading station location information ..."
end

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 1cee44fe-75b7-42a2-948e-db330cf788e8
npd = IMERGDummy(path=datadir())

# ╔═╡ 0a757a59-88af-4838-9696-f0b0ab776c5a
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC")

# ╔═╡ 5d6c3cd6-e406-461d-a226-20022060398d
begin
	lsd = getLandSea(
		geo,path=datadir(),savelsd=true,returnlsd=true,
		smooth=true,σlon=30,σlat=30
	)
end

# ╔═╡ 9c6f002e-9146-449c-9d4b-ab2415c98785
elsd = getLandSea(e5ds,ERA5Region(geo),smooth=true,σlon=2,σlat=2)

# ╔═╡ a01e1721-23cf-4a3f-b5aa-0189b1a113a3
nlsd = getLandSea(npd,geo,smooth=true,σlon=5,σlat=5)

# ╔═╡ e904cdc1-0846-4a26-a120-41b4a22b4102
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=3,axwidth=2,wspace=1.5)

	lvls = vcat(10. .^(-5:0.5:-1))
	lvls = vcat(lvls,0.25,0.5,0.9,0.95,0.99)
	lvls = log10.(lvls)

	lon = lsd.lon.-360
	lat = lsd.lat
	lsm = log10.(lsd.lsm)

	cmap_dict = Dict("right"=>0.62)
	
	c = a1[1].pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a1[2].pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a1[3].pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	
	a1[1].scatter(infomo[:,2],infomo[:,3],c="pink")
	a1[1].scatter(infody[:,2],infody[:,3],c="red")
	a1[2].scatter(infody[:,2],infody[:,3],c="red")
	a1[3].scatter(infomo[:,2],infomo[:,3],c="pink")

	a1[1].format(ultitle="(a) All")
	a1[2].format(ultitle="(b) Daily")
	a1[3].format(ultitle="(c) Monthly")

	for ax in a1
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(
			xlim=(-90,-75),ylim=(0,15),suptitle="Available Columbia Stations",
			xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
		)
	end

	f1.colorbar(c,ticklabels=["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.95","0.99"],label="Land-Sea Mask")
	f1.savefig(plotsdir("01b-stations.png"),transparent=false,dpi=300)
	load(plotsdir("01b-stations.png"))
end

# ╔═╡ e8c67f38-31aa-4855-b3fa-f24da732bc2b
md"
We see that most of the stations with daily data are located near the coast or near-coast inland (<0.99).  However, most of the monthly stations within the domain are located farther inland except for San Andreas, Cartagena and Tulenapa.
"

# ╔═╡ f252a060-111a-4e86-85dd-46257c258b77
md"
### B. Defining the Pacific and Atlantic Regions
"

# ╔═╡ ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
begin
	addGeoRegions(srcdir("OTRECGeoRegions.txt"))
	md"Defining main OTREC GeoRegion from file ..."
end

# ╔═╡ 47127cc3-442e-4fa7-b016-82ca8b6ab87a
begin
	if isGeoRegion("OTREC_PAC",throw=false)
		removeGeoRegion("OTREC_PAC")
	end
	PAC = PolyRegion(
		"OTREC_PAC","OTREC","OTREC Pacific Region",
		[-78,-90,-90,-86,-84,-82,-81,-80,-79.5,-78.5,-76,-76,-78],
		[0,0,15,14,10,8.5,8.5,9,9.3,9,7,3,0]
	)
end

# ╔═╡ 3fe4ff94-41c3-46e9-bed7-7ed6cb767d95
begin
	if isGeoRegion("OTREC_ATL",throw=false)
		removeGeoRegion("OTREC_ATL")
	end
	ATL = PolyRegion(
		"OTREC_ATL","OTREC","OTREC Atlantic Region",
		[-75,-75,-86,-86,-84,-82,-81,-80,-79.5,-78.5,-76,-75],
		[10,15,15,14,10,8.5,8.5,9,9.3,9,7,10]
	)
end

# ╔═╡ 7ca2f79d-550c-4d3f-9755-13d5f33129c9
begin
	if isGeoRegion("OTREC_CST",throw=false)
		removeGeoRegion("OTREC_CST")
	end
	CST = RectRegion(
		"OTREC_CST","OTREC","OTREC Colombia Coast",[7,3,-76.5,-78.5]
	)
end

# ╔═╡ 3de6247e-e679-49b8-aef6-967308e13d9d
begin
	if isGeoRegion("OTREC_SAN",throw=false)
			removeGeoRegion("OTREC_SAN")
	end
	SAN = RectRegion(
			"OTREC_SAN","OTREC","OTREC San Andres",
			[13.1,12.1,-81.2,-82.2],
	)
end

# ╔═╡ e51ebe7b-99eb-450d-8603-cd4811e7b311
begin
	if isGeoRegion("OTREC_OCT",throw=false)
			removeGeoRegion("OTREC_OCT")
	end
	OCT = RectRegion(
			"OTREC_OCT","OTREC","OTREC Colombia Ocean Coast",
			[7,3,-77.5,-78.5],
	)
end

# ╔═╡ eaca1f87-36f3-4f27-a1a7-9ae49e4edf8d
begin
	if isGeoRegion("OTREC_LCT",throw=false)
			removeGeoRegion("OTREC_LCT")
	end
	LCT = RectRegion(
			"OTREC_LCT","OTREC","OTREC Colombia Land Coast",
			[7,3,-76.5,-77.5],
	)
end

# ╔═╡ cab91dff-b19e-4d4e-940c-6b12e00344d4
begin
	if isGeoRegion("OTREC_PCS",throw=false)
		removeGeoRegion("OTREC_PCS")
	end
		PCS = RectRegion("OTREC_PCS","OTREC","OTREC Pacific Coast",[7,2,-75,-83])
end

# ╔═╡ 601f4ad5-3520-4555-af6b-8eb2a299609c
begin
	blon_PAC,blat_PAC,slon_PAC,slat_PAC = coordGeoRegion(PAC)
	blon_ATL,blat_ATL,slon_ATL,slat_ATL = coordGeoRegion(ATL)
	slon_PCS,slat_PCS = coordGeoRegion(PCS)
	slon_CST,slat_CST = coordGeoRegion(CST)
	slon_SAN,slat_SAN = coordGeoRegion(SAN)
	md"Loading geographic shapes of the GeoRegions ..."
end

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=3,axwidth=1.75)

	a2[1].pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a2[2].pcolormesh(elsd.lon.-360,elsd.lat,log10.(elsd.lsm)',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a2[3].pcolormesh(nlsd.lon,nlsd.lat,log10.(nlsd.lsm)',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	
	a2[3].scatter(infomo[:,2],infomo[:,3],zorder=4,c="pink",label="Monthly",legend="r",legend_kw=Dict("ncol"=>1,"frame"=>false))
	a2[3].scatter(infody[:,2],infody[:,3],zorder=4,c="red",label="Daily",legend="r")

	a2[1].format(ltitle=L"(a) ETOPO ($\sigma = 30$)")
	a2[2].format(ltitle=L"(b) ERA5 ($\sigma = 2$)")
	a2[3].format(ltitle=L"(a) IMERG ($\sigma = 5$)")

	for ax in a2
		ax.plot(x,y,lw=0.5,c="k")
		ax.scatter(infomo[:,2],infomo[:,3],zorder=4,c="pink")
		ax.scatter(infody[:,2],infody[:,3],zorder=4,c="red")
		ax.plot(slon_PAC,slat_PAC,lw=5,)
		ax.plot(slon_ATL,slat_ATL,lw=5,)
		# ax.plot(slon_CST,slat_CST,lw=5,)
		ax.plot(slon_PCS,slat_PCS,lw=3,)
		ax.plot(slon_SAN,slat_SAN,lw=3,)
		ax.text(-88.5,1.5,"PAC",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(-78,13,"ATL",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        # ax.text(-81,4.5,"CST",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(-84.7,11.1,"SAN",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(-84,1.5,"PCS",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
		ax.format(
			xlim=(-90,-75),ylim=(0,15),
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$",
			suptitle="Available Colombia Stations",
			grid=true
		)
	end

	f2.colorbar(loc="l",c,ticklabels=["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.95","0.99"],label="Land-Sea Mask")
	f2.savefig(plotsdir("01b-georegions.png"),transparent=false,dpi=300)
	load(plotsdir("01b-georegions.png"))
end

# ╔═╡ 29800e11-610b-4271-969c-26302d9d0432
md"
We also defined a coastal subGeoRegion for the Pacific GeoRegion where precipitation in the Eastern Pacific is most highly concentrated just off the coast of Colombia
"

# ╔═╡ 7c0d3874-44c1-4e17-83b5-b0aabcb317cd
md"
### C. Test Extraction of Coast in GeoRegions

First, let us define the GeoRegion to test extraction ...
"

# ╔═╡ 36b7ec91-6409-4828-9548-1eafc0559348
geo_tst = GeoRegion("OTREC_PCS")

# ╔═╡ 423a607a-4357-44df-9b7e-948aa8062baf
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo_tst,lon,lat)
	glon = ggrd.lon; ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	glat = ggrd.lat; ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	rlsm = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
			mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
	end
	# rlsm[(rlsm.>=0.9).|(rlsm.<=0.2)] .= NaN
	md"Extracting information for region and restricting to coast ..."
end

# ╔═╡ e6f6cf23-9b6b-49d2-bf75-bf2826618a4e
begin
	pplt.close(); f3,a3 = pplt.subplots(ncols=1,axwidth=2,wspace=1.5)

	a3[1].pcolormesh(glon,glat,rlsm',levels=vcat(lvls,1),cmap="delta",extend="both",cmap_kw=cmap_dict)
	a3[1].scatter(infomo[:,2],infomo[:,3],zorder=4,c="pink")
	a3[1].scatter(infody[:,2],infody[:,3],zorder=4,c="red")

	for ax in a3
		ax.plot(x,y,lw=0.5,c="k")
		# ax.plot(slon_PAC,slat_PAC,lw=5,)
		# ax.plot(slon_ATL,slat_ATL,lw=5,)
		# ax.plot(slon_CST,slat_CST,lw=5,)
		# ax.text(-88.5,1.5,"PAC",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
  #       ax.text(-78,13,"ATL",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
  #       ax.text(-81.5,4.5,"CST",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
		ax.format(
			xlim=(-90,-75),ylim=(0,15),
			suptitle="Available Columbia Stations",
			grid=true
		)
	end

	f3.colorbar(c,["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.99"],label="ETOPO Land-Sea Mask")
	f3.savefig(plotsdir("01b-testcoast.png"),transparent=false,dpi=300)
	load(plotsdir("01b-testcoast.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─1cee44fe-75b7-42a2-948e-db330cf788e8
# ╟─0a757a59-88af-4838-9696-f0b0ab776c5a
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╟─9c6f002e-9146-449c-9d4b-ab2415c98785
# ╟─a01e1721-23cf-4a3f-b5aa-0189b1a113a3
# ╟─e904cdc1-0846-4a26-a120-41b4a22b4102
# ╟─e8c67f38-31aa-4855-b3fa-f24da732bc2b
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╟─ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
# ╟─47127cc3-442e-4fa7-b016-82ca8b6ab87a
# ╟─3fe4ff94-41c3-46e9-bed7-7ed6cb767d95
# ╟─7ca2f79d-550c-4d3f-9755-13d5f33129c9
# ╟─3de6247e-e679-49b8-aef6-967308e13d9d
# ╟─e51ebe7b-99eb-450d-8603-cd4811e7b311
# ╟─eaca1f87-36f3-4f27-a1a7-9ae49e4edf8d
# ╠═cab91dff-b19e-4d4e-940c-6b12e00344d4
# ╟─601f4ad5-3520-4555-af6b-8eb2a299609c
# ╟─d7755534-3565-4011-b6e3-e131991008db
# ╟─29800e11-610b-4271-969c-26302d9d0432
# ╟─7c0d3874-44c1-4e17-83b5-b0aabcb317cd
# ╟─36b7ec91-6409-4828-9548-1eafc0559348
# ╟─423a607a-4357-44df-9b7e-948aa8062baf
# ╟─e6f6cf23-9b6b-49d2-bf75-bf2826618a4e
