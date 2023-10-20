### A Pluto.jl notebook ###
# v0.19.27

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
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using PlutoUI
	using Printf

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01a. Creating GeoRegions

In this notebook, we define additional GeoRegions of interest in addition to those given in `src/OTRECRectRegions` and `src/OTRECPolyRegions`.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Basic Datasets and Station Information
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

# ╔═╡ f252a060-111a-4e86-85dd-46257c258b77
md"
### B. Default GeoRegions for the OTREC Domain
"

# ╔═╡ ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
begin
	addGeoRegions(srcdir("OTRECRectRegions.txt"))
	addGeoRegions(srcdir("OTRECPolyRegions.txt"))
	md"Defining main OTREC GeoRegion from file ..."
end

# ╔═╡ d5bc84d7-6516-4fad-a245-703e14e22f77
begin
	PAC = GeoRegion("OTREC_PAC")
	ATL = GeoRegion("OTREC_ATL")
	PCS = GeoRegion("OTREC_PCS")
	SAN = GeoRegion("OTREC_SAN")
	md"Loading GeoRegions ..."
end

# ╔═╡ 601f4ad5-3520-4555-af6b-8eb2a299609c
begin
	blon_PAC,blat_PAC,slon_PAC,slat_PAC = coordGeoRegion(PAC)
	blon_ATL,blat_ATL,slon_ATL,slat_ATL = coordGeoRegion(ATL)
	slon_PCS,slat_PCS = coordGeoRegion(PCS)
	slon_SAN,slat_SAN = coordGeoRegion(SAN)
	md"Loading geographic shapes of the GeoRegions ..."
end

# ╔═╡ bbaa6dea-6db3-4ff9-8cd8-047c8941fd0f
md"### C. Loading and Plotting the Domain LandSea Masks"

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

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=3,axwidth=1.75)

	lvls = vcat(10. .^(-5:0.5:-1))
	lvls = vcat(lvls,0.25,0.5,0.9,0.95,0.99)
	cmap_dict = Dict("right"=>0.65)
	
	c = a2[1].pcolormesh(lsd.lon,lsd.lat,lsd.lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a2[2].pcolormesh(elsd.lon,elsd.lat,elsd.lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	a2[3].pcolormesh(nlsd.lon,nlsd.lat,nlsd.lsm',levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict)
	
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
		ax.plot(slon_PCS,slat_PCS,lw=3,)
		ax.plot(slon_SAN,slat_SAN,lw=3,)
		ax.text(271.5,1.5,"PAC",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(282,13,"ATL",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(275.3,11.1,"SAN",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
        ax.text(276,1.5,"PCS",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
		ax.format(
			xlim=(270,285),ylim=(0,15),
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$",
			suptitle="Available Colombia Stations",
			grid=true
		)
	end

	f2.colorbar(loc="l",c,ticklabels=["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.95","0.99"],label="Land-Sea Mask")
	f2.savefig(plotsdir("01a-georegions.png"),transparent=false,dpi=300)
	load(plotsdir("01a-georegions.png"))
end

# ╔═╡ e8c67f38-31aa-4855-b3fa-f24da732bc2b
md"
We see that most of the stations with daily data are located near the coast or near-coast inland (<0.99).  However, most of the monthly stations within the domain are located farther inland except for San Andreas, Cartagena and Tulenapa.
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
	ggrd = RegionGrid(geo_tst,lsd.lon,lsd.lat)
	glon = ggrd.lon; ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	glat = ggrd.lat; ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	rlsm = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
			mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		rlsm[glon,glat] = lsd.lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
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
			xlim=(270,285),ylim=(0,15),
			suptitle="Available Columbia Stations",
			grid=true
		)
	end

	f3.colorbar(c,["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.99"],label="ETOPO Land-Sea Mask")
	f3.savefig(plotsdir("01a-testcoast.png"),transparent=false,dpi=300)
	load(plotsdir("01a-testcoast.png"))
end

# ╔═╡ 58e4153d-a8a3-454e-8941-389e9e5413a6
md"
### D. Define Domain around Stations of Interest
"

# ╔═╡ ecbf4852-7d65-44ca-94aa-d811e2dafb3b
for istn = 1 : 12

	istn_lon = infody[istn,2]
	istn_lat = infody[istn,3]

	if isGeoRegion("OTREC_STN$(@sprintf("%02d",istn))",throw=false)
		removeGeoRegion("OTREC_STN$(@sprintf("%02d",istn))")
	end
	RectRegion(
		"OTREC_STN$(@sprintf("%02d",istn))","OTREC",
		"Daily Station $(istn)",[
			infody[istn,3]+0.5,infody[istn,3]-0.5,
			infody[istn,2]+0.5,infody[istn,2]-0.5
		]
	)
	
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
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╟─ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
# ╟─d5bc84d7-6516-4fad-a245-703e14e22f77
# ╟─601f4ad5-3520-4555-af6b-8eb2a299609c
# ╟─bbaa6dea-6db3-4ff9-8cd8-047c8941fd0f
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╟─9c6f002e-9146-449c-9d4b-ab2415c98785
# ╟─a01e1721-23cf-4a3f-b5aa-0189b1a113a3
# ╟─d7755534-3565-4011-b6e3-e131991008db
# ╟─e8c67f38-31aa-4855-b3fa-f24da732bc2b
# ╟─7c0d3874-44c1-4e17-83b5-b0aabcb317cd
# ╟─36b7ec91-6409-4828-9548-1eafc0559348
# ╟─423a607a-4357-44df-9b7e-948aa8062baf
# ╟─e6f6cf23-9b6b-49d2-bf75-bf2826618a4e
# ╟─58e4153d-a8a3-454e-8941-389e9e5413a6
# ╟─ecbf4852-7d65-44ca-94aa-d811e2dafb3b
