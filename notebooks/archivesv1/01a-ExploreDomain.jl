### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 8e761852-a62c-41a8-856f-7d3b3fba9c20
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 682367f2-9bef-4512-a741-20d5c43c6311
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

# ╔═╡ 8bd812fc-5c56-11ec-3514-f7b178ccb922
md"
# 01a - Exploring the OTREC Region

In this notebook, we explore the OTREC Region, specifically the topography/orography in relation to the station locations, for future reference in case we need to seperate by different land types.
"

# ╔═╡ f8dd4814-bf9d-4dbb-9d45-a0413e4710d6
md"
### A. Topography for the OTREC Region
"

# ╔═╡ b4b236b7-233c-4c6a-b3ad-4c93fcb74367
begin
	ds  = NCDataset(datadir("ETOPO1_Ice_cell.grd"))
	oro = ds["z"][:]; lon,lat = NASAPrecipitation.gpmlonlat()
	lat = reverse(lat)
	close(ds)
	md"Loading topographic data from ETOPO1 grid ..."
end

# ╔═╡ 5ef3c8c5-ef45-4695-8423-4328e43df484
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 870dbe5f-4437-40de-ba6e-581249a7bf51
geo = GeoRegion("OTREC")

# ╔═╡ 009dd0ce-a79c-4137-92b0-4c032a7b758b
begin
	soro = dropdims(mean(reshape(oro,6,3600,6,1800),dims=(1,3)),dims=(1,3))
	md"ETOPO1 data is given in units of 1-arcminute.  Averaging to 0.1º ..."
end

# ╔═╡ 2c3a1bb6-5d07-4f60-9853-3f6f96b379ec
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
	roro = zeros(nlon,nlat)
	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end
	for glat in 1 : nlat, glon in 1 : nlon
		roro[glon,glat] = soro[ilon[glon],ilat[glat]] * mask[glon,glat]
	end
	md"Extracting information for region ..."
end

# ╔═╡ 12a4b74a-ccd6-4e9b-992a-13c0fcded188
begin
	fetopo = datadir("ETOPO1_Ice_cell-OTREC.nc")
	if isfile(fetopo); rm(fetopo,force=true) end
	
	eds = NCDataset(fetopo,"c")
	eds.dim["longitude"] = nlon
	eds.dim["latitude"]  = nlat
	
	nclon = defVar(eds,"longitude", Float32, ("longitude",))
	nclat = defVar(eds,"latitude",  Float32, ("latitude",))
	ncoro = defVar(eds,"z", Float32, ("longitude", "latitude"))
	
	nclon[:] = ggrd.glon
	nclat[:] = ggrd.glat
	ncoro[:] = roro
	
	close(eds)
	md"Creating ETOPO1 file for the OTREC region ..."
end

# ╔═╡ 786f11e6-bbe5-427d-9e7d-a55e1c94cbe9
md"
### B. Station Data
"

# ╔═╡ 530286d6-f2c0-4de7-844b-9f943696d311
begin
	infoall = stninfoall(); nstn = size(infoall,1)
	infody  = stninfody();  ndy  = size(infody,1)
	infomo  = stninfomo();  nmo  = size(infomo,1)
	md"Loading station location information ..."
end

# ╔═╡ 7eadf4a3-1442-4b61-ae40-6642b11fac6c
md"
### C. Plotting the Station Data against Topography
"

# ╔═╡ 6b85fd09-7025-4eb2-9a7f-c5bdfad593f4
begin
	clr_dy = pplt.Colors("Gray",size(infody,1)+4)
	clr_mo = pplt.Colors("Phase",size(infomo,1))
	md"Loading colors to distinguish between different stations (and separate between daily and monthly stations) in the figure ..."
end

# ╔═╡ 7b1a4e9e-45c1-40d9-9fa8-1825368f63dc
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=1,axwidth=3)
	
	c = axs[1].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=-5:5,extend="both"
	)

	for ii = 1 : nmo
		axs[1].scatter(
			infomo[ii,2],infomo[ii,3],label=infomo[ii,1],c=clr_mo[ii],
			legend="b",legend_kw=Dict("ncol"=>3,"frame"=>false)
		)
	end
	
	for ii = 1 : ndy
		axs[1].scatter(
			infody[ii,2],infody[ii,3],label=infody[ii,1],c=clr_dy[ii+2],
			legend="b",legend_kw=Dict("ncol"=>3,"frame"=>false)
		)
	end
	
	axs[1].plot(x,y,lw=0.5,c="k")
	axs[1].format(xlim=(-85,-70),ylim=(0,15),xlocator=-85:2.5:-70,ylocator=0:2.5:15)
	axs[1].colorbar(c,loc="r")
	
	fig.savefig(plotsdir("01a-stations.png"),transparent=false,dpi=200)
	load(plotsdir("01a-stations.png"))
end

# ╔═╡ 3b70530b-d518-4de8-ab33-29c050452ca6
md"
### D. Defining the GeoRegions ...
"

# ╔═╡ 09e38bfc-866c-4005-bed2-436c4f22d97d
md"
It is made manifestly obvious that there would be some very different climatologies at play in this region, not in the least within Columbia, where the meteorological stations span across two mountain ranges.  Here, we define the following subregions to OTREC:
* CIS_CLB - The main domain where the Columbia weather stations are located
* CIS_CST - Coastal Columbian *subdomain* in the west, high precipitation
* CIS_MNT - Mountaineous Columbian *subdomain*, includes the valleys (can separate between different stations via height)
* CIS_INL - Inland Columbia *subdomain* to the east of the mountain ranges
"

# ╔═╡ 609a239a-7510-4e20-82bb-f7a252754950
begin
	if isGeoRegion("CIS_CLB",throw=false)
		removeGeoRegion("CIS_CLB")
	end
	CIS = RectRegion(
		"CIS_CLB","OTREC","Columbia Isotope Project Domain",
		[11,0,-72,-77.5]
	)
end

# ╔═╡ 2257ff8a-57e5-46e6-90fb-bbc99c8dedfd
begin
	if isGeoRegion("CIS_CST",throw=false)
		removeGeoRegion("CIS_CST")
	end
	CST = PolyRegion(
		"CIS_CST","CIS_CLB","Coastal Columbia",
		[-77.5,-77.5,-77,-76,-75.5,-75,-73.5,-76.3,-76.3,-77.2,-77.5],
		[2.3,6.5,8.5,9.5,10.7,11,9,7.5,5,2.7,2.3]
	)
end

# ╔═╡ 42c5a013-a21b-415e-b6b9-a37a661fcbf2
begin
	if isGeoRegion("CIS_MNT",throw=false)
		removeGeoRegion("CIS_MNT")
	end
	MNT = PolyRegion(
		"CIS_MNT","CIS_CLB","Mountaineous Columbia",
		[-72,-75,-77.5,-77.5,-77.2,-76.3,-76.3,-73.5,-72,-72],
		[6.2,2.2,0.2,2.3,2.7,5,7.5,9,7.5,6.2]
	)
end

# ╔═╡ ca22ab85-5423-47c4-adf7-367024844170
begin
	if isGeoRegion("CIS_INL",throw=false)
		removeGeoRegion("CIS_INL")
	end
	INL = PolyRegion(
		"CIS_INL","CIS_CLB","Inland Columbia",
		[-77.5,-72,-72,-75,-77.5,-77.5],
		[0,0,6.2,2.2,0.2,0]
	)
end

# ╔═╡ 68df62c9-4eed-4829-9ea5-199937435615
begin
	blon_CIS,blat_CIS = coordGeoRegion(CIS)
	blon_CST,blat_CST,slon_CST,slat_CST = coordGeoRegion(CST)
	blon_MNT,blat_MNT,slon_MNT,slat_MNT = coordGeoRegion(MNT)
	blon_INL,blat_INL,slon_INL,slat_INL = coordGeoRegion(INL)
	md"Loading geographic shapes of the GeoRegions ..."
end

# ╔═╡ 2f4194c8-c87c-4487-8b20-43e7148f353d
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=7.5/13,axwidth=2)
	
	c1 = a1[1].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=vcat(-5:0.5:-1,0,1:0.5:5),extend="both"
	)
	
	a1[1].plot(x,y,lw=0.5,c="k")
	
	a1[1].plot(slon_CST,slat_CST)
	a1[1].plot(slon_MNT,slat_MNT)
	a1[1].plot(slon_INL,slat_INL)
	a1[1].plot(blon_CIS,blat_CIS,c="k",linestyle="--")
	a1[1].scatter(infoall[:,2],infoall[:,3],c="k",s=10,zorder=4)
	
	a1[1].fill(slon_CST,slat_CST,alpha=0.4)
	a1[1].fill(slon_MNT,slat_MNT,alpha=0.4)
	a1[1].fill(slon_INL,slat_INL,alpha=0.4)
	
	a1[1].text(-76,8.7,"CST")
	a1[1].text(-75.5,6.5,"MNT")
	a1[1].text(-74,1,"INL")
	
	a1_2 = f1.add_axes([0.56, 0.811, 0.18, 0.4/3])
	a1_2.contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=-5:5,extend="both"
	)
	a1_2.plot([-77.5,-72,-72,-77.5,-77.5],[0,0,11,11,0],c="r")
	a1_2.scatter(infoall[1,2],infoall[1,3],c="k",s=10,zorder=4)
	
	a1[1].format(
		xlim=(-78.5,-71),ylim=(-1,12),
		xlocator=-85:2.5:-70,ylocator=-2.5:2.5:15,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$"
	)
	a1[1].colorbar(c1,loc="r",locator=vcat(-5:5),label="Topography / km")
	
	f1.savefig(plotsdir("01a-GeoRegions.png"),transparent=false,dpi=200)
	load(plotsdir("01a-GeoRegions.png"))
end

# ╔═╡ 919973d4-3673-4989-a98d-5087b117567d
md"
### E. Additional subGeoRegions

From results shown in later data, we see that we can and should subdivide the **Coastal** and **Mountaineous** regions into different parts.
* CIS_PAC - the Pacific Coastal region
* CIS_CRB - the Caribbean Coastal region
* CIS_WAN - the Western Andes mountain range
* CIS_EAN - the Eastern Andes mountain range
* CIS_VAL - the valley between the western and eastern Andes mountain ranges

We define them below:
"

# ╔═╡ 3fe4b855-279d-4afa-bb34-b024f22e5977
begin
	if isGeoRegion("CIS_PAC",throw=false)
		removeGeoRegion("CIS_PAC")
	end
	PAC = PolyRegion(
		"CIS_PAC","CIS_CLB","Coastal Pacific Columbia",
		[-77.5,-77.5,-76.3,-76.3,-77.2,-77.5],
		[2.3,6.5,6.5,5,2.7,2.3]
	)
end

# ╔═╡ b7083cf1-51a1-472a-9af5-ffc1b7e9b101
begin
	if isGeoRegion("CIS_CRB",throw=false)
		removeGeoRegion("CIS_CRB")
	end
	CRB = PolyRegion(
		"CIS_CRB","CIS_CLB","Coastal Caribbean Columbia",
		[-77,-76,-75.5,-75.05,-74.6,-75.2,-76.8,-77],
		[8.5,9.5,10.7,11,11,9,7.5,8.5]
	)
end

# ╔═╡ 50ebc420-ed27-4239-9aa8-5aa6e549ff2d
begin
	if isGeoRegion("CIS_WAN",throw=false)
		removeGeoRegion("CIS_WAN")
	end
	WAN = PolyRegion(
		"CIS_WAN","CIS_CLB","West Andes Columbia",
		[-77.5,-77.5,-77.2,-76.3,-76.3,-75.5,-74.8,-75,-76,-77.5],
		[0.2,2.3,2.7,5,7.5,7.5,7,4.5,1.4,0.2]
	)
end

# ╔═╡ 2633d628-1434-4727-a0ca-16c440a6eb20
begin
	if isGeoRegion("CIS_EAN",throw=false)
		removeGeoRegion("CIS_EAN")
	end
	EAN = PolyRegion(
		"CIS_EAN","CIS_CLB","East Andes Columbia",
		[-73.5,-72,-72,-75,-76,-75.8,-74.5,-74.5,-73.2,-73.5],
		[9,7.5,6.2,2.2,1.4,2.02,4.5,5.5,7.2,9]
	)
end

# ╔═╡ 8a6a805a-a80a-4f60-a574-b0c9b12d75ff
begin
	if isGeoRegion("CIS_VAL",throw=false)
		removeGeoRegion("CIS_VAL")
	end
	VAL = PolyRegion(
		"CIS_VAL","CIS_CLB","Andes Valley Columbia",
		[-75.8,-74.5,-74.5,-73.2,-74.8,-75,-75.8],
		[2.02,4.5,5.5,7.2,7,4.5,2.02]
	)
end

# ╔═╡ 19d9fd62-2cae-43ec-9d2b-d0be2c9b987a
begin
	blon_PAC,blat_PAC,slon_PAC,slat_PAC = coordGeoRegion(PAC)
	blon_CRB,blat_CRB,slon_CRB,slat_CRB = coordGeoRegion(CRB)
	blon_WAN,blat_WAN,slon_WAN,slat_WAN = coordGeoRegion(WAN)
	blon_EAN,blat_EAN,slon_EAN,slat_EAN = coordGeoRegion(EAN)
	blon_VAL,blat_VAL,slon_VAL,slat_VAL = coordGeoRegion(VAL)
	md"Loading geographic shapes of the subdivisions to the original GeoRegions ..."
end

# ╔═╡ f2c2da70-340c-4b2f-8d3c-35ba09fcc13c
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=7.5/13,axwidth=2)
	
	c2 = a2[1].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=vcat(-5:0.5:-1,0,1:0.5:5),extend="both"
	)
	
	a2[1].plot(x,y,lw=0.5,c="k")
	
	a2[1].plot(slon_PAC,slat_PAC)
	a2[1].plot(slon_CRB,slat_CRB)
	a2[1].plot(slon_WAN,slat_WAN)
	a2[1].plot(slon_EAN,slat_EAN)
	a2[1].plot(slon_VAL,slat_VAL)
	a2[1].plot(blon_CIS,blat_CIS,c="k",linestyle="--")
	a2[1].scatter(infoall[:,2],infoall[:,3],c="k",s=10,zorder=4)
	
	a2[1].fill(slon_PAC,slat_PAC,alpha=0.4)
	a2[1].fill(slon_CRB,slat_CRB,alpha=0.4)
	a2[1].fill(slon_WAN,slat_WAN,alpha=0.4)
	a2[1].fill(slon_EAN,slat_EAN,alpha=0.4)
	a2[1].fill(slon_VAL,slat_VAL,alpha=0.4)
	
	a2[1].text(-77.4,5,"PAC")
	a2[1].text(-76,8.7,"CRB")
	a2[1].text(-76.5,3.7,"WAN")
	a2[1].text(-74.7,6.5,"VAL")
	a2[1].text(-74,5.5,"EAN")
	
	a2_2 = f2.add_axes([0.56, 0.811, 0.18, 0.4/3])
	a2_2.contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=-5:5,extend="both"
	)
	a2_2.plot([-77.5,-72,-72,-77.5,-77.5],[0,0,11,11,0],c="r")
	a2_2.scatter(infoall[1,2],infoall[1,3],c="k",s=10,zorder=4)
	
	a2[1].format(
		xlim=(-78.5,-71),ylim=(-1,12),
		xlocator=-85:2.5:-70,ylocator=-2.5:2.5:15,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$"
	)
	a2[1].colorbar(c2,loc="r",locator=vcat(-5:5),label="Topography / km")
	
	f2.savefig(plotsdir("01a-subGeoRegions.png"),transparent=false,dpi=200)
	load(plotsdir("01a-subGeoRegions.png"))
end

# ╔═╡ dd844840-ce0c-4be2-b611-a2e56b152090
md"And all your domains and stations combined ..."

# ╔═╡ dd0eee10-9507-488d-ac18-171d3f04eac6
begin
	pplt.close(); f3,a3 = pplt.subplots(ncols=3,aspect=1,axwidth=2)

	a3[1].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=vcat(-5:0.5:-1,0,1:0.5:5),extend="both"
	)
	a3[1].plot(x,y,lw=0.5,c="k")

	a3[1].plot(slon_CST,slat_CST)
	a3[1].plot(slon_MNT,slat_MNT)
	a3[1].plot(slon_INL,slat_INL)
	a3[1].plot(blon_CIS,blat_CIS,c="k",linestyle="--")
	
	a3[1].fill(slon_CST,slat_CST,alpha=0.4,c="k")
	a3[1].fill(slon_MNT,slat_MNT,alpha=0.4,c="k")
	a3[1].fill(slon_INL,slat_INL,alpha=0.4,c="k")

	a3[1].scatter(
		infomo[:,2],infomo[:,3],label="Monthly",c="r",
		legend="l",legend_kw=Dict("ncol"=>1,"frame"=>false),zorder=4
	)
	a3[1].scatter(infody[:,2],infody[:,3],label="Daily",c="k",legend="l",zorder=5)
	
	a3[2].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=vcat(-5:0.5:-1,0,1:0.5:5),extend="both"
	)
	a3[2].plot(x,y,lw=0.5,c="k")
	
	a3[2].plot(slon_PAC,slat_PAC)
	a3[2].plot(slon_CRB,slat_CRB)
	a3[2].plot(slon_WAN,slat_WAN)
	a3[2].plot(slon_EAN,slat_EAN)
	a3[2].plot(slon_VAL,slat_VAL)
	a3[2].plot(blon_CIS,blat_CIS,c="k",linestyle="--")
	
	a3[2].fill(slon_PAC,slat_PAC,alpha=0.4,c="k")
	a3[2].fill(slon_CRB,slat_CRB,alpha=0.4,c="k")
	a3[2].fill(slon_WAN,slat_WAN,alpha=0.4,c="k")
	a3[2].fill(slon_EAN,slat_EAN,alpha=0.4,c="k")
	a3[2].fill(slon_VAL,slat_VAL,alpha=0.4,c="k")
	a3[2].scatter(infomo[:,2],infomo[:,3],c="r",zorder=5)
	a3[2].scatter(infody[:,2],infody[:,3],c="k",zorder=5)

	for ii = 1 : nmo
		a3[2].scatter(
			NaN,NaN,label=infomo[ii,1],c=clr_mo[ii],
			legend="b",legend_kw=Dict("ncol"=>4,"frame"=>false)
		)
	end
	
	for ii = 1 : ndy
		a3[2].scatter(
			NaN,NaN,label=infody[ii,1],c=clr_dy[ii+2],
			legend="b",legend_kw=Dict("ncol"=>4,"frame"=>false)
		)
	end
	
	c3 = a3[3].contourf(
		ggrd.glon,ggrd.glat,roro'/1000,
		cmap="delta",levels=vcat(-5:0.5:-1,0,1:0.5:5),extend="both"
	)
	a3[3].plot(x,y,lw=0.5,c="k")
	
	a3[3].plot(slon_PAC,slat_PAC)
	a3[3].plot(slon_CRB,slat_CRB)
	a3[3].plot(slon_WAN,slat_WAN)
	a3[3].plot(slon_EAN,slat_EAN)
	a3[3].plot(slon_VAL,slat_VAL)
	a3[3].plot(blon_CIS,blat_CIS,c="k",linestyle="--")
	
	a3[3].fill(slon_PAC,slat_PAC,alpha=0.4,c="k")
	a3[3].fill(slon_CRB,slat_CRB,alpha=0.4,c="k")
	a3[3].fill(slon_WAN,slat_WAN,alpha=0.4,c="k")
	a3[3].fill(slon_EAN,slat_EAN,alpha=0.4,c="k")
	a3[3].fill(slon_VAL,slat_VAL,alpha=0.4,c="k")

	for ii = 1 : nmo; a3[3].scatter(infomo[ii,2],infomo[ii,3],c=clr_mo[ii],zorder=4) end
	for ii = 1 : ndy; a3[3].scatter(infody[ii,2],infody[ii,3],c=clr_dy[ii+2],zorder=4) end

	for ax in a3
		ax.format(
			xlim=(-85,-70),ylim=(0,15),xlocator=-85:2.5:-70,ylocator=0:2.5:15,
			xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
			suptitle="Map of Stations in Columbia"
		)
	end
	
	a3[3].colorbar(c3,loc="r",label="Topography / km")
	f3.savefig(plotsdir("01a-ExploreDomain.png"),transparent=false,dpi=300)
	load(plotsdir("01a-ExploreDomain.png"))
end

# ╔═╡ 4efe19cf-69ea-4eb3-a98b-f00c16db8207
md"
In this, we therefore see that the stations where monthly and daily data are available are likely in climatologically different regions.  All the stations where data is daily in temporal in resolution are concentrated along the Pacific Coastline, while the stations that are monthly in temporal resolution are inland, with the exceptions of the San Andres and Cartagena, the latter of which is along the Caribbean Coastline.

Therefore, it is of interest to see if the climatologies of the Pacific and Caribbean Coastlines are comparable to those of the inland domains, and over San Andres.  Does San Andres have a different climatology compared to all the stations on mainland Columbia?
"

# ╔═╡ 12a8c9a8-5697-4aa3-a228-c90a5307dbe6
begin
	if isGeoRegion("CIS_SAN",throw=false)
		removeGeoRegion("CIS_SAN")
	end
	SAN = RectRegion(
		"CIS_SAN","OTREC","San Andres Columbia",
		[13.1,12.1,-81.2,-82.2],
	)
end

# ╔═╡ Cell order:
# ╟─8bd812fc-5c56-11ec-3514-f7b178ccb922
# ╟─8e761852-a62c-41a8-856f-7d3b3fba9c20
# ╟─682367f2-9bef-4512-a741-20d5c43c6311
# ╟─f8dd4814-bf9d-4dbb-9d45-a0413e4710d6
# ╟─b4b236b7-233c-4c6a-b3ad-4c93fcb74367
# ╟─5ef3c8c5-ef45-4695-8423-4328e43df484
# ╟─870dbe5f-4437-40de-ba6e-581249a7bf51
# ╟─009dd0ce-a79c-4137-92b0-4c032a7b758b
# ╟─2c3a1bb6-5d07-4f60-9853-3f6f96b379ec
# ╟─12a4b74a-ccd6-4e9b-992a-13c0fcded188
# ╟─786f11e6-bbe5-427d-9e7d-a55e1c94cbe9
# ╟─530286d6-f2c0-4de7-844b-9f943696d311
# ╟─7eadf4a3-1442-4b61-ae40-6642b11fac6c
# ╟─6b85fd09-7025-4eb2-9a7f-c5bdfad593f4
# ╟─7b1a4e9e-45c1-40d9-9fa8-1825368f63dc
# ╟─3b70530b-d518-4de8-ab33-29c050452ca6
# ╟─09e38bfc-866c-4005-bed2-436c4f22d97d
# ╟─609a239a-7510-4e20-82bb-f7a252754950
# ╟─2257ff8a-57e5-46e6-90fb-bbc99c8dedfd
# ╟─42c5a013-a21b-415e-b6b9-a37a661fcbf2
# ╟─ca22ab85-5423-47c4-adf7-367024844170
# ╟─68df62c9-4eed-4829-9ea5-199937435615
# ╟─2f4194c8-c87c-4487-8b20-43e7148f353d
# ╟─919973d4-3673-4989-a98d-5087b117567d
# ╟─3fe4b855-279d-4afa-bb34-b024f22e5977
# ╟─b7083cf1-51a1-472a-9af5-ffc1b7e9b101
# ╟─50ebc420-ed27-4239-9aa8-5aa6e549ff2d
# ╟─2633d628-1434-4727-a0ca-16c440a6eb20
# ╟─8a6a805a-a80a-4f60-a574-b0c9b12d75ff
# ╟─19d9fd62-2cae-43ec-9d2b-d0be2c9b987a
# ╟─f2c2da70-340c-4b2f-8d3c-35ba09fcc13c
# ╟─dd844840-ce0c-4be2-b611-a2e56b152090
# ╟─dd0eee10-9507-488d-ac18-171d3f04eac6
# ╟─4efe19cf-69ea-4eb3-a98b-f00c16db8207
# ╟─12a8c9a8-5697-4aa3-a228-c90a5307dbe6
