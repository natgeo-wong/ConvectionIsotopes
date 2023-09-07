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
	using DataFrames
	using DelimitedFiles
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	using XLSX
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ 8bd812fc-5c56-11ec-3514-f7b178ccb922
md"
# 01c - Exploring the Precipitation Climatology of Columbia in ERA5

In this notebook, we explore the precipitation climatology of the regions defined in the `01a-ExploreDomain` notebooks made using GPM from 2001-2020.
"

# ╔═╡ f8dd4814-bf9d-4dbb-9d45-a0413e4710d6
md"
### A. Topography for the OTREC Region
"

# ╔═╡ b4b236b7-233c-4c6a-b3ad-4c93fcb74367
begin
	ds  = NCDataset(datadir("ETOPO1_Ice_cell-OTREC.nc"))
	oro = ds["z"][:]
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	close(ds)
	md"Loading topographic data for the OTREC region ..."
end

# ╔═╡ 5ef3c8c5-ef45-4695-8423-4328e43df484
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 3b70530b-d518-4de8-ab33-29c050452ca6
md"
### B. Defining the GeoRegions ...
"

# ╔═╡ 09e38bfc-866c-4005-bed2-436c4f22d97d
md"
It is made manifestly obvious that there would be some very different climatologies at play in this region, not in the least within Columbia, where the meteorological stations span across two mountain ranges.  Here, we use the following subregions to OTREC defined in the notebook `01a-ExploreDomain`:
* CIS_CLB - The main domain where the Columbia weather stations are located
* CIS_CST - Coastal Columbian *subdomain* in the west, high precipitation
* CIS_MNT - Mountaineous Columbian *subdomain*, includes the valleys (can separate between different stations via height)
* CIS_INL - Inland Columbia *subdomain* to the east of the mountain ranges
"

# ╔═╡ 870dbe5f-4437-40de-ba6e-581249a7bf51
geo_OTREC = GeoRegion("OTREC")

# ╔═╡ 844cecf6-c06b-4178-9d2f-1c3d25a27252
begin
	geo_CIS = GeoRegion("CIS_CLB")
	geo_CST = GeoRegion("CIS_CST")
	geo_MNT = GeoRegion("CIS_MNT")
	geo_INL = GeoRegion("CIS_INL")
	md"Loading GeoRegions in the Domain ..."
end

# ╔═╡ 68df62c9-4eed-4829-9ea5-199937435615
begin
	blon_CIS,blat_CIS = coordGeoRegion(geo_CIS)
	blon_CST,blat_CST,slon_CST,slat_CST = coordGeoRegion(geo_CST)
	blon_MNT,blat_MNT,slon_MNT,slat_MNT = coordGeoRegion(geo_MNT)
	blon_INL,blat_INL,slon_INL,slat_INL = coordGeoRegion(geo_INL)
	md"Loading geographic shapes of the GeoRegions ..."
end

# ╔═╡ 409d48fb-8a91-497a-9d52-3cb98797dacc
md"
### C. Exploring Monthly Precipitation Data
"

# ╔═╡ c3a4c472-88cb-4573-82a8-30cc3ed8ad97
begin
	pds  = NCDataset(datadir("reanalysis","era5m-CIS_CLBx0.25-tp.nc"))
	plon = pds["longitude"][:]; nplon = length(plon)
	plat = pds["latitude"][:];  nplat = length(plat)
	prcp = reshape(pds["tp"][:],nplon,nplat,12,:) * 1000
	prcp_yr = dropdims(mean(prcp,dims=(3,4)),dims=(3,4))
	prcp_mo = dropdims(mean(prcp,dims=4),dims=4)
	close(pds)
	md"Loading ERA5 Precipitation data for all years"
end

# ╔═╡ 7902b269-124e-47c5-9745-be12c0b5b431
begin
	arr = [[1,1,1,2,3,4,5],[1,1,1,6,7,8,9],[1,1,1,10,11,12,13]]
	
	pplt.close(); f1,a1 = pplt.subplots(
		arr,aspect=7.5/13,axwidth=2,
		hspace=[0.1,0.1],wspace=[0,0,0.2,0.1,0.1,0.1]
	)
	
	for iax = 1 : 13
		if !isone(iax)
			c = a1[iax].contourf(
				plon,plat,prcp_mo[:,:,iax-1]',
				levels=5:2.5:35,extend="both",cmap="Blues"
			)
			if iax == 13
				f1.colorbar(c,loc="r",label=L"mm day$^{-1}$")
			end
			a1[iax].plot(slon_CST,slat_CST,lw=1)
			a1[iax].plot(slon_MNT,slat_MNT,lw=1)
			a1[iax].plot(slon_INL,slat_INL,lw=1)
			a1[iax].plot(blon_CIS,blat_CIS,lw=1,c="k",linestyle="--")
			a1[iax].contour(
				lon,lat,oro'/1000,
				lw=0.5,cmap="brown1",levels=1:5,extend="both"
			)
		else
			c = a1[iax].contourf(
				plon,plat,prcp_yr',
				levels=5:25,extend="both",cmap="Blues"
			)
			a1[iax].colorbar(c,loc="l",label=L"mm day$^{-1}$")
			a1[iax].plot(slon_CST,slat_CST,lw=1.5)
			a1[iax].plot(slon_MNT,slat_MNT,lw=1.5)
			a1[iax].plot(slon_INL,slat_INL,lw=1.5)
			a1[iax].plot(blon_CIS,blat_CIS,lw=1.5,c="k",linestyle="--")
			a1[iax].contour(lon,lat,oro'/1000,cmap="brown1",levels=1:5,extend="both")
		end
		
		
	end
	
	for ax in a1
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-78.5,-71),ylim=(-1,12))
	end
	
	a1[1].format(xlocator=-85:2.5:-70,ylocator=-2.5:2.5:15)
	for iax = 2 : 13
		a1[iax].format(xlocator=[],ylocator=[],lltitle="$(iax-1)")
	end
	
	a1[1].format(ltitle="(a) Yearly Precipitation Rate")
	a1[2].format(ltitle="(b) Monthly Precipitation Rate")
	
	f1.savefig(plotsdir("01c-era5tp.png"),transparent=false,dpi=300)
	load(plotsdir("01c-era5tp.png"))
end

# ╔═╡ 3ef9288b-c02a-4918-ac67-76036be08ceb
md"Now, we proceed by extracting and plotting the monthly-averages for each of the specified GeoRegions within our domain ..."

# ╔═╡ 2c75ada3-dff0-45f1-85a2-d8d08327ec14
begin
	ggrd_CIS = RegionGrid(geo_CIS,plon,plat)
	prcp_yr_CIS_tmp = zeros(length(ggrd_CIS.ilon),length(ggrd_CIS.ilat))
	prcp_mo_CIS_tmp = zeros(length(ggrd_CIS.ilon),length(ggrd_CIS.ilat),12)
	prcp_mo_CIS     = zeros(12)
	if typeof(ggrd_CIS) <: PolyGrid
		  mask_CIS = ggrd_CIS.mask
	else; mask_CIS = ones(length(ggrd_CIS.ilon),length(ggrd_CIS.ilat))
	end
	for glat in 1 : length(ggrd_CIS.ilat), glon in 1 : length(ggrd_CIS.ilon)
		prcp_yr_CIS_tmp[glon,glat] = prcp_yr[
			ggrd_CIS.ilon[glon],ggrd_CIS.ilat[glat]
		] * mask_CIS[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_CIS.ilat), glon in 1 : length(ggrd_CIS.ilon)
			prcp_mo_CIS_tmp[glon,glat,imo] = prcp_mo[
				ggrd_CIS.ilon[glon],ggrd_CIS.ilat[glat],imo
			] * mask_CIS[glon,glat]
		end
	end
	prcp_yr_CIS = mean(prcp_yr_CIS_tmp[.!isnan.(prcp_yr_CIS_tmp)])
	for imo = 1 : 12
		prcp_mo_CIS_ii   = @view prcp_mo_CIS_tmp[:,:,imo]
		prcp_mo_CIS[imo] = mean(prcp_mo_CIS_ii[.!isnan.(prcp_mo_CIS_ii)])
	end
	md"Extracting information for the **$(geo_CIS.name)** region ..."
end

# ╔═╡ d15b987d-fed8-49ae-a077-b1715da02d14
begin
	ggrd_CST = RegionGrid(geo_CST,plon,plat)
	prcp_yr_CST_tmp = zeros(length(ggrd_CST.ilon),length(ggrd_CST.ilat))
	prcp_mo_CST_tmp = zeros(length(ggrd_CST.ilon),length(ggrd_CST.ilat),12)
	prcp_mo_CST     = zeros(12)
	if typeof(ggrd_CST) <: PolyGrid
		  mask_CST = ggrd_CST.mask
	else; mask_CST = ones(length(ggrd_CST.ilon),length(ggrd_CST.ilat))
	end
	for glat in 1 : length(ggrd_CST.ilat), glon in 1 : length(ggrd_CST.ilon)
		prcp_yr_CST_tmp[glon,glat] = prcp_yr[
			ggrd_CST.ilon[glon],ggrd_CST.ilat[glat]
		] * mask_CST[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_CST.ilat), glon in 1 : length(ggrd_CST.ilon)
			prcp_mo_CST_tmp[glon,glat,imo] = prcp_mo[
				ggrd_CST.ilon[glon],ggrd_CST.ilat[glat],imo
			] * mask_CST[glon,glat]
		end
	end
	prcp_yr_CST = mean(prcp_yr_CST_tmp[.!isnan.(prcp_yr_CST_tmp)])
	for imo = 1 : 12
		prcp_mo_CST_ii   = @view prcp_mo_CST_tmp[:,:,imo]
		prcp_mo_CST[imo] = mean(prcp_mo_CST_ii[.!isnan.(prcp_mo_CST_ii)])
	end
	md"Extracting information for the **$(geo_CST.name)** region ..."
end

# ╔═╡ 80a51754-af0a-4756-a251-922c9e5d1501
begin
	ggrd_MNT = RegionGrid(geo_MNT,plon,plat)
	prcp_yr_MNT_tmp = zeros(length(ggrd_MNT.ilon),length(ggrd_MNT.ilat))
	prcp_mo_MNT_tmp = zeros(length(ggrd_MNT.ilon),length(ggrd_MNT.ilat),12)
	prcp_mo_MNT     = zeros(12)
	if typeof(ggrd_MNT) <: PolyGrid
		  mask_MNT = ggrd_MNT.mask
	else; mask_MNT = ones(length(ggrd_MNT.ilon),length(ggrd_MNT.ilat))
	end
	for glat in 1 : length(ggrd_MNT.ilat), glon in 1 : length(ggrd_MNT.ilon)
		prcp_yr_MNT_tmp[glon,glat] = prcp_yr[
			ggrd_MNT.ilon[glon],ggrd_MNT.ilat[glat]
		] * mask_MNT[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_MNT.ilat), glon in 1 : length(ggrd_MNT.ilon)
			prcp_mo_MNT_tmp[glon,glat,imo] = prcp_mo[
				ggrd_MNT.ilon[glon],ggrd_MNT.ilat[glat],imo
			] * mask_MNT[glon,glat]
		end
	end
	prcp_yr_MNT = mean(prcp_yr_MNT_tmp[.!isnan.(prcp_yr_MNT_tmp)])
	for imo = 1 : 12
		prcp_mo_MNT_ii   = @view prcp_mo_MNT_tmp[:,:,imo]
		prcp_mo_MNT[imo] = mean(prcp_mo_MNT_ii[.!isnan.(prcp_mo_MNT_ii)])
	end
	md"Extracting information for the **$(geo_MNT.name)** region ..."
end

# ╔═╡ 20532d37-83eb-4260-8658-d1261b8a23d4
begin
	ggrd_INL = RegionGrid(geo_INL,plon,plat)
	prcp_yr_INL_tmp = zeros(length(ggrd_INL.ilon),length(ggrd_INL.ilat))
	prcp_mo_INL_tmp = zeros(length(ggrd_INL.ilon),length(ggrd_INL.ilat),12)
	prcp_mo_INL     = zeros(12)
	if typeof(ggrd_INL) <: PolyGrid
		  mask_INL = ggrd_INL.mask
	else; mask_INL = ones(length(ggrd_INL.ilon),length(ggrd_INL.ilat))
	end
	for glat in 1 : length(ggrd_INL.ilat), glon in 1 : length(ggrd_INL.ilon)
		prcp_yr_INL_tmp[glon,glat] = prcp_yr[
			ggrd_INL.ilon[glon],ggrd_INL.ilat[glat]
		] * mask_INL[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_INL.ilat), glon in 1 : length(ggrd_INL.ilon)
			prcp_mo_INL_tmp[glon,glat,imo] = prcp_mo[
				ggrd_INL.ilon[glon],ggrd_INL.ilat[glat],imo
			] * mask_INL[glon,glat]
		end
	end
	prcp_yr_INL = mean(prcp_yr_INL_tmp[.!isnan.(prcp_yr_INL_tmp)])
	for imo = 1 : 12
		prcp_mo_INL_ii   = @view prcp_mo_INL_tmp[:,:,imo]
		prcp_mo_INL[imo] = mean(prcp_mo_INL_ii[.!isnan.(prcp_mo_INL_ii)])
	end
	md"Extracting information for the **$(geo_INL.name)** region ..."
end

# ╔═╡ b016350b-adc7-4b24-b7c5-d54c61fed723
begin
	pplt.close()
	f2,a2 = pplt.subplots([1,1,1,1,1,1,1,1,1,1,2],aspect=1.5,wspace=0.1,sharex=0)
	
	a2[1].plot(0:13,vcat(prcp_mo_CIS[end],prcp_mo_CIS,prcp_mo_CIS[1]),c="k")
	a2[1].plot(0:13,vcat(prcp_mo_CST[end],prcp_mo_CST,prcp_mo_CST[1]))
	a2[1].plot(0:13,vcat(prcp_mo_MNT[end],prcp_mo_MNT,prcp_mo_MNT[1]))
	a2[1].plot(0:13,vcat(prcp_mo_INL[end],prcp_mo_INL,prcp_mo_INL[1]))
	a2[1].format(xlabel="Month of Year",xlim=(0.5,12.5),ylim=(0,25),xlocator=1:12)
	
	a2[2].scatter(0.5,prcp_yr_CIS,s=10,label="Entire Domain",legend="r",c="k")
	a2[2].scatter(0.5,prcp_yr_CST,s=10,label="Coastal Columbia",legend="r")
	a2[2].scatter(0.5,prcp_yr_MNT,s=10,label="Mountaineous Columbia",legend="r")
	a2[2].scatter(0.5,prcp_yr_INL,s=10,label="Inland Columbia",legend="r",legend_kw=Dict("ncols"=>1,"frame"=>false))
	a2[2].format(ylabel=L"Rainfall Rate / mm day$^{-1}$",xlocator=[])
	
	f2.savefig(plotsdir("01c-era5tp-ts.png"),transparent=false,dpi=300)
	load(plotsdir("01c-era5tp-ts.png"))
end

# ╔═╡ 70fef982-f550-4ac5-b1be-21e7e06380e7
md"But what about our subGeoRegions?"

# ╔═╡ 79be324c-701c-46c3-b7e6-33e2ac5d75c8
begin
	geo_PAC = GeoRegion("CIS_PAC")
	geo_CRB = GeoRegion("CIS_CRB")
	geo_EAN = GeoRegion("CIS_EAN")
	geo_VAL = GeoRegion("CIS_VAL")
	geo_WAN = GeoRegion("CIS_WAN")
	md"Loading subGeoRegion information ..."
end

# ╔═╡ 831c1668-7228-479d-ae1e-cdd176e8010a
begin
	ggrd_PAC = RegionGrid(geo_PAC,plon,plat)
	prcp_yr_PAC_tmp = zeros(length(ggrd_PAC.ilon),length(ggrd_PAC.ilat))
	prcp_mo_PAC_tmp = zeros(length(ggrd_PAC.ilon),length(ggrd_PAC.ilat),12)
	prcp_mo_PAC     = zeros(12)
	if typeof(ggrd_PAC) <: PolyGrid
		  mask_PAC = ggrd_PAC.mask
	else; mask_PAC = ones(length(ggrd_PAC.ilon),length(ggrd_PAC.ilat))
	end
	for glat in 1 : length(ggrd_PAC.ilat), glon in 1 : length(ggrd_PAC.ilon)
		prcp_yr_PAC_tmp[glon,glat] = prcp_yr[
			ggrd_PAC.ilon[glon],ggrd_PAC.ilat[glat]
		] * mask_PAC[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_PAC.ilat), glon in 1 : length(ggrd_PAC.ilon)
			prcp_mo_PAC_tmp[glon,glat,imo] = prcp_mo[
				ggrd_PAC.ilon[glon],ggrd_PAC.ilat[glat],imo
			] * mask_PAC[glon,glat]
		end
	end
	prcp_yr_PAC = mean(prcp_yr_PAC_tmp[.!isnan.(prcp_yr_PAC_tmp)])
	for imo = 1 : 12
		prcp_mo_PAC_ii   = @view prcp_mo_PAC_tmp[:,:,imo]
		prcp_mo_PAC[imo] = mean(prcp_mo_PAC_ii[.!isnan.(prcp_mo_PAC_ii)])
	end
	md"Extracting information for the **$(geo_PAC.name)** region ..."
end

# ╔═╡ 85a91f32-5b7e-4aba-9358-3c01364bde8a
begin
	ggrd_CRB = RegionGrid(geo_CRB,plon,plat)
	prcp_yr_CRB_tmp = zeros(length(ggrd_CRB.ilon),length(ggrd_CRB.ilat))
	prcp_mo_CRB_tmp = zeros(length(ggrd_CRB.ilon),length(ggrd_CRB.ilat),12)
	prcp_mo_CRB     = zeros(12)
	if typeof(ggrd_CRB) <: PolyGrid
		  mask_CRB = ggrd_CRB.mask
	else; mask_CRB = ones(length(ggrd_CRB.ilon),length(ggrd_CRB.ilat))
	end
	for glat in 1 : length(ggrd_CRB.ilat), glon in 1 : length(ggrd_CRB.ilon)
		prcp_yr_CRB_tmp[glon,glat] = prcp_yr[
			ggrd_CRB.ilon[glon],ggrd_CRB.ilat[glat]
		] * mask_CRB[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_CRB.ilat), glon in 1 : length(ggrd_CRB.ilon)
			prcp_mo_CRB_tmp[glon,glat,imo] = prcp_mo[
				ggrd_CRB.ilon[glon],ggrd_CRB.ilat[glat],imo
			] * mask_CRB[glon,glat]
		end
	end
	prcp_yr_CRB = mean(prcp_yr_CRB_tmp[.!isnan.(prcp_yr_CRB_tmp)])
	for imo = 1 : 12
		prcp_mo_CRB_ii   = @view prcp_mo_CRB_tmp[:,:,imo]
		prcp_mo_CRB[imo] = mean(prcp_mo_CRB_ii[.!isnan.(prcp_mo_CRB_ii)])
	end
	md"Extracting information for the **$(geo_CRB.name)** region ..."
end

# ╔═╡ 56564ea6-9daf-4cf2-a58c-d93e97cbc3cf
begin
	ggrd_EAN = RegionGrid(geo_EAN,plon,plat)
	prcp_yr_EAN_tmp = zeros(length(ggrd_EAN.ilon),length(ggrd_EAN.ilat))
	prcp_mo_EAN_tmp = zeros(length(ggrd_EAN.ilon),length(ggrd_EAN.ilat),12)
	prcp_mo_EAN     = zeros(12)
	if typeof(ggrd_EAN) <: PolyGrid
		  mask_EAN = ggrd_EAN.mask
	else; mask_EAN = ones(length(ggrd_EAN.ilon),length(ggrd_EAN.ilat))
	end
	for glat in 1 : length(ggrd_EAN.ilat), glon in 1 : length(ggrd_EAN.ilon)
		prcp_yr_EAN_tmp[glon,glat] = prcp_yr[
			ggrd_EAN.ilon[glon],ggrd_EAN.ilat[glat]
		] * mask_EAN[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_EAN.ilat), glon in 1 : length(ggrd_EAN.ilon)
			prcp_mo_EAN_tmp[glon,glat,imo] = prcp_mo[
				ggrd_EAN.ilon[glon],ggrd_EAN.ilat[glat],imo
			] * mask_EAN[glon,glat]
		end
	end
	prcp_yr_EAN = mean(prcp_yr_EAN_tmp[.!isnan.(prcp_yr_EAN_tmp)])
	for imo = 1 : 12
		prcp_mo_EAN_ii   = @view prcp_mo_EAN_tmp[:,:,imo]
		prcp_mo_EAN[imo] = mean(prcp_mo_EAN_ii[.!isnan.(prcp_mo_EAN_ii)])
	end
	md"Extracting information for the **$(geo_EAN.name)** region ..."
end

# ╔═╡ bed9af31-6789-4977-811f-36bcf4d1633b
begin
	ggrd_WAN = RegionGrid(geo_WAN,plon,plat)
	prcp_yr_WAN_tmp = zeros(length(ggrd_WAN.ilon),length(ggrd_WAN.ilat))
	prcp_mo_WAN_tmp = zeros(length(ggrd_WAN.ilon),length(ggrd_WAN.ilat),12)
	prcp_mo_WAN     = zeros(12)
	if typeof(ggrd_WAN) <: PolyGrid
		  mask_WAN = ggrd_WAN.mask
	else; mask_WAN = ones(length(ggrd_WAN.ilon),length(ggrd_WAN.ilat))
	end
	for glat in 1 : length(ggrd_WAN.ilat), glon in 1 : length(ggrd_WAN.ilon)
		prcp_yr_WAN_tmp[glon,glat] = prcp_yr[
			ggrd_WAN.ilon[glon],ggrd_WAN.ilat[glat]
		] * mask_WAN[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_WAN.ilat), glon in 1 : length(ggrd_WAN.ilon)
			prcp_mo_WAN_tmp[glon,glat,imo] = prcp_mo[
				ggrd_WAN.ilon[glon],ggrd_WAN.ilat[glat],imo
			] * mask_WAN[glon,glat]
		end
	end
	prcp_yr_WAN = mean(prcp_yr_WAN_tmp[.!isnan.(prcp_yr_WAN_tmp)])
	for imo = 1 : 12
		prcp_mo_WAN_ii   = @view prcp_mo_WAN_tmp[:,:,imo]
		prcp_mo_WAN[imo] = mean(prcp_mo_WAN_ii[.!isnan.(prcp_mo_WAN_ii)])
	end
	md"Extracting information for the **$(geo_WAN.name)** region ..."
end

# ╔═╡ 3afdbf64-23de-4626-81dc-024c1e937b15
begin
	ggrd_VAL = RegionGrid(geo_VAL,plon,plat)
	prcp_yr_VAL_tmp = zeros(length(ggrd_VAL.ilon),length(ggrd_VAL.ilat))
	prcp_mo_VAL_tmp = zeros(length(ggrd_VAL.ilon),length(ggrd_VAL.ilat),12)
	prcp_mo_VAL     = zeros(12)
	if typeof(ggrd_VAL) <: PolyGrid
		  mask_VAL = ggrd_VAL.mask
	else; mask_VAL = ones(length(ggrd_VAL.ilon),length(ggrd_VAL.ilat))
	end
	for glat in 1 : length(ggrd_VAL.ilat), glon in 1 : length(ggrd_VAL.ilon)
		prcp_yr_VAL_tmp[glon,glat] = prcp_yr[
			ggrd_VAL.ilon[glon],ggrd_VAL.ilat[glat]
		] * mask_VAL[glon,glat]
	end
	for imo = 1 : 12
		for glat in 1 : length(ggrd_VAL.ilat), glon in 1 : length(ggrd_VAL.ilon)
			prcp_mo_VAL_tmp[glon,glat,imo] = prcp_mo[
				ggrd_VAL.ilon[glon],ggrd_VAL.ilat[glat],imo
			] * mask_VAL[glon,glat]
		end
	end
	prcp_yr_VAL = mean(prcp_yr_VAL_tmp[.!isnan.(prcp_yr_VAL_tmp)])
	for imo = 1 : 12
		prcp_mo_VAL_ii   = @view prcp_mo_VAL_tmp[:,:,imo]
		prcp_mo_VAL[imo] = mean(prcp_mo_VAL_ii[.!isnan.(prcp_mo_VAL_ii)])
	end
	md"Extracting information for the **$(geo_VAL.name)** region ..."
end

# ╔═╡ 8699192d-1efd-436d-9147-bd05128c7833
begin
	pplt.close()
	f3,a3 = pplt.subplots([1,1,1,1,1,1,1,1,1,1,2],aspect=1.5,wspace=0.1,sharex=0)
	
	a3[1].plot(0:13,vcat(prcp_mo_CIS[end],prcp_mo_CIS,prcp_mo_CIS[1]),c="k")
	a3[1].plot(0:13,vcat(prcp_mo_PAC[end],prcp_mo_PAC,prcp_mo_PAC[1]))
	a3[1].plot(0:13,vcat(prcp_mo_CRB[end],prcp_mo_CRB,prcp_mo_CRB[1]))
	a3[1].plot(0:13,vcat(prcp_mo_EAN[end],prcp_mo_EAN,prcp_mo_EAN[1]))
	a3[1].plot(0:13,vcat(prcp_mo_VAL[end],prcp_mo_VAL,prcp_mo_VAL[1]))
	a3[1].plot(0:13,vcat(prcp_mo_WAN[end],prcp_mo_WAN,prcp_mo_WAN[1]))
	a3[1].format(xlabel="Month of Year",xlim=(0.5,12.5),ylim=(0,40),xlocator=1:12)
	
	a3[2].scatter(0.5,prcp_yr_CIS,s=10,label="Entire Domain",legend="r",c="k")
	a3[2].scatter(0.5,prcp_yr_PAC,s=10,label="Pacific Coast",legend="r")
	a3[2].scatter(0.5,prcp_yr_CRB,s=10,label="Caribbean Cost",legend="r")
	a3[2].scatter(0.5,prcp_yr_EAN,s=10,label="East Andes",legend="r")
	a3[2].scatter(0.5,prcp_yr_VAL,s=10,label="Andes Valley",legend="r")
	a3[2].scatter(0.5,prcp_yr_WAN,s=10,label="West Andes",legend="r",legend_kw=Dict("ncols"=>1,"frame"=>false))
	a3[2].format(ylabel=L"Rainfall Rate / mm day$^{-1}$",xlocator=[])
	
	f3.savefig(plotsdir("01c-era5tp-ts-subgeo.png"),transparent=false,dpi=300)
	load(plotsdir("01c-era5tp-ts-subgeo.png"))
end

# ╔═╡ Cell order:
# ╟─8bd812fc-5c56-11ec-3514-f7b178ccb922
# ╟─8e761852-a62c-41a8-856f-7d3b3fba9c20
# ╟─682367f2-9bef-4512-a741-20d5c43c6311
# ╟─f8dd4814-bf9d-4dbb-9d45-a0413e4710d6
# ╟─b4b236b7-233c-4c6a-b3ad-4c93fcb74367
# ╟─5ef3c8c5-ef45-4695-8423-4328e43df484
# ╟─3b70530b-d518-4de8-ab33-29c050452ca6
# ╟─09e38bfc-866c-4005-bed2-436c4f22d97d
# ╟─870dbe5f-4437-40de-ba6e-581249a7bf51
# ╟─844cecf6-c06b-4178-9d2f-1c3d25a27252
# ╟─68df62c9-4eed-4829-9ea5-199937435615
# ╟─409d48fb-8a91-497a-9d52-3cb98797dacc
# ╟─c3a4c472-88cb-4573-82a8-30cc3ed8ad97
# ╟─7902b269-124e-47c5-9745-be12c0b5b431
# ╟─3ef9288b-c02a-4918-ac67-76036be08ceb
# ╟─2c75ada3-dff0-45f1-85a2-d8d08327ec14
# ╟─d15b987d-fed8-49ae-a077-b1715da02d14
# ╟─80a51754-af0a-4756-a251-922c9e5d1501
# ╟─20532d37-83eb-4260-8658-d1261b8a23d4
# ╟─b016350b-adc7-4b24-b7c5-d54c61fed723
# ╟─70fef982-f550-4ac5-b1be-21e7e06380e7
# ╟─79be324c-701c-46c3-b7e6-33e2ac5d75c8
# ╟─831c1668-7228-479d-ae1e-cdd176e8010a
# ╟─85a91f32-5b7e-4aba-9358-3c01364bde8a
# ╟─56564ea6-9daf-4cf2-a58c-d93e97cbc3cf
# ╟─bed9af31-6789-4977-811f-36bcf4d1633b
# ╟─3afdbf64-23de-4626-81dc-024c1e937b15
# ╟─8699192d-1efd-436d-9147-bd05128c7833
