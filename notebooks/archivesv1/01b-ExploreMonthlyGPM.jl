### A Pluto.jl notebook ###
# v0.17.4

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

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ 8bd812fc-5c56-11ec-3514-f7b178ccb922
md"
# 01b - Exploring the Precipitation Climatology of Columbia

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

# ╔═╡ 2ddc5f1b-685b-4584-8aa4-9e2fad3885f4
npdmo = IMERGMonthly(dtbeg=Date(2001,1,1),dtend=Date(2020,1,1),sroot=datadir())

# ╔═╡ c3a4c472-88cb-4573-82a8-30cc3ed8ad97
begin
	pfnc = npdmo.npdID * "-" * geo_OTREC.regID * "-2020.nc"
	pds  = NCDataset(joinpath(npdmo.sroot,geo_OTREC.regID,"raw",pfnc))
	plon = pds["longitude"][:]; nplon = length(plon)
	plat = pds["latitude"][:];  nplat = length(plat)
	prcp_yr = zeros(nplon,nplat)
	prcp_mo = zeros(nplon,nplat,12)
	
	for yr = 1 : 20
		ifnc = npdmo.npdID * "-" * geo_OTREC.regID * "-$(yr+2000).nc"
		ids  = NCDataset(joinpath(npdmo.sroot,geo_OTREC.regID,"raw",ifnc))
		prcp_yr[:,:] += dropdims(mean(ids["prcp_rate"][:],dims=3),dims=3) * 86400
		for mo = 1 : 12
			prcp_mo[:,:,mo] += ids["prcp_rate"][:,:,mo] * 86400
		end
	end
	close(pds)
	
	prcp_yr = prcp_yr / 20
	prcp_mo = prcp_mo / 20
	md"Loading GPM Precipitation data for all years"
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
	
	f1.savefig(plotsdir("01b-gpmclimatology.png"),transparent=false,dpi=300)
	load(plotsdir("01b-gpmclimatology.png"))
end

# ╔═╡ 3ef9288b-c02a-4918-ac67-76036be08ceb
md"Now, we proceed by extracting and plotting the monthly-averages for each of the specified GeoRegions within our domain ..."

# ╔═╡ 2c75ada3-dff0-45f1-85a2-d8d08327ec14
begin
	prcp_yr_CIS,prcp_mo_CIS = extractregclimate(prcp_yr,prcp_mo,geo_CIS,plon,plat)
	prcp_yr_CST,prcp_mo_CST = extractregclimate(prcp_yr,prcp_mo,geo_CST,plon,plat)
	prcp_yr_MNT,prcp_mo_MNT = extractregclimate(prcp_yr,prcp_mo,geo_MNT,plon,plat)
	prcp_yr_INL,prcp_mo_INL = extractregclimate(prcp_yr,prcp_mo,geo_INL,plon,plat)
	md"Extracting information for the large GeoRegions ..."
end

# ╔═╡ b016350b-adc7-4b24-b7c5-d54c61fed723
begin
	pplt.close()
	f2,a2 = pplt.subplots([1,1,1,1,1,1,1,1,1,1,2],aspect=1.5,wspace=0.1,sharex=0)
	
	a2[1].plot(0:13,vcat(prcp_mo_CIS[end],prcp_mo_CIS,prcp_mo_CIS[1]),c="k")
	a2[1].plot(0:13,vcat(prcp_mo_CST[end],prcp_mo_CST,prcp_mo_CST[1]))
	a2[1].plot(0:13,vcat(prcp_mo_MNT[end],prcp_mo_MNT,prcp_mo_MNT[1]))
	a2[1].plot(0:13,vcat(prcp_mo_INL[end],prcp_mo_INL,prcp_mo_INL[1]))
	a2[1].format(xlabel="Month of Year",xlim=(0.5,12.5),ylim=(0,15),xlocator=1:12)
	
	a2[2].scatter(0.5,prcp_yr_CIS,s=10,label="Entire Domain",legend="r",c="k")
	a2[2].scatter(0.5,prcp_yr_CST,s=10,label="Coastal Columbia",legend="r")
	a2[2].scatter(0.5,prcp_yr_MNT,s=10,label="Mountaineous Columbia",legend="r")
	a2[2].scatter(0.5,prcp_yr_INL,s=10,label="Inland Columbia",legend="r",legend_kw=Dict("ncols"=>1,"frame"=>false))
	a2[2].format(ylabel=L"Rainfall Rate / mm day$^{-1}$",xlocator=[])
	
	f2.savefig(plotsdir("01b-prcp_ts-largegeo.png"),transparent=false,dpi=300)
	load(plotsdir("01b-prcp_ts-largegeo.png"))
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
	geo_SAN = GeoRegion("CIS_SAN")
	md"Loading subGeoRegion information ..."
end

# ╔═╡ 23fef6f6-f4f2-4b31-a0b2-c75e23f41499
begin
	prcp_yr_PAC,prcp_mo_PAC = extractregclimate(prcp_yr,prcp_mo,geo_PAC,plon,plat)
	prcp_yr_CRB,prcp_mo_CRB = extractregclimate(prcp_yr,prcp_mo,geo_CRB,plon,plat)
	prcp_yr_EAN,prcp_mo_EAN = extractregclimate(prcp_yr,prcp_mo,geo_EAN,plon,plat)
	prcp_yr_VAL,prcp_mo_VAL = extractregclimate(prcp_yr,prcp_mo,geo_VAL,plon,plat)
	prcp_yr_WAN,prcp_mo_WAN = extractregclimate(prcp_yr,prcp_mo,geo_WAN,plon,plat)
	prcp_yr_SAN,prcp_mo_SAN = extractregclimate(prcp_yr,prcp_mo,geo_SAN,plon,plat)
	md"Extracting information for the subGeoRegions ..."
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
	a3[1].plot(0:13,vcat(prcp_mo_SAN[end],prcp_mo_SAN,prcp_mo_SAN[1]))
	a3[1].format(xlabel="Month of Year",xlim=(0.5,12.5),ylim=(0,20),xlocator=1:12)
	
	a3[2].scatter(0.5,prcp_yr_CIS,s=10,label="Entire Domain",legend="r",c="k")
	a3[2].scatter(0.5,prcp_yr_PAC,s=10,label="Pacific Coast",legend="r")
	a3[2].scatter(0.5,prcp_yr_CRB,s=10,label="Caribbean Cost",legend="r")
	a3[2].scatter(0.5,prcp_yr_EAN,s=10,label="East Andes",legend="r")
	a3[2].scatter(0.5,prcp_yr_VAL,s=10,label="Andes Valley",legend="r")
	a3[2].scatter(0.5,prcp_yr_WAN,s=10,label="West Andes",legend="r")
	a3[2].scatter(0.5,prcp_yr_SAN,s=10,label="San Andres",legend="r",legend_kw=Dict("ncols"=>1,"frame"=>false))
	a3[2].format(ylabel=L"Rainfall Rate / mm day$^{-1}$",xlocator=[])
	
	f3.savefig(plotsdir("01b-prcp_ts_subgeo.png"),transparent=false,dpi=300)
	load(plotsdir("01b-prcp_ts_subgeo.png"))
end

# ╔═╡ 61f6bf0b-5b36-44e2-b7d6-d7de92e64cad
md"
It seems that breaking up the Mountaineous GeoRegion into the various subGeoRegions of the East, West and Valley of the Andes doesn't make much difference to the overall climatology - they show roughly the same patterns, of course with some modulation from the Pacific coastline climatology and inland Columbia climatology as one approaches those regions.
"

# ╔═╡ fc02dcf8-bee3-4c68-93da-aa09f7ec68d9
begin
	blon_PAC,blat_PAC,slon_PAC,slat_PAC = coordGeoRegion(geo_PAC)
	blon_CRB,blat_CRB,slon_CRB,slat_CRB = coordGeoRegion(geo_CRB)
	md"Loading geographic shapes of the subGeoRegions ..."
end

# ╔═╡ 1a0586a9-d0ca-4e76-8ef0-360cdb9d941f
begin
	pplt.close(); fig1,axs1 = pplt.subplots(
		arr,aspect=7.5/13,axwidth=2,
		hspace=[0.1,0.1],wspace=[0,0,0.2,0.1,0.1,0.1]
	)
	
	for iax = 1 : 13
		if !isone(iax)
			c = axs1[iax].contourf(
				plon,plat,prcp_mo[:,:,iax-1]',
				levels=5:2.5:35,extend="both",cmap="Blues"
			)
			if iax == 13
				fig1.colorbar(c,loc="r",label=L"mm day$^{-1}$")
			end
			axs1[iax].plot(slon_PAC,slat_PAC,lw=1)
			axs1[iax].plot(slon_CRB,slat_CRB,lw=1)
			axs1[iax].plot(slon_MNT,slat_MNT,lw=1)
			axs1[iax].plot(slon_INL,slat_INL,lw=1)
			axs1[iax].plot(blon_CIS,blat_CIS,lw=1,c="k",linestyle="--")
			axs1[iax].contour(
				lon,lat,oro'/1000,
				lw=0.5,cmap="brown1",levels=1:5,extend="both"
			)
		else
			c = axs1[iax].contourf(
				plon,plat,prcp_yr',
				levels=5:25,extend="both",cmap="Blues"
			)
			axs1[iax].colorbar(c,loc="l",label=L"mm day$^{-1}$")
			axs1[iax].plot(slon_PAC,slat_PAC,lw=1.5)
			axs1[iax].plot(slon_CRB,slat_CRB,lw=1.5)
			axs1[iax].plot(slon_MNT,slat_MNT,lw=1.5)
			axs1[iax].plot(slon_INL,slat_INL,lw=1.5)
			axs1[iax].plot(blon_CIS,blat_CIS,lw=1.5,c="k",linestyle="--")
			axs1[iax].contour(
				lon,lat,oro'/1000,
				cmap="brown1",levels=1:5,extend="both"
			)
		end
		
		
	end
	
	for ax in axs1
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-78.5,-71),ylim=(-1,12))
	end
	
	axs1[1].format(xlocator=-85:2.5:-70,ylocator=-2.5:2.5:15)
	for iax = 2 : 13
		axs1[iax].format(xlocator=[],ylocator=[],lltitle="$(iax-1)")
	end
	
	axs1[1].format(ltitle="(a) Yearly Precipitation Rate")
	axs1[2].format(ltitle="(b) Monthly Precipitation Rate")
	
	fig1.savefig(plotsdir("01b-gpmclimatologyfinal.png"),transparent=false,dpi=300)
	load(plotsdir("01b-gpmclimatologyfinal.png"))
end

# ╔═╡ b8f3728c-4478-4c42-b27b-a244bc7e3348
begin
	pplt.close()
	fig2,axs2 = pplt.subplots([1,1,1,1,1,1,1,1,1,1,2],aspect=1.5,wspace=0.1,sharex=0)
	
	axs2[1].plot(0:13,vcat(prcp_mo_CIS[end],prcp_mo_CIS,prcp_mo_CIS[1]),c="k")
	axs2[1].plot(0:13,vcat(prcp_mo_PAC[end],prcp_mo_PAC,prcp_mo_PAC[1]))
	axs2[1].plot(0:13,vcat(prcp_mo_CRB[end],prcp_mo_CRB,prcp_mo_CRB[1]))
	axs2[1].plot(0:13,vcat(prcp_mo_MNT[end],prcp_mo_MNT,prcp_mo_MNT[1]))
	axs2[1].plot(0:13,vcat(prcp_mo_INL[end],prcp_mo_INL,prcp_mo_INL[1]))
	axs2[1].plot(0:13,vcat(prcp_mo_SAN[end],prcp_mo_SAN,prcp_mo_SAN[1]))
	axs2[1].format(xlabel="Month of Year",xlim=(0.5,12.5),ylim=(0,20),xlocator=1:12)
	
	axs2[2].scatter(0.5,prcp_yr_CIS,s=10,label="Entire Domain",legend="r",c="k")
	axs2[2].scatter(0.5,prcp_yr_PAC,s=10,label="Pacific Coast",legend="r")
	axs2[2].scatter(0.5,prcp_yr_CRB,s=10,label="Caribbean Cost",legend="r")
	axs2[2].scatter(0.5,prcp_yr_MNT,s=10,label="Mountaineous Columbia",legend="r")
	axs2[2].scatter(0.5,prcp_yr_INL,s=10,label="Inland Columbia",legend="r")
	axs2[2].scatter(0.5,prcp_yr_SAN,s=10,label="San Andres",legend="r",legend_kw=Dict("ncols"=>1,"frame"=>false))
	axs2[2].format(ylabel=L"Rainfall Rate / mm day$^{-1}$",xlocator=[])
	
	fig2.savefig(plotsdir("01b-prcp_ts.png"),transparent=false,dpi=300)
	load(plotsdir("01b-prcp_ts.png"))
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
# ╟─2ddc5f1b-685b-4584-8aa4-9e2fad3885f4
# ╠═c3a4c472-88cb-4573-82a8-30cc3ed8ad97
# ╟─7902b269-124e-47c5-9745-be12c0b5b431
# ╟─3ef9288b-c02a-4918-ac67-76036be08ceb
# ╟─2c75ada3-dff0-45f1-85a2-d8d08327ec14
# ╟─b016350b-adc7-4b24-b7c5-d54c61fed723
# ╟─70fef982-f550-4ac5-b1be-21e7e06380e7
# ╟─79be324c-701c-46c3-b7e6-33e2ac5d75c8
# ╟─23fef6f6-f4f2-4b31-a0b2-c75e23f41499
# ╟─8699192d-1efd-436d-9147-bd05128c7833
# ╟─61f6bf0b-5b36-44e2-b7d6-d7de92e64cad
# ╟─fc02dcf8-bee3-4c68-93da-aa09f7ec68d9
# ╟─1a0586a9-d0ca-4e76-8ef0-360cdb9d941f
# ╟─b8f3728c-4478-4c42-b27b-a244bc7e3348
