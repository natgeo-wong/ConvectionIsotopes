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
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 02a. Exploring the GPM IMERG Dataset for the OTREC Domain

In this notebook, we explore the diurnal cycle of the GPM IMERG dataset for the OTREC domain using a gaussian-filtered/smoothed version of the IMERG landsea mask.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Datasets and LandSea Masks ...
"

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 1cee44fe-75b7-42a2-948e-db330cf788e8
npd = IMERGFinalHH(start=Date(2001),stop=Date(2020,12,31),path=datadir())

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC")

# ╔═╡ a01e1721-23cf-4a3f-b5aa-0189b1a113a3
lsd = getLandSea(npd,geo,smooth=true,σlon=5,σlat=5)

# ╔═╡ 2cb56837-ba5f-476d-aba6-7654c28cb3a8
begin
	ds = NCDataset(datadir("imergfinalhh-OTREC-compiled.nc"))
	prcp = ds["precipitation"][:] * 3600
	close(ds)
	md"Loading compiled half-hourly precipitation ..."
end

# ╔═╡ 14ddca58-cd31-4b52-902e-7d7131d97cff
begin
	ggrd_PAC = RegionGrid(GeoRegion("OTREC_PAC"),lsd.lon,lsd.lat)
	ggrd_ATL = RegionGrid(GeoRegion("OTREC_ATL"),lsd.lon,lsd.lat)
	ggrd_PCS = RegionGrid(GeoRegion("OTREC_PCS"),lsd.lon,lsd.lat)
	ggrd_SAN = RegionGrid(GeoRegion("OTREC_SAN"),lsd.lon,lsd.lat)
end

# ╔═╡ 08c34e26-adbb-4d37-978d-eae09323c53a
begin
	prcp_PAC = extractGrid(prcp,ggrd_PAC)
	prcp_ATL = extractGrid(prcp,ggrd_ATL)
	prcp_PCS = extractGrid(prcp,ggrd_PCS)
	prcp_SAN = extractGrid(prcp,ggrd_SAN)
	md"Extracting precipitation grid for subGeoRegions"
end

# ╔═╡ 65a60484-777d-437d-96ca-42ec4f04aba3
begin
	slon_PCS,slat_PCS = coordGeoRegion(GeoRegion("OTREC_PCS"))
	slon_SAN,slat_SAN = coordGeoRegion(GeoRegion("OTREC_SAN"))
	md"Loading georegion coordinates for San Andres and Pacific Colombia Coast ..."
end

# ╔═╡ 54cb64f0-7de1-4101-8c4f-6757f9423ec6
for it = 1 : 48
	pplt.close(); f2,a2 = pplt.subplots(axwidth=3)
	
	c2 = a2[1].pcolormesh(
		lsd.lon,lsd.lat,prcp[:,:,it]',
		levels=10. .^(-1:0.1:0.5),cmap="blues",extend="both"
	)
	a2[1].plot(x,y,lw=1,c="k")
	a2[1].plot(slon_PCS,slat_PCS,c="r",lw=3)
	a2[1].format(
		xlim=(270,285),xlocator=270:2.5:285,xlabel=L"Longitude / $\degree$",
		ylim=(0,15),ylocator=0:2.5:15,ylabel=L"Latitude / $\degree$",
		suptitle="Colombia Local Time: Hour $(@sprintf("%04.1f",mod(it/2-5,24)))"
	)

	f2.colorbar(c2,label=L"Rainfall Rate / mm hr$^{-1}$")
	mkpath(plotsdir("02a-prcpanim"))
	f2.savefig(plotsdir("02a-prcpanim/$it.png"),transparent=false,dpi=200)
end

# ╔═╡ a24911dc-7b4a-4988-a63b-5797e8b9e283
load(plotsdir("02a-prcpanim/25.png"))

# ╔═╡ 3889238e-bf9c-41a5-a835-bfc638fee76f
md"
### C. Binning Precipitation according to Land-Sea Mask
"

# ╔═╡ 25c9f007-d2f8-40f4-9e32-6da66b292424
lsmvec = vcat(0,collect(10. .^(-4:0.125:-1)),0.15,0.2,0.25,0.5,0.9,0.95,0.99,1)

# ╔═╡ d4c749db-f107-4ee7-a51f-be5f0f9db5ed
begin
	prcpbin_PAC = zeros(48,length(lsmvec)-1); lsm_PAC = extractGrid(lsd.lsm,ggrd_PAC)
	prcpbin_ATL = zeros(48,length(lsmvec)-1); lsm_ATL = extractGrid(lsd.lsm,ggrd_ATL)
	prcpbin_PCS = zeros(48,length(lsmvec)-1); lsm_PCS = extractGrid(lsd.lsm,ggrd_PCS)
	prcpbin_SAN = zeros(48,length(lsmvec)-1); lsm_SAN = extractGrid(lsd.lsm,ggrd_SAN)
	for ilsm = 1 : (length(lsmvec)-1)
		
		ii_PAC = (lsm_PAC .< lsmvec[ilsm+1]) .& (lsm_PAC .> lsmvec[ilsm])
		ii_ATL = (lsm_ATL .< lsmvec[ilsm+1]) .& (lsm_ATL .> lsmvec[ilsm])
		ii_PCS = (lsm_PCS .< lsmvec[ilsm+1]) .& (lsm_PCS .> lsmvec[ilsm])
		ii_SAN = (lsm_SAN .< lsmvec[ilsm+1]) .& (lsm_SAN .> lsmvec[ilsm])
	
		for it = 1 : 48
	
			prcpbin_PAC[it,ilsm] = mean((prcp_PAC[:,:,it])[ii_PAC])
			prcpbin_ATL[it,ilsm] = mean((prcp_ATL[:,:,it])[ii_ATL])
			prcpbin_PCS[it,ilsm] = mean((prcp_PCS[:,:,it])[ii_PCS])
			prcpbin_SAN[it,ilsm] = mean((prcp_SAN[:,:,it])[ii_SAN])
	
		end
	
	end
	prcpbin_PAC = circshift(prcpbin_PAC,(-10,0))
	prcpbin_ATL = circshift(prcpbin_ATL,(-10,0))
	prcpbin_PCS = circshift(prcpbin_PCS,(-10,0))
	prcpbin_SAN = circshift(prcpbin_SAN,(-10,0))
	md"Binning precipitation by landsea mask ..."
end

# ╔═╡ 4a0c88d8-f366-44e8-836c-a55691c49790
begin
	pplt.close(); f1,a1 = pplt.subplots(nrows=4,aspect=4,axwidth=5)

	lsmplot = vcat(0,0.5:0.25:6.5,6.5+1/3,6.5+2/3,7.5,8,8.5,9,9.5,10)
	
	c1 = a1[1].pcolormesh(
		lsmplot,0:0.5:24,prcpbin_PAC,levels=10. .^(-1:0.1:0.5),
		cmap="blues",extend="both"
	); a1[1].format(ultitle="(a) PAC")
	
	a1[2].pcolormesh(
		lsmplot,0:0.5:24,prcpbin_ATL,levels=10. .^(-1:0.1:0.5),
		cmap="blues",extend="both"
	); a1[2].format(ultitle="(b) ATL")
	
	a1[3].pcolormesh(
		lsmplot,0:0.5:24,prcpbin_PCS,levels=10. .^(-1:0.1:0.5),
		cmap="blues",extend="both"
	); a1[3].format(ultitle="(c) PCS")
	
	a1[4].pcolormesh(
		lsmplot,0:0.5:24,prcpbin_SAN,levels=10. .^(-1:0.1:0.5),
		cmap="blues",extend="both"
	); a1[4].format(ultitle="(d) SAN")

	for ax in a1
		ax.plot([8,8],[0,24],c="k",linestyle="-")
		ax.text(6.75,20,"Coastline",c="k")
		ax.format(
			ylocator=0:3:24,ylabel="Hour of Day",
			yticklabels=["0","3","6","9","12","15","18","21","0"],
			xlim=(0,10),xlocator=vcat(0,0.5:9.5,10),xlabel="Smoothed Land-Sea Mask",
			xticklabels=["0","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","0.1","0.25","0.9","0.99","1"]
		)
	end
	
	f1.colorbar(
		c1,label=L"Rain Rate / mm hr$^{-1}$",ticklabels=[L"10^{-1}",L"10^{-0.5}","1",L"10^{0.5}"],
		locator=10. .^(-1:0.5:1),minorlocator=10. .^(-1:0.1:1),length=0.5
	)
	
	f1.savefig(plotsdir("02a-precipitationbin.png"),transparent=false,dpi=400)
	load(plotsdir("02a-precipitationbin.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─1cee44fe-75b7-42a2-948e-db330cf788e8
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─a01e1721-23cf-4a3f-b5aa-0189b1a113a3
# ╟─2cb56837-ba5f-476d-aba6-7654c28cb3a8
# ╟─14ddca58-cd31-4b52-902e-7d7131d97cff
# ╟─08c34e26-adbb-4d37-978d-eae09323c53a
# ╟─65a60484-777d-437d-96ca-42ec4f04aba3
# ╟─54cb64f0-7de1-4101-8c4f-6757f9423ec6
# ╠═a24911dc-7b4a-4988-a63b-5797e8b9e283
# ╟─3889238e-bf9c-41a5-a835-bfc638fee76f
# ╟─25c9f007-d2f8-40f4-9e32-6da66b292424
# ╟─d4c749db-f107-4ee7-a51f-be5f0f9db5ed
# ╟─4a0c88d8-f366-44e8-836c-a55691c49790
