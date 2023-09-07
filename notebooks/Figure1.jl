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
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using PlutoUI

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# Figure 1. Regions of Interest
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
	infocr  = stninfocostarica()
	md"Loading station location information ..."
end

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 60d38ccd-6538-426f-985e-9d782834dce9
begin
	xinset = deepcopy(x); xinset[(x.>-75).|(x.<-90).|(y.>15).|(y.<0)].=NaN
	yinset = deepcopy(y); yinset[(x.>-75).|(x.<-90).|(y.>15).|(y.<0)].=NaN
	md"Loading coastlines data ..."
end

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC")

# ╔═╡ 194bb86d-a3bb-4585-ba66-067dd4488189
geo_large = GeoRegion("OTREC_LRG")

# ╔═╡ 5d6c3cd6-e406-461d-a226-20022060398d
begin
	lsd = getLandSea(
		geo,path=datadir(),savelsd=true,returnlsd=true,
		smooth=true,σlon=30,σlat=30
	)
end

# ╔═╡ 76e35851-71a0-4b89-98bb-9a82bf34bd34
begin
	lsd_large = getLandSea(
		geo_large,path=datadir(),savelsd=true,returnlsd=true,
	)
end

# ╔═╡ f252a060-111a-4e86-85dd-46257c258b77
md"
### B. Loading GeoRegion Information
"

# ╔═╡ 47127cc3-442e-4fa7-b016-82ca8b6ab87a
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

# ╔═╡ 1a643e58-39c1-4c6b-b340-978056871b6b
md"
### C. Plotting Regions of Interest
"

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=27/17,axwidth=5)

	lvls = vcat(10. .^(-5:0.5:-1))
	lvls = vcat(lvls,0.25,0.5,0.9,0.95,0.99)

	cmap_dict = Dict("right"=>0.65)
	
	axs[1].pcolormesh(
		lsd_large.lon.-360,lsd_large.lat,lsd_large.lsm',alpha=0.6,
		levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict
	)
	
	c = axs[1].pcolormesh(
		lsd.lon.-360,lsd.lat,lsd.lsm',
		levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict
	)
	
	axs[1].plot(xinset,yinset,lw=0.5,c="k")
	axs[1].scatter(infomo[:,2],infomo[:,3],zorder=4,c="blue")
	axs[1].scatter(infody[:,2],infody[:,3],zorder=4,c="pink")
	axs[1].scatter(-infocr[:,3],infocr[:,2],zorder=4,c="red")
	axs[1].plot([-90,-75,-75,-90,-90],[0,0,15,15,0],lw=5,c="w")
	axs[1].plot(slon_PAC,slat_PAC,lw=2,)
	axs[1].plot(slon_ATL,slat_ATL,lw=2,)
	axs[1].plot(slon_PCS,slat_PCS,lw=2,)
	axs[1].plot(slon_SAN,slat_SAN,lw=2,)
	axs[1].text(-88.5,1.5,"PAC",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
	axs[1].text(-78,13,"ATL",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
	axs[1].text(-81.7,10.8,"SAN",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
	axs[1].text(-84,1.5,"PCS",bbox=Dict("facecolor"=>"gray2","edgecolor"=>"none"))
	axs[1].text(-74.5,14.5,"(a)",c="w",size=10)
	axs[1].text(-66,7.6,"(b)",c="k",size=10)
	axs[1].format(
		xlim=(-91,-64),xlocator=-105:5:-60,xminorlocator=-120:2.5:-45,
		ylim=(-1,16),ylocator=-15:5:30,yminorlocator=-20:2.5:35,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
		grid=true,gridcolor="w",#suptitle="Available Colombia Stations",
	)

	ix = fig.add_axes([0.617,0.175,0.2,0.36])
	ix.pcolormesh(
		lsd_large.lon.-360,lsd_large.lat,lsd_large.lsm',
		levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict
	)
	ix.pcolormesh(
		lsd.lon.-360,lsd.lat,lsd.lsm',
		levels=lvls,cmap="delta",extend="both",cmap_kw=cmap_dict
	)
	ix.plot(xinset,yinset,lw=0.5,c="k")
	ix.plot([-90,-75,-75,-90,-90],[0,0,15,15,0],lw=1.5,c="w")
	ix.plot([-105,-60,-60,-105,-105],[-15,-15,30,30,-15],lw=1.5,c="w",linestyle="--")
	ix.text(-100,0,"d02",c="w",size=8)
	ix.text(-70,23,"d01",c="w",size=8)
	ix.format(
		xlim=(-110,-55),ylim=(-20,35),xtickloc="none",ytickloc="none",
		xlocator=-105:15:-60,xminorlocator=-120:5:-45,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	fig.colorbar(c,ticklabels=["1e-5","1e-4.5","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1","0.25","0.5","0.9","0.95","0.99"],label="Land-Sea Mask")
	fig.savefig(plotsdir("Figure1.png"),transparent=false,dpi=300)
	load(plotsdir("Figure1.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╠═189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─60d38ccd-6538-426f-985e-9d782834dce9
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─194bb86d-a3bb-4585-ba66-067dd4488189
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╟─76e35851-71a0-4b89-98bb-9a82bf34bd34
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╟─47127cc3-442e-4fa7-b016-82ca8b6ab87a
# ╟─601f4ad5-3520-4555-af6b-8eb2a299609c
# ╟─1a643e58-39c1-4c6b-b340-978056871b6b
# ╠═d7755534-3565-4011-b6e3-e131991008db
