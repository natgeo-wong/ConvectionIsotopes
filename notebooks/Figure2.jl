### A Pluto.jl notebook ###
# v0.19.36

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
# Figure 2. OTREC Region of Interest
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
	xinset = deepcopy(x); xinset[(x.>285).|(x.<270).|(y.>15).|(y.<0)].=NaN
	yinset = deepcopy(y); yinset[(x.>285).|(x.<270).|(y.>15).|(y.<0)].=NaN
	md"Filtering out coastlines outside domain ..."
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

# ╔═╡ 1a643e58-39c1-4c6b-b340-978056871b6b
md"
### C. Plotting Regions of Interest
"

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=27/17,axwidth=5)

	lvls = -6 : 6
	textdict = Dict("fc"=>"grey3","ec"=>"none","alpha"=>0.6)
	
	axs[1].pcolormesh(
		lsd_large.lon,lsd_large.lat,lsd_large.z'/1000,alpha=0.3,
		levels=lvls,cmap="delta",extend="both"
	)
	
	c = axs[1].pcolormesh(
		lsd.lon,lsd.lat,lsd.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)
	
	axs[1].plot(xinset,yinset,lw=0.5,c="k")
	axs[1].plot([275,275,278,278,275],[8,10.5,10.5,8,8],lw=1,c="k",linestyle="--")
	axs[1].scatter(infody[:,2],infody[:,3],zorder=4,c="pink")
	axs[1].scatter(infocr[:,2],infocr[:,3],zorder=4,c="red")
	axs[1].plot([275,288],[10.5,15],lw=1,c="k",linestyle=":")
	axs[1].plot([278,295],[8,9.03],lw=1,c="k",linestyle=":")
	axs[1].plot([270,285,285,270,270],[0,0,15,15,0],lw=5,c="k")
	axs[1].text(285.5,14.5,"(a)",c="k",size=10)
	axs[1].text(286.8,14.5,"(b)",c="k",size=10)
	axs[1].text(286.8,6.5,"(c)",c="k",size=10)

	axs[1].text(272.3,10,"Liberia",c="k",size=8,bbox=textdict)
	axs[1].text(277.5,13,"San Andres",c="k",size=8,bbox=textdict)
	axs[1].text(278.6,6,"Bahía Solano",c="k",size=8,bbox=textdict)
	axs[1].text(282,4.75,"Quibdó",c="k",size=8,bbox=textdict)
	axs[1].text(278.8,3.5,"Buenaventura",c="k",size=8,bbox=textdict)
	
	axs[1].format(
		xlim=(269,296),xlocator=255:5:300,xminorlocator=240:2.5:315,
		ylim=(-1,16),ylocator=-15:5:30,yminorlocator=-20:2.5:35,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
		grid=true,gridcolor="w",#suptitle="Available Colombia Stations",
	)

	ix = fig.add_axes([0.635,0.632,0.21,0.30])
	ix.pcolormesh(
		lsd.lon,lsd.lat,lsd.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)
	ix.plot(xinset,yinset,lw=0.5,c="k")
	ix.scatter(infocr[:,2],infocr[:,3],s=10,zorder=4,c="red")

	ix.text(275.2,10.15,"EEFMB",c="k",size=7,bbox=textdict)
	ix.text(275.4,9.65,"CGFI",c="k",size=7,bbox=textdict)
	ix.text(276.1,9.64,"AMDQ",c="k",size=7,bbox=textdict)
	ix.text(276.1,8.3,"OSA",c="k",size=7,bbox=textdict)
	ix.text(276.2,10.22,"Bataan",c="k",size=7,bbox=textdict)
	ix.text(277,10.1,"Limon",c="k",size=7,bbox=textdict)
	ix.text(277,9.45,"Cahuita",c="k",size=7,bbox=textdict)
	
	ix.format(
		xlim=(275,278),ylim=(8,10.5),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	ix = fig.add_axes([0.635,0.175,0.21,0.36])
	ix.pcolormesh(
		lsd_large.lon,lsd_large.lat,lsd_large.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)
	ix.pcolormesh(
		lsd.lon,lsd.lat,lsd.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)
	ix.plot(xinset,yinset,lw=0.5,c="k")
	ix.plot([270,285,285,270,270],[0,0,15,15,0],lw=1.5,c="k")
	ix.plot([255,300,300,255,255],[-15,-15,30,30,-15],lw=1.5,c="k",linestyle="--")
	ix.text(260,0,"d02",c="k",size=8)
	ix.text(290,-12,"d01",c="k",size=8)
	ix.format(
		xlim=(250,305),ylim=(-20,35),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	fig.colorbar(c,label="Topographic Height/ km")
	fig.savefig(plotsdir("fig2-otrecdomain.png"),transparent=false,dpi=400)
	load(plotsdir("fig2-otrecdomain.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─60d38ccd-6538-426f-985e-9d782834dce9
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─194bb86d-a3bb-4585-ba66-067dd4488189
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╟─76e35851-71a0-4b89-98bb-9a82bf34bd34
# ╟─1a643e58-39c1-4c6b-b340-978056871b6b
# ╟─d7755534-3565-4011-b6e3-e131991008db
