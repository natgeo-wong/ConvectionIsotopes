### A Pluto.jl notebook ###
# v0.19.46

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
	using NASAPrecipitation
	using PlutoUI
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
npdv7 = IMERGMonthly(start=Date(2001),stop=Date(2020,12,31),path=datadir())

# ╔═╡ 379d74ef-e4e4-4879-8b9a-92cacb26dc3a
npdv6 = IMERGMonthly(start=Date(2001),stop=Date(2020,12,31),path=datadir(),v6=true)

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ a01e1721-23cf-4a3f-b5aa-0189b1a113a3
lsd = getLandSea(npdv7,geo)

# ╔═╡ 7c466ad9-7c99-4af0-bbc2-c81d8575be6e
begin
	dsv7 = read(npdv7,geo,Date(2020))
	prcpv7 = dsv7["precipitation"][:,:,:] * 86400
	close(dsv7)
end

# ╔═╡ 707fa69c-1eaf-46d4-b0c1-dc810c01e777
begin
	dsv6 = read(npdv6,geo,Date(2020))
	prcpv6 = dsv6["precipitation"][:,:,:] * 86400
	close(dsv6)
end

# ╔═╡ 545bc46f-07e9-422f-a639-5f0a29b3798c
imo = 1

# ╔═╡ d7364f79-67ad-4155-b37a-6c5d59d83a52
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=1.5)

	δprcp = 1 .- (prcpv6[:,:,imo]./prcpv7[:,:,imo])
	δprcp[prcpv7[:,:,imo].<5] .= NaN
	
	axs[1].pcolormesh(lsd.lon,lsd.lat,δprcp',cmap="drywet",levels=vcat(-10:10)/10,extend="both")
	axs[2].pcolormesh(lsd.lon,lsd.lat,prcpv6[:,:,imo]',levels=2:2:40,cmap="blues",extend="both")
	axs[3].pcolormesh(lsd.lon,lsd.lat,prcpv7[:,:,imo]',levels=2:2:40,cmap="blues",extend="both")

	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-95,-65),ylim=(-10,20))
	end
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 64dcc985-e4b1-4583-b502-eda08850522d
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=3,axwidth=1.5)

	δprcpyr = 1 .- dropdims(mean(prcpv6./prcpv7,dims=3),dims=3)
	δprcpyr[dropdims(mean(prcpv7,dims=3),dims=3).<5] .= NaN
	
	a2[1].pcolormesh(lsd.lon,lsd.lat,δprcpyr',cmap="drywet",levels=vcat(-10:10)/10,extend="both")
	a2[2].pcolormesh(lsd.lon,lsd.lat,dropdims(mean(prcpv6,dims=3),dims=3)',levels=2:2:40,cmap="blues",extend="both")
	a2[3].pcolormesh(lsd.lon,lsd.lat,dropdims(mean(prcpv7,dims=3),dims=3)',levels=2:2:40,cmap="blues",extend="both")

	for ax in a2
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-95,-65),ylim=(-10,20))
	end
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─1cee44fe-75b7-42a2-948e-db330cf788e8
# ╟─379d74ef-e4e4-4879-8b9a-92cacb26dc3a
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─a01e1721-23cf-4a3f-b5aa-0189b1a113a3
# ╟─7c466ad9-7c99-4af0-bbc2-c81d8575be6e
# ╟─707fa69c-1eaf-46d4-b0c1-dc810c01e777
# ╠═545bc46f-07e9-422f-a639-5f0a29b3798c
# ╟─d7364f79-67ad-4155-b37a-6c5d59d83a52
# ╟─64dcc985-e4b1-4583-b502-eda08850522d
