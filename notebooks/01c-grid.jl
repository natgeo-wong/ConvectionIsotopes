### A Pluto.jl notebook ###
# v0.19.45

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
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01a. Creating GeoRegions

In this notebook, we define additional GeoRegions of interest for plotting and for analysis based on WRF modelling output and as necessary for figures.
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

# ╔═╡ f252a060-111a-4e86-85dd-46257c258b77
md"
### B. Loading WRF Grid
"

# ╔═╡ f494290b-a2a3-47c9-80ed-44a1842e5c42
begin
	ds  = NCDataset(datadir("wrf2","grid.nc"))
	lon = ds["longitude"][:,:]
	lat = ds["latitude"][:,:]
	close(ds)
end

# ╔═╡ f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
md"
### C. Finding the Centroid and Tilt of the Grid
"

# ╔═╡ 7a2645c6-e78b-476e-862a-b13b0e58bace
begin
	clon = mean(lon)
	clat = mean(lat)
	lon1 = lon[1,1]; lon2 = lon[end,1]; lon3 = lon[1,end]; lon4 = lon[end,end]
	lat1 = lat[1,1]; lat2 = lat[end,1]; lat3 = lat[1,end]; lat4 = lat[end,end]
	tilt = (atand(lat2-lat1,lon2-lon1) .+ atand(lat4-lat3,lon4-lon3))/2
end

# ╔═╡ 951b4841-257e-41d6-b264-a06a27fe8062
begin
	r = sqrt.((lon .- clon).^2 .+ (lat.-clat).^2)
	θ = atand.(lat.-clat,lon.-clon) .- tilt
	rotX = r .* cosd.(θ)
	rotY = r .* sind.(θ)
end

# ╔═╡ 13afaa2d-05ed-406e-a9ca-d8e6cef40e23
tst = TiltRegion("","","",0,0,12,12,-11,save=false)

# ╔═╡ 4b8484e3-8945-42f4-860a-b16fb4266e73
tlon,tlat = coordGeoRegion(tst)

# ╔═╡ cadcbc38-0c61-493a-85f6-bbba012f44eb
ii = argmin((mod.(lon,360).-infody[1,2]).^2 .+ (lat.-infody[1,3]).^2)

# ╔═╡ 081339ef-8c6f-4333-bdfd-47ac78c02dce
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2)
	
	axs[1].pcolormesh(lon.-clon,lat.-clat,lon)
	axs[1].plot(tlon,tlat)
	axs[1].scatter(lon[ii].-clon,lat[ii].-clat)
	axs[1].pcolormesh(
		lon[ii[1].+(-50:50),ii[2].+(-50:50)].-clon,
		lat[ii[1].+(-50:50),ii[2].+(-50:50)].-clat,
		ones(101,101)
	)
	axs[2].pcolormesh(rotX,rotY,lon)
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╠═f494290b-a2a3-47c9-80ed-44a1842e5c42
# ╟─f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
# ╠═7a2645c6-e78b-476e-862a-b13b0e58bace
# ╠═951b4841-257e-41d6-b264-a06a27fe8062
# ╠═13afaa2d-05ed-406e-a9ca-d8e6cef40e23
# ╠═4b8484e3-8945-42f4-860a-b16fb4266e73
# ╠═081339ef-8c6f-4333-bdfd-47ac78c02dce
# ╠═cadcbc38-0c61-493a-85f6-bbba012f44eb
