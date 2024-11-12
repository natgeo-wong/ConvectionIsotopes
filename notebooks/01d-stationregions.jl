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
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
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
	infody  = stninfody(); ndystn = size(infody,1)
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
for istn = 1 : ndystn

	lon_stn = infody[istn,2]
	lat_stn = infody[istn,3]

	ii = argmin((mod.(lon,360).-lon_stn).^2 .+ (lat.-lat_stn).^2)
	ilon_stn = ii[1]
	ilat_stn = ii[2]

	ilon = vcat(
		lon[ilon_stn.+(-20:19),ilat_stn-20],
		lon[ilon_stn+20,ilat_stn.+(-20:19)],
		reverse(lon[ilon_stn.+(-19:20),ilat_stn+20]),
		reverse(lon[ilon_stn-20,ilat_stn.+(-19:20)])
	)

	ilat = vcat(
		lat[ilon_stn.+(-20:19),ilat_stn-20],
		lat[ilon_stn+20,ilat_stn.+(-20:19)],
		reverse(lat[ilon_stn.+(-19:20),ilat_stn+20]),
		reverse(lat[ilon_stn-20,ilat_stn.+(-19:20)])
	)

	if isGeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir(),throw=false)
		removeGeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir())
	end
	geo_d01 = PolyRegion(
		"OTREC_wrf_stn$(@sprintf("%02d",istn))","GLB","WRF Station $(@sprintf("%02d",istn))",
		ilon,ilat,path=srcdir()
	)

	for iboxlat = 1 : 2, iboxlon = 1 : 2

		ibox = iboxlon + (iboxlat-1)*2
		ilon = vcat(
			lon[ilon_stn.+(-20:-1).+(iboxlon-1)*20,ilat_stn-20+(iboxlat-1)*20],
			lon[ilon_stn+(iboxlon-1)*20,ilat_stn.+(-20:-1).+(iboxlat-1)*20],
			reverse(lon[ilon_stn.+(-19:0).+(iboxlon-1)*20,ilat_stn+(iboxlat-1)*20]),
			reverse(lon[ilon_stn-20+(iboxlon-1)*20,ilat_stn.+(-19:0).+(iboxlat-1)*20])
		)
	
		ilat = vcat(
			lat[ilon_stn.+(-20:-1).+(iboxlon-1)*20,ilat_stn-20+(iboxlat-1)*20],
			lat[ilon_stn+(iboxlon-1)*20,ilat_stn.+(-20:-1).+(iboxlat-1)*20],
			reverse(lat[ilon_stn.+(-19:0).+(iboxlon-1)*20,ilat_stn+(iboxlat-1)*20]),
			reverse(lat[ilon_stn-20+(iboxlon-1)*20,ilat_stn.+(-19:0).+(iboxlat-1)*20])
		)

		if isGeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox",path=srcdir(),throw=false)
			removeGeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox",path=srcdir())
		end
		geo_d01 = PolyRegion(
			"OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox","GLB","WRF Station $(@sprintf("%02d",istn)) | Box $ibox",
			ilon,ilat,path=srcdir()
		)

		
	end
	
end

# ╔═╡ 912a101b-7eb9-4322-a713-031aeffff20d
begin
	geo_stn1 =      GeoRegion("OTREC_wrf_stn02",path=srcdir())
	geo_stn1_box1 = GeoRegion("OTREC_wrf_stn02_box1",path=srcdir())
	geo_stn1_box2 = GeoRegion("OTREC_wrf_stn02_box2",path=srcdir())
	geo_stn1_box3 = GeoRegion("OTREC_wrf_stn02_box3",path=srcdir())
	geo_stn1_box4 = GeoRegion("OTREC_wrf_stn02_box4",path=srcdir())
end

# ╔═╡ 25bf5b00-0084-4241-909e-b48a7c9baaaa
begin
	lon1,lat1 = coordGeoRegion(geo_stn1)
	lon2,lat2 = coordGeoRegion(geo_stn1_box1)
	lon3,lat3 = coordGeoRegion(geo_stn1_box2)
	lon4,lat4 = coordGeoRegion(geo_stn1_box3)
	lon5,lat5 = coordGeoRegion(geo_stn1_box4)
end

# ╔═╡ cb306a1f-affd-432c-8c46-977758531654
begin
	pplt.close(); fig,axs = pplt.subplots()
	
	axs[1].plot(lon1,lat1,lw=5)
	axs[1].plot(lon2,lat2)
	axs[1].plot(lon3,lat3)
	axs[1].plot(lon4,lat4)
	axs[1].plot(lon5,lat5)

	axs[1].format(xlim=(-77.5,-75.5),ylim=(4.5,6.5))
	
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
# ╠═32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╠═f494290b-a2a3-47c9-80ed-44a1842e5c42
# ╟─f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
# ╠═7a2645c6-e78b-476e-862a-b13b0e58bace
# ╠═912a101b-7eb9-4322-a713-031aeffff20d
# ╠═25bf5b00-0084-4241-909e-b48a7c9baaaa
# ╠═cb306a1f-affd-432c-8c46-977758531654
