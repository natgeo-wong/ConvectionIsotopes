### A Pluto.jl notebook ###
# v0.20.4

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

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

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
### B. Domain Shape for WRF Modelling of OTREC
"

# ╔═╡ ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
begin
	ds_d01  = NCDataset(datadir("wrf3","wrfout","2020","JAN_P01","wrfout3D_d01_2020-01-02_00:00:00"))
	nlon_01 = ds_d01.dim["west_east"]
	nlat_01 = ds_d01.dim["south_north"]
	lon_d01 = vcat(
		ds_d01["XLONG"][1:(nlon_01-1),1,1],
		ds_d01["XLONG"][nlon_01,1:(nlat_01-1),1],
		reverse(ds_d01["XLONG"][2:nlon_01,nlat_01,1]),
		reverse(ds_d01["XLONG"][1,1:nlat_01,1])
	)
	lat_d01 = vcat(
		ds_d01["XLAT"][1:(nlon_01-1),1,1],
		ds_d01["XLAT"][nlon_01,1:(nlat_01-1),1],
		reverse(ds_d01["XLAT"][2:nlon_01,nlat_01,1]),
		reverse(ds_d01["XLAT"][1,1:nlat_01,1])
	)
	close(ds_d01)
end

# ╔═╡ d5bc84d7-6516-4fad-a245-703e14e22f77
begin
	ds_d02  = NCDataset(datadir("wrf3","wrfout","2020","JAN_P01","wrfout3D_d02_2020-01-02_00:00:00"))
	nlon_02 = ds_d02.dim["west_east"]
	nlat_02 = ds_d02.dim["south_north"]
	lon_d02 = vcat(
		ds_d02["XLONG"][1:(nlon_02-1),1,1],
		ds_d02["XLONG"][nlon_02,1:(nlat_02-1),1],
		reverse(ds_d02["XLONG"][2:nlon_02,nlat_02,1]),
		reverse(ds_d02["XLONG"][1,1:nlat_02,1])
	)
	lat_d02 = vcat(
		ds_d02["XLAT"][1:(nlon_02-1),1,1],
		ds_d02["XLAT"][nlon_02,1:(nlat_02-1),1],
		reverse(ds_d02["XLAT"][2:nlon_02,nlat_02,1]),
		reverse(ds_d02["XLAT"][1,1:nlat_02,1])
	)
	close(ds_d02)
end

# ╔═╡ 7853fab4-b431-4872-9645-581e06728d19
begin
	pplt.close(); f1,a1 = pplt.subplots()
	
	a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].plot(lon_d01,lat_d01)
	a1[1].plot(lon_d02,lat_d02)
	a1[1].format(xlim=(-105,-55),ylim=(-20,30))
	
	f1.savefig(plotsdir("01a-tiltedwrfdomain.png"),transparent=false,dpi=200)
	load(plotsdir("01a-tiltedwrfdomain.png"))
end

# ╔═╡ b2dc9507-7794-4445-b7f2-4351cba54af0
if !isID("OTREC_wrf_d01",path=srcdir(),throw=false)
	geo_d01 = GeoRegion(
		lon_d01,lat_d01, save = true,
		ID = "OTREC_wrf_d01", pID = "GLB", name = "WRF d01", path = srcdir()
	)
else
	geo_d01 = GeoRegion("OTREC_wrf_d01",path=srcdir())
end

# ╔═╡ 12766ec3-8b45-49f5-908f-697322059df5
if !isID("OTREC_wrf_d02",path=srcdir(),throw=false)
	geo_d02 = GeoRegion(
		lon_d02,lat_d02, save = true,
		ID = "OTREC_wrf_d02", pID = "GLB", name = "WRF d02", path = srcdir()
	)
else
	geo_d02 = GeoRegion("OTREC_wrf_d02",path=srcdir())
end

# ╔═╡ 9fba6b0f-5d17-4e3a-b9da-de815b12a6cf
md"
### C. Defining Regions for Plotting
"

# ╔═╡ 29fa97df-30dc-4258-b1e7-22662151919f
if !isID("Fig2_small",path=srcdir(),throw=false)
	geo_fig2sml = GeoRegion(
		ID = "Fig2_small", pID = "GLB", name = "Fig2 Small Domain",
		[275,275,278,278,275], [8,10.5,10.5,8,8],
		path=srcdir(), save=true
	)
else
	geo_fig2sml = GeoRegion("Fig2_small",path=srcdir())
end


# ╔═╡ a8f2c79d-d79d-4163-ad98-68578937629b
if !isID("Fig2_big",path=srcdir(),throw=false)
	geo_fig2big = GeoRegion(
		ID = "Fig2_big", pID = "GLB", name = "Fig2 Big Domain",
		[250,250,305,305,250], [-20,35,35,-20,-20],
		path=srcdir(), save=true
	)
else
	geo_fig2big = GeoRegion("Fig2_big",path=srcdir())
end

# ╔═╡ dcc32199-77ac-4b8a-9fec-8d83cbf00cf2
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=27/17,axwidth=5)

	lon1,lat1 = coordinates(geo_d01)
	lon2,lat2 = coordinates(geo_d02)
	a2[1].plot(lon1.+360,lat1)
	a2[1].plot(lon2.+360,lat2)
	
	a2[1].format(xlim=(269,296),ylim=(-1,16))

	ix = f2.add_axes([0.635,0.632,0.21,0.30])
	ix.format(
		xlim=(275,278),ylim=(8,10.5),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	ix = f2.add_axes([0.635,0.175,0.21,0.36])
	ix.plot(lon1.+360,lat1)
	ix.plot(lon2.+360,lat2)
	ix.format(
		xlim=(250,305),ylim=(-20,35),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)
	
	f2.savefig(plotsdir("01a-domain.png"),transparent=false,dpi=400)
	load(plotsdir("01a-domain.png"))
end

# ╔═╡ 7bb22689-ade3-4cc6-b0bf-a30d5a2fe9f7
md"
### D. Retrieving Land-Sea Mask
"

# ╔═╡ e6bf2538-b4d6-4e2c-b5f9-99464b6df643
# ╠═╡ disabled = true
#=╠═╡
lsd = getLandSea(geo_d02,path=datadir(),returnlsd=true,savelsd=true)
  ╠═╡ =#

# ╔═╡ 7df3ca55-c981-473a-9edb-7e0f7f957caf
# ╠═╡ disabled = true
#=╠═╡
lsd_sml = getLandSea(geo_fig2sml,path=datadir(),returnlsd=true,savelsd=true)
  ╠═╡ =#

# ╔═╡ c104d6ef-79ac-47dc-acc6-81242dc935af
# ╠═╡ disabled = true
#=╠═╡
lsd_big = getLandSea(geo_fig2big,path=datadir(),returnlsd=true,savelsd=true)
  ╠═╡ =#

# ╔═╡ b93e089e-40d0-440f-8ac9-d5f1152272e1
# ╠═╡ disabled = true
#=╠═╡
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=27/17,axwidth=5)

	a3[1].plot(lon1.+360,lat1)
	a3[1].plot(lon2.+360,lat2)
	a3[1].pcolormesh(lsd_big.lon[1:10:end],lsd_big.lat[1:10:end],lsd_big.z[1:10:end,1:10:end]')
	c = a3[1].pcolormesh(lsd.lon.+360,lsd.lat,lsd.z'./1000,levels=-6:6,cmap="delta",extend="both")
	
	a3[1].format(
		xlim=(269,296),xlocator=255:5:300,xminorlocator=240:2.5:315,
		ylim=(-1,16),ylocator=-15:5:30,yminorlocator=-20:2.5:35,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
		grid=true,gridcolor="w",#suptitle="Available Colombia Stations",
	)

	ix3 = f3.add_axes([0.635,0.632,0.21,0.30])
	ix3.format(
		xlim=(275,278),ylim=(8,10.5),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)
	ix3.pcolormesh(lsd_sml.lon,lsd_sml.lat,lsd_sml.z')

	ix3 = f3.add_axes([0.635,0.175,0.21,0.36])
	ix3.plot(lon1.+360,lat1)
	ix3.plot(lon2.+360,lat2)
	ix3.format(
		xlim=(250,305),ylim=(-20,35),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)
	ix3.pcolormesh(lsd_big.lon[1:10:end],lsd_big.lat[1:10:end],lsd_big.z[1:10:end,1:10:end]')
	ix3.pcolormesh(lsd.lon.+360,lsd.lat,lsd.z')

	f3.colorbar(c,label="Topographic Height/ km",length=0.8)
	f3.savefig(plotsdir("01a-domain.png"),transparent=false,dpi=400)
	load(plotsdir("01a-domain.png"))
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─f252a060-111a-4e86-85dd-46257c258b77
# ╟─ab5bab9c-fc33-4a6a-aa20-b9b1a44c4da8
# ╟─d5bc84d7-6516-4fad-a245-703e14e22f77
# ╟─7853fab4-b431-4872-9645-581e06728d19
# ╟─b2dc9507-7794-4445-b7f2-4351cba54af0
# ╟─12766ec3-8b45-49f5-908f-697322059df5
# ╟─9fba6b0f-5d17-4e3a-b9da-de815b12a6cf
# ╟─29fa97df-30dc-4258-b1e7-22662151919f
# ╠═a8f2c79d-d79d-4163-ad98-68578937629b
# ╠═dcc32199-77ac-4b8a-9fec-8d83cbf00cf2
# ╟─7bb22689-ade3-4cc6-b0bf-a30d5a2fe9f7
# ╟─e6bf2538-b4d6-4e2c-b5f9-99464b6df643
# ╟─7df3ca55-c981-473a-9edb-7e0f7f957caf
# ╟─c104d6ef-79ac-47dc-acc6-81242dc935af
# ╠═b93e089e-40d0-440f-8ac9-d5f1152272e1
