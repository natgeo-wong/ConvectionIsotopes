### A Pluto.jl notebook ###
# v0.20.5

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
	using Dates
	using DelimitedFiles
	using GeoRegions
	using RegionGrids
	using NCDatasets
	using PlutoUI
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01d. Creating GeoRegions for Stations

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
	pds = NCDataset(datadir("era5mo","OTREC_wrf_d02x0.25","tp","era5mo-OTREC_wrf_d02x0.25-tp-2019.nc"))
	pln = pds["longitude"][:]
	plt = pds["latitude"][:]
	H2O = dropdims(mean(pds["tp"][:,:,:],dims=3),dims=3) * 1000
	close(pds)
	pds = NCDataset(datadir("wrf3","2D","RAINNC-daily-20190801_20201231.nc"))
	pln = pds["longitude"][:,:]
	plt = pds["latitude"][:,:]
	H2O = dropdims(mean(pds["RAINNC"][:,:,vcat(1:152,155:end-1)],dims=3),dims=3)
	close(pds)
end

# ╔═╡ 7c9a5dd7-9fcd-4ef5-8694-2c6358c9b467
ndt = length(Date(2019,8,1) : Date(2020,1,1))

# ╔═╡ f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
md"
### C. Finding the Centroid and Tilt of the Grid
"

# ╔═╡ dccc851d-9639-4d0e-89c9-13fc03d962b9
setupGeoRegions(path=srcdir())

# ╔═╡ ed11fd66-cc0a-4623-a77c-bbbc7686f255
docreategeoregions = true

# ╔═╡ 7a2645c6-e78b-476e-862a-b13b0e58bace
if docreategeoregions
	ds  = NCDataset(datadir("wrf2","grid.nc"))
	lon = ds["longitude"][:,:]
	lat = ds["latitude"][:,:]
	close(ds)
	for istn = 1 : ndystn
	
		lon_stn = infody[istn,2]
		lat_stn = infody[istn,3]
	
		ii = argmin((mod.(lon,360).-lon_stn).^2 .+ (lat.-lat_stn).^2)
		ilon_stn = ii[1]
		ilat_stn = ii[2]
	
		ilon = vcat(
			lon[ilon_stn.+(-9:8),ilat_stn-9],
			lon[ilon_stn+9,ilat_stn.+(-9:8)],
			reverse(lon[ilon_stn.+(-8:9),ilat_stn+9]),
			reverse(lon[ilon_stn-9,ilat_stn.+(-8:9)])
		)
	
		ilat = vcat(
			lat[ilon_stn.+(-9:8),ilat_stn-9],
			lat[ilon_stn+9,ilat_stn.+(-9:8)],
			reverse(lat[ilon_stn.+(-8:9),ilat_stn+9]),
			reverse(lat[ilon_stn-9,ilat_stn.+(-8:9)])
		)
	
		if isID("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir(),throw=false)
			rmID("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir())
		end
		geo_d01 = GeoRegion(
			ilon, ilat, rotation = -10,
			ID = "OTREC_wrf_stn$(@sprintf("%02d",istn))", pID = "GLB",
			name = "WRF Station $(@sprintf("%02d",istn))",
			path = srcdir(), save = true
		)
	
		for iboxlat = 1 : 2, iboxlon = 1 : 2
	
			ibox = iboxlon + (iboxlat-1)*2
			ilon = vcat(
				lon[ilon_stn.+(-18:-1).+(iboxlon-1)*18,ilat_stn-18+(iboxlat-1)*18],
				lon[ilon_stn+(iboxlon-1)*18,ilat_stn.+(-18:-1).+(iboxlat-1)*18],
				reverse(lon[ilon_stn.+(-17:0).+(iboxlon-1)*18,ilat_stn+(iboxlat-1)*18]),
				reverse(lon[ilon_stn-18+(iboxlon-1)*18,ilat_stn.+(-17:0).+(iboxlat-1)*18])
			)
		
			ilat = vcat(
				lat[ilon_stn.+(-18:-1).+(iboxlon-1)*18,ilat_stn-18+(iboxlat-1)*18],
				lat[ilon_stn+(iboxlon-1)*18,ilat_stn.+(-18:-1).+(iboxlat-1)*18],
				reverse(lat[ilon_stn.+(-17:0).+(iboxlon-1)*18,ilat_stn+(iboxlat-1)*18]),
				reverse(lat[ilon_stn-18+(iboxlon-1)*18,ilat_stn.+(-17:0).+(iboxlat-1)*18])
			)
	
			if isID("OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox",path=srcdir(),throw=false)
				rmID("OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox",path=srcdir())
			end
			geo_d01 = GeoRegion(
				ilon, ilat, rotation = -10,
				ID = "OTREC_wrf_stn$(@sprintf("%02d",istn))_box$ibox", pID = "GLB",
				name = "WRF Station $(@sprintf("%02d",istn)) | Box $ibox",
				path = srcdir(), save = true
			)
	
			
		end
		
	end
end

# ╔═╡ 0f265e38-253f-49f2-9dc6-0d062ccb875d
if docreategeoregions
	for istn = 1 : ndystn
	
		lon_stn = infody[istn,2]
		lat_stn = infody[istn,3]
	
		ii = argmin((mod.(lon,360).-lon_stn).^2 .+ (lat.-lat_stn).^2)
		ilon_stn = ii[1]
		ilat_stn = ii[2]
	
		ilon = vcat(
			lon[ilon_stn.+(-18:17),ilat_stn-18],
			lon[ilon_stn+18,ilat_stn.+(-18:17)],
			reverse(lon[ilon_stn.+(-17:18),ilat_stn+18]),
			reverse(lon[ilon_stn-18,ilat_stn.+(-17:18)])
		)
	
		ilat = vcat(
			lat[ilon_stn.+(-18:17),ilat_stn-18],
			lat[ilon_stn+18,ilat_stn.+(-18:17)],
			reverse(lat[ilon_stn.+(-17:18),ilat_stn+18]),
			reverse(lat[ilon_stn-18,ilat_stn.+(-17:18)])
		)
	
		if isID("OTREC_wrf_stn$(@sprintf("%02d",istn))b",path=srcdir(),throw=false)
			rmID("OTREC_wrf_stn$(@sprintf("%02d",istn))b",path=srcdir())
		end
		geo_d01 = GeoRegion(
			ilon, ilat, rotation = -10,
			ID = "OTREC_wrf_stn$(@sprintf("%02d",istn))b", pID = "GLB",
			name = "WRF Station $(@sprintf("%02d",istn))b",
			path = srcdir(), save = true
		)
		
	end
end

# ╔═╡ 912a101b-7eb9-4322-a713-031aeffff20d
begin
	geo_stn =      GeoRegion("OTREC_wrf_stn04",path=srcdir())
	geo_stn_box1 = GeoRegion("OTREC_wrf_stn04_box1",path=srcdir())
	geo_stn_box2 = GeoRegion("OTREC_wrf_stn04_box2",path=srcdir())
	geo_stn_box3 = GeoRegion("OTREC_wrf_stn04_box3",path=srcdir())
	geo_stn_box4 = GeoRegion("OTREC_wrf_stn04_box4",path=srcdir())
	geo_stnb =     GeoRegion("OTREC_wrf_stn04b",path=srcdir())
end

# ╔═╡ 25bf5b00-0084-4241-909e-b48a7c9baaaa
begin
	lon1,lat1 = coordinates(geo_stn)
	lon2,lat2 = coordinates(geo_stn_box1)
	lon3,lat3 = coordinates(geo_stn_box2)
	lon4,lat4 = coordinates(geo_stn_box3)
	lon5,lat5 = coordinates(geo_stn_box4)
	lonb,latb = coordinates(geo_stnb)
end

# ╔═╡ 6175e264-2450-4074-b875-2cd8014a116f
ggrdstn = RegionGrid(geo_stn,Point2.(pln,plt),sigdigits=15)

# ╔═╡ d21bddbb-8f55-4a55-b3d7-01d98639eb6e
H2Ostn = extract(H2O,ggrdstn);

# ╔═╡ 4eb75767-dfb3-44e0-b4c0-3016b1d2a33d
ggrdbox = RegionGrid(geo_stn_box2,Point2.(pln,plt),sigdigits=15)

# ╔═╡ c883a08d-8c7d-42e6-a988-da8b354c8612
H2Obox = extract(H2O,ggrdbox);

# ╔═╡ 52521082-b677-44b7-8ee9-2cc56a697c90
it = 30

# ╔═╡ cb306a1f-affd-432c-8c46-977758531654
begin
	pplt.close(); fig,axs = pplt.subplots(axwidth=1.5,ncols=3)

	c = axs[1].pcolormesh(pln,plt,H2O,levels=2:2:20,extend="both",cmap="Blues")
	axs[2].pcolormesh(ggrdstn.lon,ggrdstn.lat,H2Ostn,levels=2:2:20,extend="both",cmap="Blues")
	axs[3].pcolormesh(ggrdbox.lon,ggrdbox.lat,H2Obox,levels=2:2:20,extend="both",cmap="Blues")

	xc,yc = geo_stn.geometry.centroid
	for ax in axs
		ax.plot(lon1,lat1,c="k")
		ax.plot(lon2,lat2,c="k")
		ax.plot(lon3,lat3,c="k")
		ax.plot(lon4,lat4,c="k")
		ax.plot(lon5,lat5,c="k")
		# ax.plot(lonb,latb,lw=5,c="k")
		ax.plot(x,y,c="k",lw=0.5)
		ax.scatter(xc,yc,c=:r)
		ax.format(
			xlim=(geo_stnb.W-0.2,geo_stnb.E+0.2),xlabel=L"Longitude / $\degree$",
			ylim=(geo_stnb.S-0.2,geo_stnb.N+0.2),ylabel=L"Latitude / $\degree$"
		)
	end

	axs[1].format(ltitle="(a) Entire WRF grid")
	axs[2].format(ltitle="(b) OTREC_wrf_stn04")
	axs[3].format(ltitle="(c) OTREC_wrf_stn04_box2")

	fig.colorbar(c,label=L"Precipitation Rate / mm day$^{-1}$")
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ b0086e10-a75e-4256-b847-6d083b44ee9d
begin
	pplt.close(); f2,a2 = pplt.subplots(axwidth=1)
	
	a2[1].pcolormesh(
		ggrdstn.X./ 1e3,ggrdstn.Y./ 1e3,H2Ostn,
		levels=2:2:20,extend="both",cmap="Blues"
	)

	r1clon,r1clat = derotatecoordinates(x,y,geo_stn) ./ 1e3
	rlon,rlat = coordinates(geo_stn,derotate=true) ./ 1e3
	a2[1].plot(r1clon,r1clat)
	a2[1].plot(rlon,rlat)
	a2[1].format(xlim=(-1,1).*1e2,ylim=(-1,1).*1e2)
	
	f2.savefig("test.png",transparent=false,dpi=150)
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
# ╠═7c9a5dd7-9fcd-4ef5-8694-2c6358c9b467
# ╟─f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
# ╠═dccc851d-9639-4d0e-89c9-13fc03d962b9
# ╠═ed11fd66-cc0a-4623-a77c-bbbc7686f255
# ╟─7a2645c6-e78b-476e-862a-b13b0e58bace
# ╟─0f265e38-253f-49f2-9dc6-0d062ccb875d
# ╠═912a101b-7eb9-4322-a713-031aeffff20d
# ╠═25bf5b00-0084-4241-909e-b48a7c9baaaa
# ╠═6175e264-2450-4074-b875-2cd8014a116f
# ╠═d21bddbb-8f55-4a55-b3d7-01d98639eb6e
# ╠═4eb75767-dfb3-44e0-b4c0-3016b1d2a33d
# ╠═c883a08d-8c7d-42e6-a988-da8b354c8612
# ╠═52521082-b677-44b7-8ee9-2cc56a697c90
# ╠═cb306a1f-affd-432c-8c46-977758531654
# ╠═b0086e10-a75e-4256-b847-6d083b44ee9d
