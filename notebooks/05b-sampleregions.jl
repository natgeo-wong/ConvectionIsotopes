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
	using NASAPrecipitation
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

# ╔═╡ e94be26b-b4d6-43ff-aa14-6c1c04c4a4c1
begin
	ds  = NCDataset(datadir("wrf3","grid.nc"))
	lon = ds["longitude"][:,:]
	lat = ds["latitude"][:,:]
	close(ds)
end

# ╔═╡ 04442b46-6181-4180-8c07-a73636f0e915
npd = IMERGFinalHH(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ 6705ac2e-27a7-4709-8a42-b1b4e7ff3939
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ dc10f521-772a-4833-99f8-b4b8eac22ca3
lsd = getLandSea(npd,geo)

# ╔═╡ f494290b-a2a3-47c9-80ed-44a1842e5c42
begin
	gH2O = zeros(length(lsd.lon),length(lsd.lat))
	for idt = Date(2019,8) : Day(1) : Date(2020,12,30)
		gds = read(npd,geo,idt)
		gH2O[:,:] += dropdims(sum(gds["precipitation"][:,:,:],dims=3),dims=3) * 1800
		close(gds)
	end
	gH2O ./= length(Date(2019,8) : Day(1) : Date(2020,12,30))
end

# ╔═╡ 95d17be8-7caa-4f1b-b64e-85af9afc935e
begin
	wds = NCDataset(datadir("wrf3","2D","RAINNC-daily-20190801_20201231.nc"))
	wH2O = dropdims(mean(wds["RAINNC"][:,:,vcat(1:152,155:end-1)],dims=3),dims=3)
	close(wds)
end

# ╔═╡ f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
md"
### C. Finding the Centroid and Tilt of the Grid
"

# ╔═╡ 7a2645c6-e78b-476e-862a-b13b0e58bace
function pnts2geo(pnts :: Vector{Point2{FT}}, ID, name) where FT <: Real
	for ipnt in 1 : length(pnts)
		ii = argmin((mod.(lon,360).-mod.(pnts[ipnt][1],360)).^2 .+ (lat.-pnts[ipnt][2]).^2)
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
	
		if isID("OTREC_wrf_$ID$(@sprintf("%02d",ipnt))",path=srcdir(),throw=false)
			rmID("OTREC_wrf_$ID$(@sprintf("%02d",ipnt))",path=srcdir())
		end
		geo_d01 = GeoRegion(
			ilon, ilat, rotation = -10,
			ID = "OTREC_wrf_$ID$(@sprintf("%02d",ipnt))", pID = "GLB",
			name = "WRF $(name) Sample Location $(@sprintf("%02d",ipnt))",
			path = srcdir(), save = true, checkshape = false
		)
	end
		
end

# ╔═╡ 51cc4262-48f1-461e-8618-a385ad3809c8
begin
	ITCZpnts = Point2.(collect(-90:0.5:-78),ones(25)*6.5)
	pnts2geo(ITCZpnts,"ITCZ","ITCZ")
end

# ╔═╡ 71744f0e-cc94-4312-89aa-efb7c937f48f
begin
	Pac2Atlpnts = Point2.(collect(-90:0.5:-78),collect(0:0.5:12).+2)
	pnts2geo(Pac2Atlpnts,"PAC2ATL","Pacific-Atlantic")
end

# ╔═╡ 91f67301-3472-4f69-966c-f2b63262b315
begin
	crossITCZpnts = Point2.(ones(25)*-88,collect(0:0.5:12).+1)
	pnts2geo(crossITCZpnts,"CrossITCZ","Cross ITCZ")
end

# ╔═╡ cb306a1f-affd-432c-8c46-977758531654
begin
	pplt.close(); fig,axs = pplt.subplots(axwidth=1.5,ncols=2)

	axs[1].pcolormesh(lsd.lon,lsd.lat,gH2O',levels=5:3:20,extend="both",cmap="Blues")
	axs[2].pcolormesh(lon,lat,wH2O,levels=5:3:20,extend="both",cmap="Blues")

	for ax in axs
		ax.plot(x,y,c="k",lw=0.5)
		for ipnt = 1 : length(ITCZpnts)
			geosample = GeoRegion("OTREC_wrf_ITCZ$(@sprintf("%02d",ipnt))",path=srcdir())
			ilon,ilat = coordinates(geosample)
			ax.plot(ilon,ilat,lw=0.5,c="pink",alpha=0.5.+(ipnt/length(ITCZpnts)/2))
		end
		for ipnt = 1 : length(Pac2Atlpnts)
			geosample = GeoRegion("OTREC_wrf_PAC2ATL$(@sprintf("%02d",ipnt))",path=srcdir())
			ilon,ilat = coordinates(geosample)
			ax.plot(ilon,ilat,lw=0.5,c="g",alpha=0.5.+(ipnt/length(Pac2Atlpnts)/2))
		end
		for ipnt = 1 : length(crossITCZpnts)
			geosample = GeoRegion("OTREC_wrf_CrossITCZ$(@sprintf("%02d",ipnt))",path=srcdir())
			ilon,ilat = coordinates(geosample)
			ax.plot(ilon,ilat,lw=0.5,c="b",alpha=0.5.+(ipnt/length(crossITCZpnts)/2))
		end
		for istn = 1 : 12
			geostn = GeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir())
			xc,yc = geostn.geometry.centroid
			ax.scatter(xc,yc,lw=0.5,c="red",s=2)
		end
		ax.format(xlim=(geo.W-0.2,geo.E+0.2),ylim=(geo.S-0.2,geo.N+0.2))
	end
	
	fig.savefig("test.png",transparent=false,dpi=250)
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
# ╠═e94be26b-b4d6-43ff-aa14-6c1c04c4a4c1
# ╠═04442b46-6181-4180-8c07-a73636f0e915
# ╠═6705ac2e-27a7-4709-8a42-b1b4e7ff3939
# ╠═dc10f521-772a-4833-99f8-b4b8eac22ca3
# ╠═f494290b-a2a3-47c9-80ed-44a1842e5c42
# ╠═95d17be8-7caa-4f1b-b64e-85af9afc935e
# ╟─f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
# ╠═7a2645c6-e78b-476e-862a-b13b0e58bace
# ╠═51cc4262-48f1-461e-8618-a385ad3809c8
# ╠═71744f0e-cc94-4312-89aa-efb7c937f48f
# ╟─91f67301-3472-4f69-966c-f2b63262b315
# ╠═cb306a1f-affd-432c-8c46-977758531654
