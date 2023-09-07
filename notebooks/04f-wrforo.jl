### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ c4157dda-6de3-4e17-8ba8-e2fc91e0141f
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ f49ad78f-f2b1-49e6-ba70-8e0e9c840fa1
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	using Interpolations
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ d56818f6-0893-11ed-2763-bdc0465266e9
md"
# 04f. WRF Topography/Orography
"

# ╔═╡ 3ddcbbf1-4772-4546-b811-15d69a91975c
begin
	fid = datadir("GLB-h.txt")
	if !isfile(fid)
		download("https://raw.githubusercontent.com/natgeo-wong/GeoPlottingData/main/coastline_resh.txt",fid)
	end
	coast = readdlm(fid,comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	icst = (x.>=270).&(x.<=285).&(y.>=0).&(y.<=15)
	x[.!icst] .= NaN
	y[.!icst] .= NaN
	md"Loading coastlines data ..."
end

# ╔═╡ a0a022ce-a1f6-4522-9eed-0a3f426506ed
geo = RectRegion("OTREC_lsd","GLB","Bigger OTREC Region",[15.5,-0.5,285.5,269.5],savegeo=false)

# ╔═╡ 0f18cfae-7d57-45b3-8b4c-2321d576474c
lsd = getLandSea(
	geo,path=datadir(),savelsd=true,returnlsd=true,
	smooth=true,σlon=30,σlat=30
)

# ╔═╡ aa1083c4-e014-487c-b5d3-b9707316145d
md"
### A. Binning WRF topography to ERA5 grid
"

# ╔═╡ d6338ae9-78c8-42b0-a025-16e7090c3432
begin
	ds  = NCDataset(datadir("wrf","wrfout","AUGp01","wrfout_d02_2019-08-01_00:00:00"))
	wln = ds["XLONG"][:,:,1] .+ 360; nwlon = size(wln,1)
	wlt = ds["XLAT"][:,:,1]; 		 nwlat = size(wlt,2)
	const oro = ds["HGT"][:,:,1] / 1000
	close(ds)
end

# ╔═╡ e0984fe0-b872-40e4-a301-4ffb5e4d582c
begin
	const ipnt_lon = zeros(Int,nwlon,nwlat)
	const ipnt_lat = zeros(Int,nwlon,nwlat)
	for ilat = 1 : nwlat, ilon = 1 : nwlon
	
		ipnt_lon[ilon,ilat] = argmin(abs.(wln[ilon,ilat].-lsd.lon))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlt[ilon,ilat].-lsd.lat))
	
	end
	md"Finding closest ERA5 points to each of the WRF points ..."
end

# ╔═╡ b40fb4a5-5e22-4396-b4e0-2835ed1007e4
begin
	const boro = zeros(length(lsd.lon),length(lsd.lat))
	md"Preallocating array for binned orography ..."
end

# ╔═╡ dc6a5a01-b03e-4aec-aeda-7bbe7e660040
for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
	ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
	iioro = @view oro[ind]
	if sum(iszero.(iioro)) > (sum(ind)/2)
		boro[ilon,ilat] = 0
	else
		boro[ilon,ilat] = mean(iioro[.!isnan.(iioro)])
	end
end

# ╔═╡ 98fb392c-0051-49f9-92b8-f19335e59364
lsd.z[lsd.lsm.<0.5] .= 0

# ╔═╡ 825997b1-445f-48fc-b344-87209400b060
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,nrows=2,axwidth=2)
	
	c1 = axs[2].pcolormesh(
		wln,wlt,oro,levels=vcat(-4:-1,-0.1,-0.01,0,0.01,0.1,1:4),
		extend="both",cmap="delta"
	)
	
	axs[1].pcolormesh(
		lsd.lon,lsd.lat,boro',
		levels=vcat(-4:-1,-0.1,-0.01,0,0.01,0.1,1:4),
		extend="both",cmap="delta"
	)
	
	axs[3].pcolormesh(
		lsd.lon,lsd.lat,lsd.z'/9.81/1000,
		levels=vcat(-4:-1,-0.1,-0.01,0,0.01,0.1,1:4),
		extend="both",cmap="delta"
	)
	
	c2 = axs[4].pcolormesh(
		lsd.lon,lsd.lat,boro' .- lsd.z'/9.81/1000,
		levels=vcat(-10:2:-2,-1,1,2:2:10)/10,
		extend="both",cmap="RdBu_r"
	)

	for ax in axs
		ax.plot(x,y,c="k",lw=1)
		ax.format(
			xlim=(270,285),xlocator=(270:3:285),xlabel=L"Longitude / $\degree$",
			ylim=(0,15),ylocator=0:3:15,ylabel=L"Latitude / $\degree$",
			xminorlocator=270:0.5:285,yminorlocator=0:0.5:15
		)
	end

	axs[1].format(ltitle="(a) WRF Topography (ETOPO)")
	axs[2].format(ltitle="(b) WRF Topography (Raw)")
	axs[3].format(ltitle="(c) ERA5 Topography")
	axs[4].format(ltitle="(d) WRF - ERA5 Topography")

	axs[2].colorbar(c1,label="km")
	axs[4].colorbar(c2,label="km")
	fig.savefig(plotsdir("04f-wrforo.png"),transparent=false,dpi=400)
	load(plotsdir("04f-wrforo.png"))
end

# ╔═╡ 1c0a436a-8d2c-47df-8f64-11e3e8ec1414
md"
### B. Interpolating Filtered Land-Sea Mask to WRF Grid
"

# ╔═╡ 5f13d09f-5ead-4139-b187-eadb52bdb8cd
begin
	itp = interpolate(lsd.lsm, BSpline(Linear()))
	sitp = scale(itp,lsd.lon[1]:(1/60):lsd.lon[end],lsd.lat[1]:(1/60):lsd.lat[end])
	md"Making Interpolations and scaling to grid"
end

# ╔═╡ d7937289-6ef3-4459-bdda-b3ed3e4483f4
begin
	wrf_flsm = zeros(nwlon,nwlat)
	for ilat = 1 : nwlat, ilon = 1 : nwlon
		wrf_flsm[ilon,ilat] = sitp(wln[ilon,ilat],wlt[ilon,ilat])
	end
end

# ╔═╡ 80339c17-f090-41e1-8b58-51590c89eee3
PCS = GeoRegion("OTREC_PAC")

# ╔═╡ 7a74e495-f5bf-41e0-bb00-d54f80d933cd
begin
	slon_PCS,slat_PCS = coordGeoRegion(PCS)
	md"Loading geographic shapes of the GeoRegions ..."
end

# ╔═╡ be932b20-c8db-4865-9c38-0f692838ff0f
ginfo = RegionGrid(PCS,wln,wlt)

# ╔═╡ 3fce6f0a-0b80-45ec-bfb1-0cbab2e11ee0
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,nrows=2,axwidth=2)

	llvls = -5:0.5:-1; llvls = llvls[llvls.<log10(0.5)]
	llvls = vcat(llvls,log10(0.25),log10(0.5))
	olvls = vcat(-4:-1,-0.1,-0.01,0,0.01,0.1,1:4)
	
	c2_1 = a2[1].pcolormesh(wln,wlt,oro,levels=olvls,cmap="delta",extend="both")
	c2_2 = a2[3].pcolormesh(wln,wlt,log10.(wrf_flsm),levels=llvls,cmap="delta",extend="both")
	a2[2].pcolormesh(wln,wlt,oro.*ginfo.mask,levels=olvls,cmap="delta",extend="both")
	a2[4].pcolormesh(wln,wlt,log10.(wrf_flsm.*ginfo.mask),levels=llvls,cmap="delta",extend="both")
	
	for ax in a2
		ax.plot(x,y,c="k",lw=1)
		ax.plot(slon_PCS.+360,slat_PCS,c="r")
		ax.format(
			xlim=(270,285),xlocator=(270:3:285),xlabel=L"Longitude / $\degree$",
			ylim=(0,15),ylocator=0:3:15,ylabel=L"Latitude / $\degree$",
			xminorlocator=270:0.5:285,yminorlocator=0:0.5:15
		)
	end

	a2[2].colorbar(c2_1)
	a2[4].colorbar(c2_2)

	a2[1].format(ltitle="(a) Topography / km")
	a2[2].format(ltitle="(b) Filtered Land-Sea Mask / 0-1")
	
	f2.savefig(plotsdir("04f-wrfflsm.png"),transparent=false,dpi=400)
	load(plotsdir("04f-wrfflsm.png"))
end

# ╔═╡ 21365092-62d7-4e18-90a0-00c3c9e0f5af
begin
	fnc = datadir("flsm","flsm_wrf.nc")
	if isfile(fnc)
		rm(fnc,force=true)
	end

	gds = NCDataset(fnc,"c")

	defDim(gds,"longitude",size(wln,1))
	defDim(gds,"latitude", size(wlt,2))

	nclon = defVar(gds,"longitude",Float64,("longitude","latitude",),attrib = Dict(
		"units"     => "degrees_east",
		"long_name" => "longitude",
	))

	nclat = defVar(gds,"latitude",Float64,("longitude","latitude",),attrib = Dict(
		"units"     => "degrees_north",
		"long_name" => "latitude",
	))

	ncvar = defVar(gds,"flsm",Float64,("longitude","latitude",),attrib = Dict(
		"long_name"     => "land_sea_mask_filtered",
		"full_name"     => "Filtered Land-Sea Mask",
		"units"         => "0-1",
	))

	nclon[:] = wln
	nclat[:] = wlt
	ncvar[:] = wrf_flsm
	
	close(gds)
end

# ╔═╡ Cell order:
# ╟─d56818f6-0893-11ed-2763-bdc0465266e9
# ╟─c4157dda-6de3-4e17-8ba8-e2fc91e0141f
# ╟─f49ad78f-f2b1-49e6-ba70-8e0e9c840fa1
# ╟─3ddcbbf1-4772-4546-b811-15d69a91975c
# ╟─a0a022ce-a1f6-4522-9eed-0a3f426506ed
# ╟─0f18cfae-7d57-45b3-8b4c-2321d576474c
# ╟─aa1083c4-e014-487c-b5d3-b9707316145d
# ╟─d6338ae9-78c8-42b0-a025-16e7090c3432
# ╟─e0984fe0-b872-40e4-a301-4ffb5e4d582c
# ╟─b40fb4a5-5e22-4396-b4e0-2835ed1007e4
# ╟─dc6a5a01-b03e-4aec-aeda-7bbe7e660040
# ╟─98fb392c-0051-49f9-92b8-f19335e59364
# ╟─825997b1-445f-48fc-b344-87209400b060
# ╟─1c0a436a-8d2c-47df-8f64-11e3e8ec1414
# ╟─5f13d09f-5ead-4139-b187-eadb52bdb8cd
# ╟─d7937289-6ef3-4459-bdda-b3ed3e4483f4
# ╠═80339c17-f090-41e1-8b58-51590c89eee3
# ╟─7a74e495-f5bf-41e0-bb00-d54f80d933cd
# ╟─be932b20-c8db-4865-9c38-0f692838ff0f
# ╠═3fce6f0a-0b80-45ec-bfb1-0cbab2e11ee0
# ╠═21365092-62d7-4e18-90a0-00c3c9e0f5af
