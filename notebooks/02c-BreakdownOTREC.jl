### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 1743ac77-eddb-4d08-8841-eaf63a451b8f
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ f6017eb1-5b2f-4113-ab6b-91e973659aa2
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using PlutoUI
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("vertprofile_region.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ 8585baa8-f995-11ec-2f6b-052638e59e43
md"
# 02c. OTREC Domain Breakdown
"

# ╔═╡ 9974fed3-9299-421f-8bc1-05b6a3f738b0
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ ec7c7485-09f3-410a-b5f2-a7768e3b9ba6
md"
### A. Loading Dataset Information ...
"

# ╔═╡ 1431f296-2715-47ad-af22-0342d5654f4a
e5ds = ERA5Monthly(dtbeg=Date(2013,1,1),dtend=Date(2021,12,31),eroot=datadir())

# ╔═╡ bb48e553-6c6f-4a68-b0a2-c81d5b8d51d3
ereg = ERA5Region(GeoRegion("OTREC"))

# ╔═╡ 657ca136-c1ce-4c3a-a7f0-aa0dce82348b
evar = SingleVariable("p_wwgt")

# ╔═╡ 8d79e82b-69dc-4a77-8c6b-d7c79c6d233f
begin
	ds  = read(e5ds,evar,ereg,Date(2013))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	var = nomissing(ds[evar.varID][:,:,1],0) / 100
	cnt = Int.(.!iszero.(var))
	close(ds)
end

# ╔═╡ 2d3ec34f-d4a2-4700-8561-6c89a880628a
for dt in Date(2014) : Year(1) : Date(2021)
	ids = read(e5ds,evar,ereg,dt)
	var[:,:] += nomissing(ids[evar.varID][:,:,1],0) / 100
	cnt[:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:,:,1],0)))
	close(ids)	
end

# ╔═╡ 1e368291-a98f-4719-b598-f0c355f09e16
begin
	wwgt_pre = var ./ cnt
	
	lds = NCDataset(datadir("flsm","flsm-$(ereg.geoID).nc"))
	lsm = lds["flsm"][:]; nlon = size(lsm,1); nlat = size(lsm,2)
	close(lds)
	
	for ilat = 1 : nlat, ilon = 1 : nlon
		if lsm[ilon,ilat] > 0.9
			wwgt_pre[ilon,ilat] = NaN
		end
	end
	md"Filtering out data based on the filtered land-sea mask ..."
end

# ╔═╡ 00b914a4-5c64-4b5f-b770-ead2608cd9d0
md"
### B. Loading NASA Precipitation Datasets
"

# ╔═╡ c4388c8c-784d-4185-a7f2-dea7c0198140
npd = IMERGMonthly(dtbeg=Date(2013),dtend=Date(2020),sroot=datadir())

# ╔═╡ d59a7856-3f0c-461b-9fdf-3eee8e03ac7b
begin
	pds = read(npd,ereg.geo,Date(2013))
	pln = pds["longitude"][:]
	plt = pds["latitude"][:]
	prc = pds["prcp_rate"][:] / 8 * 86400
	close(pds)
end

# ╔═╡ 098646e1-831f-41fb-99e6-7d798ff9e5ec
for dt in Date(2014) : Year(1) : Date(2020)
	ids = read(npd,ereg.geo,dt)
	prc[:,:,:] += ids["prcp_rate"][:] / 8 * 86400
	close(ids)	
end

# ╔═╡ fd15563d-5ac5-49c3-9f5a-a97e386cc1c5
md"
### C. Loading and Calculating Regional-Mean Vertical Profile
"

# ╔═╡ 3f588b23-571d-4e92-b107-4b130c1fb63e
geo_btm = RectRegion(
	"OTREC_BTM","OTREC","OTREC Bottom Heavy",
	[3.5,1,278.5,276],savegeo=false
)

# ╔═╡ cb43b21b-7e67-426e-be12-ef1e32394ab5
geo_top = RectRegion(
	"OTREC_TOP","OTREC","OTREC Top Heavy",
	[7.8,7.2,282.1,281.5],savegeo=false
)

# ╔═╡ c20aff23-5b3f-4d5f-b347-83f79797cf54
wprf_μ_BTM,wprf_σ_BTM = vertprofile(e5ds,ERA5Region(geo_btm))

# ╔═╡ 2c62e5b4-82bc-4028-9417-f824a3e6880a
begin
	evar_wp = SingleVariable("p_wwgt")
	extract(e5ds,evar_wp,ERA5Region(geo_top))
end

# ╔═╡ 76d50655-6ecc-4246-b189-9ec7ccb06a5b
begin
	p = era5Pressures(); p = p[p.>=10]
	for ip in p
		evar_w = PressureVariable("w",hPa=ip)
		extract(e5ds,evar_w,ERA5Region(geo_top))
	end
end

# ╔═╡ fe3cad6f-a9b8-4f05-8057-d4279a17a6d1
wprf_μ_TOP,wprf_σ_TOP = vertprofile(e5ds,ERA5Region(geo_top))

# ╔═╡ d75ac010-7be3-40fe-9982-d0037f9e155d
begin
	blon_btm,blat_btm = coordGeoRegion(geo_btm)
	blon_top,blat_top = coordGeoRegion(geo_top)
end

# ╔═╡ c9ffefa3-e51a-406a-ab98-b48bc50e3bc5
begin
	pplt.close();
	asp = (ereg.geo.E-ereg.geo.W)/(ereg.geo.N-ereg.geo.S)
	f1,a1 = pplt.subplots([1,1,1,1,2],axwidth=2.5,aspect=asp,sharex=0,sharey=0)
	
	c = a1[1].pcolormesh(
		lon,lat,wwgt_pre',levels=250:50:850,
		extend="both",cmap="delta"
	)
	a1[1].contour(pln.+360,plt,prc[:,:,1]',levels=[5,10,20],c="k")
	a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].plot(blon_btm,blat_btm)
	a1[1].plot(blon_top,blat_top)
	a1[1].colorbar(c,loc="l",label=L"$p_w$ / hPa",length=0.8)
	a1[1].format(
		xlim=(minimum(lon),maximum(lon)),xlabel=L"Longitude / $\degree$",
		ylim=(minimum(lat),maximum(lat)),ylabel=L"Latitude / $\degree$",
		suptitle="W-weighted Mean Pressure (2013-2021, Jan)",
		ytickloc="right"
	)

	a1[2].plot(wprf_μ_BTM,p)
	a1[2].fill_betweenx(p,wprf_μ_BTM-wprf_σ_BTM,wprf_μ_BTM+wprf_σ_BTM,alpha=0.3)
	a1[2].plot(wprf_μ_TOP,p)
	a1[2].fill_betweenx(p,wprf_μ_TOP-wprf_σ_TOP,wprf_μ_TOP+wprf_σ_TOP,alpha=0.3)
	a1[2].format(
		ylim=(1000,100),ylabel="Pressure / hPa",
		xlim=(0.05,-0.15),xlabel=L"$\omega$ / Pa s$^{-1}$"
	)

	f1.savefig(plotsdir("02c-BreakdownOTREC.png"),transparent=false,dpi=400)
	load(plotsdir("02c-BreakdownOTREC.png"))
end

# ╔═╡ Cell order:
# ╟─8585baa8-f995-11ec-2f6b-052638e59e43
# ╟─1743ac77-eddb-4d08-8841-eaf63a451b8f
# ╟─f6017eb1-5b2f-4113-ab6b-91e973659aa2
# ╟─9974fed3-9299-421f-8bc1-05b6a3f738b0
# ╟─ec7c7485-09f3-410a-b5f2-a7768e3b9ba6
# ╟─1431f296-2715-47ad-af22-0342d5654f4a
# ╟─bb48e553-6c6f-4a68-b0a2-c81d5b8d51d3
# ╟─657ca136-c1ce-4c3a-a7f0-aa0dce82348b
# ╟─8d79e82b-69dc-4a77-8c6b-d7c79c6d233f
# ╟─2d3ec34f-d4a2-4700-8561-6c89a880628a
# ╟─1e368291-a98f-4719-b598-f0c355f09e16
# ╟─00b914a4-5c64-4b5f-b770-ead2608cd9d0
# ╟─c4388c8c-784d-4185-a7f2-dea7c0198140
# ╟─d59a7856-3f0c-461b-9fdf-3eee8e03ac7b
# ╟─098646e1-831f-41fb-99e6-7d798ff9e5ec
# ╟─fd15563d-5ac5-49c3-9f5a-a97e386cc1c5
# ╠═3f588b23-571d-4e92-b107-4b130c1fb63e
# ╠═cb43b21b-7e67-426e-be12-ef1e32394ab5
# ╠═c20aff23-5b3f-4d5f-b347-83f79797cf54
# ╠═2c62e5b4-82bc-4028-9417-f824a3e6880a
# ╠═76d50655-6ecc-4246-b189-9ec7ccb06a5b
# ╠═fe3cad6f-a9b8-4f05-8057-d4279a17a6d1
# ╠═d75ac010-7be3-40fe-9982-d0037f9e155d
# ╟─c9ffefa3-e51a-406a-ab98-b48bc50e3bc5
