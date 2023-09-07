### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ e32a00ee-5f32-47a1-a983-91fb77bc5d18
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using ERA5Reanalysis
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 02a. Exploring W-Weighted Mean Pressure
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 75925166-d1e8-450e-bae6-29affd25d635
md"
### A. Defining Datasets, Variables and Regions
"

# ╔═╡ d3f4a6a9-0915-470b-83af-e05e00dff5d2
npd = IMERGMonthly(start=Date(2013),stop=Date(2020),path=datadir())

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 9473b30b-787a-495b-bdfe-b386c9995761
evar = SingleVariable("p_wwgt")

# ╔═╡ 444e1d6a-6edb-4511-a9aa-e251c5b8013b
pvar = SingleVariable("tp")

# ╔═╡ c7d0ae59-8d78-4946-8348-9fa570197b0b
geo = GeoRegion("OTREC")

# ╔═╡ 3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
ereg = ERA5Region(geo)

# ╔═╡ fecbd37b-f74a-4d70-9bc7-540adec0c23d
mo = 6

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
begin
	ds  = read(e5ds,evar,ereg,Date(2013))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	var = nomissing(ds[evar.varID][:],0) / 100
	cnt = Int.(.!iszero.(var))
	close(ds)
end

# ╔═╡ b2833181-0706-4f0b-be14-85ef6b682945
for dt in Date(2014) : Year(1) : Date(2021)
	ids = read(e5ds,evar,ereg,dt)
	var[:,:,:] += nomissing(ids[evar.varID][:],0) / 100
	cnt[:,:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:],0)))
	close(ids)	
end

# ╔═╡ 5fc05bef-09ab-4b71-b3b8-52824a4132e4
begin
	pds = read(npd,geo,Date(2013))
	pln = pds["longitude"][:]
	plt = pds["latitude"][:]
	prc = pds["precipitation"][:] / 8 * 86400
	close(pds)
end

# ╔═╡ cdb15f83-4235-447c-95c1-8753d46ece27
for dt in Date(2014) : Year(1) : Date(2020)
	ids = read(npd,geo,dt)
	prc[:,:,:] += ids["precipitation"][:] / 8 * 86400
	close(ids)	
end

# ╔═╡ 09a04b03-d779-4277-85df-77e5448a533b
begin
	wwgt_pre = var ./ cnt
	
	lds = NCDataset(datadir("flsm","flsm-$(ereg.geoID).nc"))
	lsm = lds["flsm"][:]; nlon = size(lsm,1); nlat = size(lsm,2)
	close(lds)
	
	for ilat = 1 : nlat, ilon = 1 : nlon
		if lsm[ilon,ilat] > 0.9
			wwgt_pre[ilon,ilat,:] .= NaN
		end
	end
	md"Filtering out data based on the filtered land-sea mask ..."
end

# ╔═╡ 14663694-56f8-4651-ab52-8560374be699
begin
	pplt.close(); asp = (ereg.geo.E-ereg.geo.W)/(ereg.geo.N-ereg.geo.S)
	if asp > 2
		f1,a1 = pplt.subplots(nrows=4,axwidth=asp*1.2,aspect=asp)
	else
		f1,a1 = pplt.subplots(nrows=2,ncols=2,axwidth=2.5,aspect=asp)
	end
	
	c = a1[1].pcolormesh(
		lon,lat,wwgt_pre[:,:,1]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# a1[1].contour(pln.+360,plt,prc[:,:,1]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].format(ultitle="(a) Jan")

	a1[2].pcolormesh(
		lon,lat,wwgt_pre[:,:,4]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# a1[2].contour(pln.+360,plt,prc[:,:,4]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	a1[2].plot(x,y,c="k",lw=0.5)
	a1[2].format(ultitle="(b) Apr")

	a1[3].pcolormesh(
		lon,lat,wwgt_pre[:,:,7]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# a1[3].contour(pln.+360,plt,prc[:,:,7]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	a1[3].plot(x,y,c="k",lw=0.5)
	a1[3].format(ultitle="(c) Jul")

	a1[4].pcolormesh(
		lon,lat,wwgt_pre[:,:,10]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# a1[4].contour(pln.+360,plt,prc[:,:,10]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	a1[4].plot(x,y,c="k",lw=0.5)
	a1[4].format(ultitle="(d) Oct")

	for ax in a1
		ax.format(
			xlim=(minimum(lon),maximum(lon)),xlabel=L"Longitude / $\degree$",
			ylim=(minimum(lat),maximum(lat)),ylabel=L"Latitude / $\degree$",
			suptitle="W-weighted Mean Pressure (2013-2020)"
		)
	end

	f1.colorbar(c,label=L"$p_w$ / hPa",length=0.5)
	f1.savefig(plotsdir("02a-wwgt_pre-$(ereg.geoID).png"),transparent=false,dpi=400)
	load(plotsdir("02a-wwgt_pre-$(ereg.geoID).png"))
end

# ╔═╡ c2904ec4-24af-4771-bd94-31c020078fcf
md"
### C. W-Weighted Mean Pressure in Different Regions
"

# ╔═╡ 88c6924d-1bbb-4769-a6e0-2d8f4fe1456c
ereg_WPC = ERA5Region(GeoRegion("CIS_WPC"))

# ╔═╡ 25efdeaf-e22f-495a-843f-e05201f9dfa3
ereg_CST = ERA5Region(GeoRegion("OTREC_CST"))

# ╔═╡ cd0d4b80-6b22-4831-8d37-1a9489d3a9df
ereg_PAC = ERA5Region(GeoRegion("OTREC_PAC"))

# ╔═╡ c07fc5fe-9c39-4813-84ab-77016de261a9
begin
	ds_WPC  = read(e5ds,evar,ereg_WPC,Date(2013))
	lon_WPC = ds_WPC["longitude"][:]
	lat_WPC = ds_WPC["latitude"][:]
	var_WPC = nomissing(ds_WPC[evar.varID][:],0) / 100
	cnt_WPC = Int.(.!iszero.(var_WPC))
	close(ds_WPC)
	ds_PAC  = read(e5ds,evar,ereg_PAC,Date(2013))
	lon_PAC = ds_PAC["longitude"][:]
	lat_PAC = ds_PAC["latitude"][:]
	var_PAC = nomissing(ds_PAC[evar.varID][:],0) / 100
	cnt_PAC = Int.(.!iszero.(var_PAC))
	close(ds_PAC)
	ds_CST  = read(e5ds,evar,ereg_CST,Date(2013))
	lon_CST = ds_CST["longitude"][:]
	lat_CST = ds_CST["latitude"][:]
	var_CST = nomissing(ds_CST[evar.varID][:],0) / 100
	cnt_CST = Int.(.!iszero.(var_CST))
	close(ds_CST)
end

# ╔═╡ 08f9dfd2-3d70-400f-8be4-dcc92348e6ea
for dt in Date(2014) : Year(1) : Date(2021)
	ids = read(e5ds,evar,ereg_WPC,dt)
	var_WPC[:,:,:] += nomissing(ids[evar.varID][:],0) / 100
	cnt_WPC[:,:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:],0)))
	close(ids)
	ids = read(e5ds,evar,ereg_PAC,dt)
	var_PAC[:,:,:] += nomissing(ids[evar.varID][:],0) / 100
	cnt_PAC[:,:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:],0)))
	close(ids)
	ids = read(e5ds,evar,ereg_CST,dt)
	var_CST[:,:,:] += nomissing(ids[evar.varID][:],0) / 100
	cnt_CST[:,:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:],0)))
	close(ids)
end

# ╔═╡ ea75bee7-74ac-42ef-83fe-35544b02eaaf
begin
	wwgt_pre_WPC = var_WPC ./ cnt_WPC
	wwgt_pre_PAC = var_PAC ./ cnt_PAC
	wwgt_pre_CST = var_CST ./ cnt_CST

	lds_WPC  = NCDataset(datadir("flsm","flsm-$(ereg_WPC.geoID).nc"))
	lsm_WPC  = lds_WPC["flsm"][:]
	nlon_WPC = size(lsm_WPC,1)
	nlat_WPC = size(lsm_WPC,2)
	close(lds_WPC)
	for ilat = 1 : nlat_WPC, ilon = 1 : nlon_WPC
		if (lsm_WPC[ilon,ilat] > 0.99) && (lsm_WPC[ilon,ilat] < 0.2)
			wwgt_pre_WPC[ilon,ilat,:] .= NaN
		end
	end

	lds_PAC  = NCDataset(datadir("flsm","flsm-$(ereg_PAC.geoID).nc"))
	lsm_PAC  = lds_PAC["flsm"][:]
	nlon_PAC = size(lsm_PAC,1)
	nlat_PAC = size(lsm_PAC,2)
	close(lds_PAC)
	for ilat = 1 : nlat_PAC, ilon = 1 : nlon_PAC
		if (lsm_PAC[ilon,ilat] > 0.99) && (lsm_PAC[ilon,ilat] < 0.2)
			wwgt_pre_PAC[ilon,ilat,:] .= NaN
		end
	end

	lds_CST  = NCDataset(datadir("flsm","flsm-$(ereg_CST.geoID).nc"))
	lsm_CST  = lds_CST["flsm"][:]
	nlon_CST = size(lsm_CST,1)
	nlat_CST = size(lsm_CST,2)
	close(lds_CST)
	for ilat = 1 : nlat_CST, ilon = 1 : nlon_CST
		if (lsm_CST[ilon,ilat] > 0.99) && (lsm_CST[ilon,ilat] < 0.2)
			wwgt_pre_CST[ilon,ilat,:] .= NaN
		end
	end
	
	md"Filtering out data based on the filtered land-sea mask ..."
end

# ╔═╡ 099fef7c-59bf-45d1-92ef-96959229b93f
begin
	wwgt_pre_WPC_μ = zeros(12)
	wwgt_pre_PAC_μ = zeros(12)
	wwgt_pre_CST_μ = zeros(12)

	for it = 1 : 12
		wwgt = @view wwgt_pre_WPC[:,:,it]
		wwgt_pre_WPC_μ[it] = mean(wwgt[.!isnan.(wwgt)])
		wwgt = @view wwgt_pre_PAC[:,:,it]
		wwgt_pre_PAC_μ[it] = mean(wwgt[.!isnan.(wwgt)])
		wwgt = @view wwgt_pre_CST[:,:,it]
		wwgt_pre_CST_μ[it] = mean(wwgt[.!isnan.(wwgt)])
	end
	
	md"Finding the monthly means ..."
end

# ╔═╡ 7a28d1e2-d6d1-4447-abf1-c5bbd904fe5c
begin
	pplt.close(); f2,a2 = pplt.subplots([1,1,1,1,1,1,1,1,2],aspect=2,axwidth=3,sharex=0)
	
	a2[1].plot(1:12,wwgt_pre_WPC_μ)
	a2[1].plot(1:12,wwgt_pre_PAC_μ)
	a2[1].plot(1:12,wwgt_pre_CST_μ)
	a2[1].format(xlocator=1:12,ylabel=L"$p_w$ / hPa",xlabel="Month")

	a2[2].scatter(0,mean(wwgt_pre_WPC_μ),label="WPC",legend="r",legend_kw=Dict("frame"=>false,"ncol"=>1))
	a2[2].scatter(0,mean(wwgt_pre_PAC_μ),label="PAC",legend="r")
	a2[2].scatter(0,mean(wwgt_pre_CST_μ),label="CST",legend="r")
	a2[2].format(xlocator=0:0,xticklabels=[L"$\mu$"])
	
	f2.savefig(plotsdir("02a-ExploreWWgtPre.png"),transparent=false,dpi=150)
	load(plotsdir("02a-ExploreWWgtPre.png"))
end

# ╔═╡ eebc6c53-2aa2-412a-8724-780d65e05ef9
md"
### D. Tropical and OTREC Regions
"

# ╔═╡ e79ca9bf-790b-4744-87a6-95fc86d1f040
begin
	tds = read(e5ds,evar,ERA5Region(GeoRegion("TRP")),Date(2013))
	ipd = read(e5ds,pvar,ERA5Region(GeoRegion("TRP")),Date(2013))
	tln = tds["longitude"][:]
	tlt = tds["latitude"][:]
	tvr = nomissing(tds[evar.varID][:],0) / 100
	prcpii = nomissing(ipd[pvar.varID][:],0)
	ip = (prcpii*30*1000 .> 6) .& (prcpii/30*1000 .< 5)
	tvr[ip] .= 0
	tcn = Int.(.!iszero.(tvr))
	close(tds)
	close(ipd)
end

# ╔═╡ f2c2c5f6-2a4e-41d1-8a6b-b8884932fbc7
for dt in Date(2014) : Year(1) : Date(2021)
	ids = read(e5ds,evar,ERA5Region(GeoRegion("TRP")),dt)
	pds = read(e5ds,pvar,ERA5Region(GeoRegion("TRP")),dt)
	tvrii = nomissing(ids[evar.varID][:],0)
	prcpii = nomissing(pds[pvar.varID][:],0)
	ii = (prcpii*30*1000 .> 6) .& (prcpii*30*1000 .< 5)
	tvrii[ii] .= 0
	tvr[:,:,:] += tvrii / 100
	tcn[:,:,:] += Int.(.!iszero.(tvrii))
	close(ids)
	close(pds)
end

# ╔═╡ 0f3192f2-35aa-484e-bd66-d3eb9d714db1
begin
	wwgt_pre_TRP = tvr ./ tcn
	
	tfds = NCDataset(datadir("flsm","flsm-$(ERA5Region(GeoRegion("TRP")).geoID).nc"))
	tlsm = tfds["flsm"][:]; ntlon = size(tlsm,1); ntlat = size(tlsm,2)
	close(tfds)
	
	for ilat = 1 : ntlat, ilon = 1 : ntlon
		if tlsm[ilon,ilat] > 0.9
			wwgt_pre_TRP[ilon,ilat,:] .= NaN
		end
	end
	md"Filtering out data based on the filtered land-sea mask ..."
end

# ╔═╡ ad184f2f-fc45-4376-b700-ebbc9f2910ef
blon,blat = coordGeoRegion(GeoRegion("OTREC"))

# ╔═╡ f562b4cc-93a2-41a5-85d7-ec47c390751e
begin
	arr = [[1,1,1,1,1,1,5],[2,2,2,2,2,2,6],[3,3,3,3,3,3,7],[4,4,4,4,4,4,8]]
	pplt.close();
	fig,axs = pplt.subplots(arr,aspect=6,axwidth=7,sharey=1,hspace=1.5)
	
	axs[1].pcolormesh(
		tln,tlt,wwgt_pre_TRP[:,:,1]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	axs[1].plot(blon.+360,blat,c="r")
	axs[1].format(ultitle="(a) Jan")

	axs[2].pcolormesh(
		tln,tlt,wwgt_pre_TRP[:,:,4]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	axs[2].plot(blon.+360,blat,c="r")
	axs[2].format(ultitle="(b) Apr")

	axs[3].pcolormesh(
		tln,tlt,wwgt_pre_TRP[:,:,7]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	axs[3].plot(blon.+360,blat,c="r")
	axs[3].format(ultitle="(c) Jul")

	axs[4].pcolormesh(
		tln,tlt,wwgt_pre_TRP[:,:,10]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	axs[4].plot(blon.+360,blat,c="r")
	axs[4].format(ultitle="(d) Oct")
	
	cxs = axs[5].pcolormesh(
		lon,lat,wwgt_pre[:,:,1]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# axs[5].contour(pln.+360,plt,prc[:,:,1]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	axs[5].format(ultitle="(e) Jan")

	axs[6].pcolormesh(
		lon,lat,wwgt_pre[:,:,4]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# axs[6].contour(pln.+360,plt,prc[:,:,4]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	axs[6].format(ultitle="(f) Apr")

	axs[7].pcolormesh(
		lon,lat,wwgt_pre[:,:,7]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# axs[7].contour(pln.+360,plt,prc[:,:,7]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	axs[7].format(ultitle="(g) Jul")

	axs[8].pcolormesh(
		lon,lat,wwgt_pre[:,:,10]',levels=250:50:850,
		extend="both",cmap="delta"
	)
	# axs[8].contour(pln.+360,plt,prc[:,:,10]',levels=5:5:30,cmap="blues",extend="both",cmap_kw=Dict("left"=>0.3))
	axs[8].format(ultitle="(h) Oct")

	for iax in 5 : 8
		axs[iax].format(
			xlim=(minimum(lon),maximum(lon)),xlabel=L"Longitude / $\degree$",
			ylim=(minimum(lat),maximum(lat)),ylabel=L"Latitude / $\degree$",
			suptitle="W-weighted Mean Pressure"
		)
	end

	for iax in 1 : 4
		axs[iax].format(
			xlim=(0,360),xlabel=L"Longitude / $\degree$",
			ylim=(-30,30),ylabel=L"Latitude / $\degree$",
			xlocator=30:30:330,
			suptitle="W-weighted Mean Pressure"
		)
	end

	for ax in axs
		ax.plot(x,y,c="k",lw=0.5)
	end

	fig.colorbar(cxs,label=L"$p_w$ / hPa",length=0.5)
	fig.savefig(plotsdir("02a-wwgt_pre.png"),transparent=false,dpi=400)
	load(plotsdir("02a-wwgt_pre.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─d3f4a6a9-0915-470b-83af-e05e00dff5d2
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─9473b30b-787a-495b-bdfe-b386c9995761
# ╟─444e1d6a-6edb-4511-a9aa-e251c5b8013b
# ╟─c7d0ae59-8d78-4946-8348-9fa570197b0b
# ╟─3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
# ╠═fecbd37b-f74a-4d70-9bc7-540adec0c23d
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╟─acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
# ╟─b2833181-0706-4f0b-be14-85ef6b682945
# ╟─5fc05bef-09ab-4b71-b3b8-52824a4132e4
# ╟─cdb15f83-4235-447c-95c1-8753d46ece27
# ╟─09a04b03-d779-4277-85df-77e5448a533b
# ╟─14663694-56f8-4651-ab52-8560374be699
# ╟─c2904ec4-24af-4771-bd94-31c020078fcf
# ╟─88c6924d-1bbb-4769-a6e0-2d8f4fe1456c
# ╟─25efdeaf-e22f-495a-843f-e05201f9dfa3
# ╟─cd0d4b80-6b22-4831-8d37-1a9489d3a9df
# ╟─c07fc5fe-9c39-4813-84ab-77016de261a9
# ╟─08f9dfd2-3d70-400f-8be4-dcc92348e6ea
# ╟─ea75bee7-74ac-42ef-83fe-35544b02eaaf
# ╟─099fef7c-59bf-45d1-92ef-96959229b93f
# ╟─7a28d1e2-d6d1-4447-abf1-c5bbd904fe5c
# ╟─eebc6c53-2aa2-412a-8724-780d65e05ef9
# ╠═e79ca9bf-790b-4744-87a6-95fc86d1f040
# ╠═f2c2c5f6-2a4e-41d1-8a6b-b8884932fbc7
# ╠═0f3192f2-35aa-484e-bd66-d3eb9d714db1
# ╟─ad184f2f-fc45-4376-b700-ebbc9f2910ef
# ╟─f562b4cc-93a2-41a5-85d7-ec47c390751e
