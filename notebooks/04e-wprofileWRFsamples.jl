### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 3c4b037f-c854-4148-ae6e-fff50bbb51cf
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using NumericalIntegration
	using Statistics
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
md"
# 04e. Isotope vs W-Weighted Pressure
"

# ╔═╡ c853d30e-5a51-47f4-a6be-88ab857253c1
md"
### A. Loading some Sample Data
"

# ╔═╡ 3922aee7-0e55-4074-af98-9f09919df7df
begin
	geo_1  = GeoRegion("OTREC_CST")
	geo_2  = GeoRegion("OTREC_SAN")
end

# ╔═╡ 31016a69-296c-47e7-a440-bde58d68d759
begin
	geovec = Vector{GeoRegion}(undef,2)
	geovec[1] = geo_1
	geovec[2] = geo_2
end

# ╔═╡ 14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
begin
	rainb = 5  : 5 : 95; nr = length(rainb)
	rainu = 10 : 5 : 100
	rainm = 5  : 5 : 100
	pwgtb = 100 : 25 : 925; np = length(pwgtb)
	pwgtu = 125 : 25 : 950
	pwgtm = 100 : 25 : 950

	binmat = zeros(nr,np,12)
	binstd = zeros(nr,np,12)
	binnum = zeros(nr,np,12)
	
	for igeo = 1 : length(geovec)
		rain = []
		rHDO = []
		pwgt = []
		for imo = 8 : 12
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-RAINNC"))
			lon  = ds["longitude"][:]
			lat  = ds["latitude"][:]
			irin = ds["RAINNC"][(lon.<(geovec[igeo].E)).&(lon.>(geovec[igeo].W)).&(lat .>geovec[igeo].S).&(lat.<geovec[igeo].N)]
			close(ds)
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-O18_RAINNC"))
			iHDO = ds["O18_RAINNC"][(lon.<(geovec[igeo].E)).&(lon.>(geovec[igeo].W)).&(lat .>geovec[igeo].S).&(lat.<geovec[igeo].N)]
			iHDO = (iHDO./irin .-1) *1000
			close(ds)
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-p_wwgt"))
			ipwg = ds["p_wwgt"][(lon.<(geovec[igeo].E)).&(lon.>(geovec[igeo].W)).&(lat .>geovec[igeo].S).&(lat.<geovec[igeo].N)]
			close(ds)
			rain = vcat(rain,irin[:])
			rHDO = vcat(rHDO,iHDO[:])
			pwgt = vcat(pwgt,ipwg[:])
		end
		rainl = rain[:]
		pwgtl = pwgt[:] ./ 100
		rHDOl = rHDO[:]
	
		for ip = 1 : np, ir = 1 : nr
			try
				binmat[ir,ip,igeo] = mean(rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir])])
				binnum[ir,ip,igeo] = sum((pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir]))
			catch
				binmat[ir,ip,igeo] = NaN
				binnum[ir,ip,igeo] = NaN
			end
	
		end
		for ip = 1 : np, ir = 1 : nr
			try
				binstd[ir,ip,igeo] = sum((rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir])] .- binmat[ir,ip,igeo]).^2)
			catch
				binstd[ir,ip,igeo] = NaN
			end
		end
	end
end

# ╔═╡ e884be27-b9ec-45ba-b800-afccbfb54a78
begin
	pplt.close()
	f1,a1 = pplt.subplots(ncols=4,axwidth=1.2,aspect=1/2,wspace=[0,1.5,0])
	
	c1_1 = a1[1].pcolormesh(
		rainm,pwgtm,binmat[:,:,1]'*NaN,
		cmap="viridis",levels=-15:0.5:-7.5,extend="both"
	)
	c1_2 = a1[1].pcolormesh(
		rainm,pwgtm,sqrt.(binstd[:,:,1])'.*NaN,
		cmap="fire",levels=0:0.2:2,extend="both"
	)
	
	for ii = 1 : 2 : 3
		jj = Int((ii+1)/2)
		binni = binnum[:,:,jj]
		binii = binmat[:,:,jj]
		binii[binni.<=5] .= NaN
		a1[ii].pcolormesh(
			rainm,pwgtm,binii',
			cmap="viridis",levels=-15:0.5:-7.5,extend="both"
		)
		a1[ii].contour(
			rainm,pwgtm,binni'./sum(binni[.!isnan.(binni)])*np*nr,linestyle="--",
			cmap="fire",levels=[1,2,5,10,20],extend="both",
		)
	end

	for ii = 2 : 2 : 4
		jj = Int((ii)/2)
		binii = sqrt.(binstd[:,:,jj]./binnum[:,:,jj])
		binni = binnum[:,:,jj]
		binii[binni.<=5] .= NaN
		a1[ii].pcolormesh(
			rainm,pwgtm,binii',levels=0:0.2:2,extend="both"
		)
	end
	a1[2].format(rtitle="(a) CST")
	a1[4].format(rtitle="(b) SAN")

	for ax = a1
		ax.format(
			ylim=(maximum(pwgtm),minimum(pwgtm)),ylabel=L"p$_w$ / hPa",
			xlim=(minimum(rainm),maximum(rainm)),xlabel=L"Rain Rate / mm day$^{-1}$",
			xlocator=0:20:90
		)
	end

	lbl = L"$\delta^{18}$O / $\perthousand$"
	f1.colorbar(c1_1,loc="b",cols=(1,2),length=0.75,label= L"$\mu$ " * lbl)
	f1.colorbar(c1_2,loc="b",cols=(3,4),length=0.75,label=L"$\sigma$ " * lbl)
	f1.savefig(plotsdir("04e-IsotopeWRFsample.png"),transparent=false,dpi=200)
	load(plotsdir("04e-IsotopeWRFsample.png"))
end

# ╔═╡ dd6c0275-f939-4bdd-9cce-df59acebd1af
md"
### B. Bin separately according to Longitude (which is an analogue to land/sea)
"

# ╔═╡ fa2f0411-6a08-4f9a-bdab-0303d5d26291
begin
	geolist = Vector{GeoRegion}(undef,10)
	for ii = 1 : 10
		geolist[ii] = RectRegion(
			"TMP_$ii","OTREC","OTREC Colombia Ocean Coast",
			[7,3,-79.5+(ii-1)*0.5,-80+(ii-1)*0.5],savegeo=false
		)
	end
end

# ╔═╡ a8c95b81-2ac0-4c05-86a5-6d11616bf458
begin
	bin2 = zeros(nr,np,10)
	binn = zeros(nr,np,10)
		
	for igeo = 1 : length(geolist)
		rain = []
		rHDO = []
		pwgt = []
		for imo = 8 : 12
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-RAINNC"))
			lon  = ds["longitude"][:]
			lat  = ds["latitude"][:]
			irin = ds["RAINNC"][(lon.<(geolist[igeo].E)).&(lon.>(geolist[igeo].W)).&(lat .>geolist[igeo].S).&(lat.<geolist[igeo].N)]
			close(ds)
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-O18_RAINNC"))
			iHDO = ds["O18_RAINNC"][(lon.<(geolist[igeo].E)).&(lon.>(geolist[igeo].W)).&(lat .>geolist[igeo].S).&(lat.<geolist[igeo].N)]
			iHDO = (iHDO./irin .-1) *1000
			close(ds)
			ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-p_wwgt"))
			ipwg = ds["p_wwgt"][(lon.<(geolist[igeo].E)).&(lon.>(geolist[igeo].W)).&(lat .>geolist[igeo].S).&(lat.<geolist[igeo].N)]
			close(ds)
			rain = vcat(rain,irin[:])
			rHDO = vcat(rHDO,iHDO[:])
			pwgt = vcat(pwgt,ipwg[:])
		end
		rainl = rain[:]
		pwgtl = pwgt[:] ./ 100
		rHDOl = rHDO[:]
	
		for ip = 1 : np, ir = 1 : nr
			try
				bin2[ir,ip,igeo] = mean(rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir])])
				binn[ir,ip,igeo] = sum((pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir]))
			catch
				bin2[ir,ip,igeo] = NaN
				binn[ir,ip,igeo] = NaN
			end
	
		end
	end
end

# ╔═╡ 08391c66-72f6-44a2-9f8d-e1c12f954538
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=5,nrows=2,axwidth=1.2,aspect=2/3,wspace=1.5)

	c2 = a2[1].pcolormesh(
		rainm,pwgtm,bin2[:,:,1]',
		cmap="viridis",levels=-15:0.5:-7.5,extend="both"
	)
	aa = binn[:,:,1]
	a2[1].contour(
		rainm,pwgtm,aa'./sum(aa[.!isnan.(aa)])*np*nr,linestyle="--",
		cmap="fire",levels=[1,2,5,10,20],extend="both",
	)
	a2[1].text(50,825,"E = $(-79.5)")
	a2[1].text(50,900,"W = $(-80)")
	for ii = 1 : 10
		a2[ii].pcolormesh(
			rainm,pwgtm,bin2[:,:,ii]',
			cmap="viridis",levels=-15:0.5:-7.5,extend="both"
		)
		aii = binn[:,:,ii]
		a2[ii].contour(
			rainm,pwgtm,aii'./sum(aii[.!isnan.(aii)])*np*nr,linestyle="--",
			cmap="fire",levels=[1,2,5,10,20],extend="both",
		)
        a2[ii].text(50,825,"E = $(-79.5+(ii-1)*0.5)")
        a2[ii].text(50,900,"W = $(-80+(ii-1)*0.5)")
	end

	for ax = a2
		ax.format(
			ylim=(maximum(pwgtm),minimum(pwgtm)),ylabel=L"p$_w$ / hPa",
			xlim=(minimum(rainm),maximum(rainm)),xlabel=L"Rain Rate / mm day$^{-1}$",
			xlocator=0:20:100
		)
	end

	f2.colorbar(c2,length=0.6,label=L"$\delta^{18}$O / $\perthousand$")
	f2.savefig(plotsdir("04e-IsotopeWRFsample2.png"),transparent=false,dpi=400)
	load(plotsdir("04e-IsotopeWRFsample2.png"))
end

# ╔═╡ 8a5bb14a-035b-470a-9efe-8ec483711a79
md"
### C. Bin accordingly to Land-Sea Mask
"

# ╔═╡ Cell order:
# ╟─c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
# ╟─3c4b037f-c854-4148-ae6e-fff50bbb51cf
# ╟─4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
# ╟─c853d30e-5a51-47f4-a6be-88ab857253c1
# ╟─3922aee7-0e55-4074-af98-9f09919df7df
# ╟─31016a69-296c-47e7-a440-bde58d68d759
# ╟─14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
# ╟─e884be27-b9ec-45ba-b800-afccbfb54a78
# ╟─dd6c0275-f939-4bdd-9cce-df59acebd1af
# ╟─fa2f0411-6a08-4f9a-bdab-0303d5d26291
# ╟─a8c95b81-2ac0-4c05-86a5-6d11616bf458
# ╟─08391c66-72f6-44a2-9f8d-e1c12f954538
# ╟─8a5bb14a-035b-470a-9efe-8ec483711a79
