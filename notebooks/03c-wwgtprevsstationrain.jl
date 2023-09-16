### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ConvectionIsotopes"
	using DataFrames
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
	using Statistics
	using XLSX
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 04c. W-Weighted Column Mean Pressure vs Precipitation
"

# ╔═╡ 5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ c8221168-ddb1-40f9-bc3c-c27393ce7568
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 07c9b363-5455-4a6d-bd60-c53899b4d929
e5ds = ERA5Monthly(start=Date(2013),stop=Date(2021),path=datadir())

# ╔═╡ 25498e3c-a6f0-4686-8e5c-980894f70e38
e5dy = ERA5Daily(start=Date(2013),stop=Date(2021),path=datadir())

# ╔═╡ bd5c89ef-d3ba-4843-909a-1a7768316364
ereg = ERA5Region(GeoRegion("OTREC"))

# ╔═╡ ad40ce5d-cf0f-40d8-b2b4-05672d1d1f95
evar = SingleVariable("p_wwgt")

# ╔═╡ 849f8f15-a3c7-48a0-a1db-376cbb943a1b
lds = getLandSea(e5ds,ereg)

# ╔═╡ 487c65b1-17f7-46fe-8c07-20c6e5f1c3f0
md"
### B. Loading and Comparing Monthly Station Data
"

# ╔═╡ 73ab7928-b428-4792-af0e-1f44b42beee9
begin
	infomo  = stninfomo()
	md"Loading monthly station location information ..."
end

# ╔═╡ 8927283a-c4e2-4c3a-b32f-850391eef9c7
md"
### C. Loading and Comparing Daily Station Data
"

# ╔═╡ f4597370-4bc8-4763-8030-3a4dc89533b6
begin
	infody  = stninfody(); ndy = size(infody,1)
	md"Loading daily station location information ..."
end

# ╔═╡ 2c4f6363-7b5a-40a0-aa47-251c53ae7367
begin
	ds = NCDataset(datadir("processed.nc"))
	time  = ds["time"][:]; it = time .> Date(2019,1,31)
	time  = time[it]
	prcps = ds["prcps"][it,:]
	δ18Oμ = ds["δ2Hμ"][it,:]
	nstn = size(prcps,2)
	close(ds)
end

# ╔═╡ f0b67cf7-7159-44d2-bc9b-77306bddacf8
md"We do month-by-month analysis first ..."

# ╔═╡ b2dc480b-b6fc-471a-9724-1034415f4cfc
begin
	rbin = 0 : 10 : 100;    rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 100 : 50 : 950; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
	abin = zeros(length(rpnt),length(ppnt)); anum = zeros(length(rpnt),length(ppnt))
	bbin = zeros(length(rpnt),length(ppnt)); bnum = zeros(length(rpnt),length(ppnt))
	cbin = zeros(length(rpnt),length(ppnt)); cnum = zeros(length(rpnt),length(ppnt))
	dbin = zeros(length(rpnt),length(ppnt)); dnum = zeros(length(rpnt),length(ppnt))
	ebin = zeros(length(rpnt),length(ppnt)); enum = zeros(length(rpnt),length(ppnt))
	fbin = zeros(length(rpnt),length(ppnt)); fnum = zeros(length(rpnt),length(ppnt))
	aprc = zeros(length(rpnt),length(ppnt))
	bprc = zeros(length(rpnt),length(ppnt))
	cprc = zeros(length(rpnt),length(ppnt))
	dprc = zeros(length(rpnt),length(ppnt))
	eprc = zeros(length(rpnt),length(ppnt))
	fprc = zeros(length(rpnt),length(ppnt))
	md"Preallocation of arrays ..."
end

# ╔═╡ 9411b053-e4b8-420d-9e78-a91aab4fdcc1
ithr_30dy = 20

# ╔═╡ cc4af5ad-ad08-425a-92ea-2284f9d03b1a
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	abin .= 0; bbin .= 0; cbin .= 0; dbin .= 0; ebin .= 0; fbin .= 0
	anum .= 0; bnum .= 0; cnum .= 0; dnum .= 0; enum .= 0; fnum .= 0
	aprc .= 0; bprc .= 0; cprc .= 0; dprc .= 0; eprc .= 0; fprc .= 0

	lvls = -90 : 5 : 0
	for idy = 16 : (length(time)-15)

		dt  = time[idy]; dayii = day(dt)
		ids = read(e5dy,evar,ereg,dt,smooth=true,smoothtime=7,quiet=true)

		prcps_ii = prcps[idy.+(-15:15),:]
		δ18Oμ_ii = δ18Oμ[idy.+(-15:15),:]

		for istn in [1]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_30dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					bbin[rind,pind] += δ18Oμii; bnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; bprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [3,4]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_30dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					cbin[rind,pind] += δ18Oμii; cnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; cprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [2]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_30dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					dbin[rind,pind] += δ18Oμii; dnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; dprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 5 : 12

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_30dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					ebin[rind,pind] += δ18Oμii; enum[rind,pind] += 1
					aprc[rind,pind] += prcptot; eprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 1 : 4

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_30dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					fbin[rind,pind] += δ18Oμii; fnum[rind,pind] += 1
					fprc[rind,pind] += prcptot
				end
			end

		end

		close(ids)

	end

	c2_1 = a2[1].pcolormesh(rbin,pbin,(abin./aprc)',cmap="viridis",levels=lvls,extend="both")
	a2[3].pcolormesh(rbin,pbin,(fbin./fprc)',cmap="viridis",levels=lvls,extend="both")
	a2[5].pcolormesh(rbin,pbin,(ebin./eprc)',cmap="viridis",levels=lvls,extend="both")
	a2[7].pcolormesh(rbin,pbin,(bbin./bprc)',cmap="viridis",levels=lvls,extend="both")
	a2[9].pcolormesh(rbin,pbin,(cbin./cprc)',cmap="viridis",levels=lvls,extend="both")
	a2[11].pcolormesh(rbin,pbin,(dbin./dprc)',cmap="viridis",levels=lvls,extend="both")

	anum[iszero.(anum)] .= NaN
	bnum[iszero.(bnum)] .= NaN
	cnum[iszero.(cnum)] .= NaN
	dnum[iszero.(dnum)] .= NaN
	enum[iszero.(enum)] .= NaN
	fnum[iszero.(fnum)] .= NaN
	
	c2_2 = a2[2].pcolormesh(rbin,pbin,anum',levels=0:5:50,extend="both")
	a2[4].pcolormesh(rbin,pbin,fnum',levels=0:5:50,extend="both")
	a2[6].pcolormesh(rbin,pbin,enum',levels=0:5:50,extend="both")
	a2[8].pcolormesh(rbin,pbin,bnum',levels=0:5:50,extend="both")
	a2[10].pcolormesh(rbin,pbin,cnum',levels=0:5:50,extend="both")
	a2[12].pcolormesh(rbin,pbin,dnum',levels=0:5:50,extend="both")

	a2[2].format(urtitle="(a) All Stations",suptitle="30-Day Moving Average")
	a2[4].format(urtitle="(b) Colombia")
	a2[6].format(urtitle="(c) Costa Rica")
	a2[8].format(urtitle="(d) San Andres")
	a2[10].format(urtitle="(e) Buenaventura")
	a2[10].text(-23,300,"Bahia Solano",fontsize=10)
	a2[12].format(urtitle="(f) Quibdo")
	for ax in a2
		ax.format(
			ylabel=L"$p_\omega$ / hPa",
			xlabel=L"Rainfall Rate / mm day$^{-1}$",
			ylim=(950,100),xlim=(0,100),ylocator=200:200:800,xlocator=0:50:100
		)
	end
	f2.colorbar(c2_1,label=L"$\delta^{2}$H / $\perthousand$",locator=-75:15:0,minorlocator=-90:5:0,row=[1])
	f2.colorbar(c2_2,label="Number of Observations",row=[2])

	f2.savefig(plotsdir("03c-30dyprcpvpwwgt-dailystn.png"),transparent=false,dpi=400)
	load(plotsdir("03c-30dyprcpvpwwgt-dailystn.png"))
end

# ╔═╡ 492b98ea-e5d8-40e3-ae89-df68d183c100
md"Now, we do a 7-day moving average ..."

# ╔═╡ 189bb0b7-d3ae-4aa7-b9e3-1ee5dcc63e4d
ithr_7dy = 4

# ╔═╡ bbcf49bb-4be0-463d-87ea-83f674a5729d
begin
	pplt.close()
	f3,a3 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	abin .= 0; bbin .= 0; cbin .= 0; dbin .= 0; ebin .= 0; fbin .= 0
	anum .= 0; bnum .= 0; cnum .= 0; dnum .= 0; enum .= 0; fnum .= 0
	aprc .= 0; bprc .= 0; cprc .= 0; dprc .= 0; eprc .= 0; fprc .= 0

	for idy = 4 : (length(time)-3)

		dt  = time[idy]; dayii = day(dt)
		ids = read(e5dy,evar,ereg,dt,smooth=true,smoothtime=7,quiet=true)

		prcps_ii = prcps[idy.+(-3:3),:]
		δ18Oμ_ii = δ18Oμ[idy.+(-3:3),:]

		for istn in [1]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					bbin[rind,pind] += δ18Oμii; bnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; bprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [3,4]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					cbin[rind,pind] += δ18Oμii; cnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; cprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [2]

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					dbin[rind,pind] += δ18Oμii; dnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; dprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 5 : 12

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					ebin[rind,pind] += δ18Oμii; enum[rind,pind] += 1
					aprc[rind,pind] += prcptot; eprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 1 : 4

			lon = infody[istn,2]; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					fbin[rind,pind] += δ18Oμii; fnum[rind,pind] += 1
					fprc[rind,pind] += prcptot
				end
			end

		end

		close(ids)

	end

	c3_1 = a3[1].pcolormesh(rbin,pbin,(abin./aprc)',cmap="viridis",levels=lvls,extend="both")
	a3[3].pcolormesh(rbin,pbin,(fbin./fprc)',cmap="viridis",levels=lvls,extend="both")
	a3[5].pcolormesh(rbin,pbin,(ebin./eprc)',cmap="viridis",levels=lvls,extend="both")
	a3[7].pcolormesh(rbin,pbin,(bbin./bprc)',cmap="viridis",levels=lvls,extend="both")
	a3[9].pcolormesh(rbin,pbin,(cbin./cprc)',cmap="viridis",levels=lvls,extend="both")
	a3[11].pcolormesh(rbin,pbin,(dbin./dprc)',cmap="viridis",levels=lvls,extend="both")

	anum[iszero.(anum)] .= NaN
	bnum[iszero.(bnum)] .= NaN
	cnum[iszero.(cnum)] .= NaN
	dnum[iszero.(dnum)] .= NaN
	enum[iszero.(enum)] .= NaN
	fnum[iszero.(fnum)] .= NaN
	
	c3_2 = a3[2].pcolormesh(rbin,pbin,anum',levels=0:5:50,extend="both")
	a3[4].pcolormesh(rbin,pbin,fnum',levels=0:5:50,extend="both")
	a3[6].pcolormesh(rbin,pbin,enum',levels=0:5:50,extend="both")
	a3[8].pcolormesh(rbin,pbin,bnum',levels=0:5:50,extend="both")
	a3[10].pcolormesh(rbin,pbin,cnum',levels=0:5:50,extend="both")
	a3[12].pcolormesh(rbin,pbin,dnum',levels=0:5:50,extend="both")

	a3[2].format(urtitle="(a) All Stations",suptitle="7-Day Moving Average")
	a3[4].format(urtitle="(b) Colombia")
	a3[6].format(urtitle="(c) Costa Rica")
	a3[8].format(urtitle="(d) San Andres")
	a3[10].format(urtitle="(e) Buenaventura")
	a3[10].text(-23,300,"Bahia Solano",fontsize=10)
	a3[12].format(urtitle="(f) Quibdo")
	for ax in a3
		ax.format(
			ylabel=L"$p_\omega$ / hPa",
			xlabel=L"Rainfall Rate / mm day$^{-1}$",
			ylim=(950,100),xlim=(0,100),ylocator=200:200:800,xlocator=0:50:100
		)
	end
	f3.colorbar(c3_1,label=L"$\delta^{2}$H / $\perthousand$",locator=-75:15:0,minorlocator=-90:5:0,row=[1])
	f3.colorbar(c3_2,label="Number of Observations",row=[2])

	f3.savefig(plotsdir("03c-7dyprcpvpwwgt-dailystn.png"),transparent=false,dpi=400)
	load(plotsdir("03c-7dyprcpvpwwgt-dailystn.png"))
end

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─c8221168-ddb1-40f9-bc3c-c27393ce7568
# ╟─07c9b363-5455-4a6d-bd60-c53899b4d929
# ╟─25498e3c-a6f0-4686-8e5c-980894f70e38
# ╟─bd5c89ef-d3ba-4843-909a-1a7768316364
# ╟─ad40ce5d-cf0f-40d8-b2b4-05672d1d1f95
# ╟─849f8f15-a3c7-48a0-a1db-376cbb943a1b
# ╟─487c65b1-17f7-46fe-8c07-20c6e5f1c3f0
# ╟─73ab7928-b428-4792-af0e-1f44b42beee9
# ╟─8927283a-c4e2-4c3a-b32f-850391eef9c7
# ╟─f4597370-4bc8-4763-8030-3a4dc89533b6
# ╟─2c4f6363-7b5a-40a0-aa47-251c53ae7367
# ╟─f0b67cf7-7159-44d2-bc9b-77306bddacf8
# ╟─b2dc480b-b6fc-471a-9724-1034415f4cfc
# ╠═9411b053-e4b8-420d-9e78-a91aab4fdcc1
# ╟─cc4af5ad-ad08-425a-92ea-2284f9d03b1a
# ╟─492b98ea-e5d8-40e3-ae89-df68d183c100
# ╠═189bb0b7-d3ae-4aa7-b9e3-1ee5dcc63e4d
# ╟─bbcf49bb-4be0-463d-87ea-83f674a5729d
