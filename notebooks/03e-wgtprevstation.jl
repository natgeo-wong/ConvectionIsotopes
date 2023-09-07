### A Pluto.jl notebook ###
# v0.19.26

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
	@quickactivate "ColombiaIsotope"
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
	
	md"Loading modules for the ColumbiaIsotope project..."
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

# ╔═╡ 4bcf89d6-d9a2-4777-9e2f-1cf984b88712
begin
	fields = [
		:"Station",:"DateMid",
		:"δ2H, in ‰",:"δ18O, in ‰",:"Totalizer precipitation (mm)",
	]
	mdf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Monthly"))
	mdf = mdf[:,fields]; mgdf = groupby(mdf,"Station"); nmo = length(mgdf)
	md"Loading Isotope and Precipitation Data information ..."
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

# ╔═╡ 667fe7de-db5d-4649-b14a-0493f81c29ea
begin
	pplt.close(); fig,axs = pplt.subplots(axwidth=2,sharey=0,sharex=0)
	
	for istn = [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16]
	
		stn = mgdf[istn][1,:"Station"]
		lon = infomo[istn,2] + 360; ilon = argmin(abs.(lon .- lds.lon))
		lat = infomo[istn,3]; ilat = argmin(abs.(lat .- lds.lat))
		ndt = size(mgdf[istn],1)
		
		prc = mgdf[istn][:,:"Totalizer precipitation (mm)"]
		δOμ = mgdf[istn][:,:"δ18O, in ‰"]
		p_wwgt = zeros(Float32,ndt)
	
		for idt = 1 : ndt
	
			dt = mgdf[istn][idt,:"DateMid"]
			yr = year(dt)
			mo = month(dt)
	
			ids = read(e5ds,evar,ereg,Date(yr))
			p_wwgtii = ids["p_wwgt"][ilon,ilat,mo] / 100
			if ismissing(p_wwgtii); p_wwgtii = NaN end
			p_wwgt[idt] = p_wwgtii
			close(ids)

			if ismissing(prc[idt])
				prc[idt] = NaN
			end

			if ismissing(δOμ[idt])
				δOμ[idt] = NaN
			end
	
		end
	
		axs[1].scatter(prc/30,p_wwgt,c=δOμ,levels=-25:5)
	
	end

	axs[1].format(
		ylabel=L"$p_{wwtg}$ / hPa",xlim=(0,25),ylim=(800,200),
		xlabel=L"$\delta^{18}$O / $\perthousand$"
	)
	
	fig.savefig(plotsdir("03e-moprcpvpwwgt-monthlystn.png"),transparent=false,dpi=150)
	load(plotsdir("03e-moprcpvpwwgt-monthlystn.png"))
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
	prcps = ds["prcps"][:]
	δ18Oμ = ds["δ2Hμ"][:]
	time  = ds["time"][:]
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
	md"Preallocation of arrays ..."
end

# ╔═╡ cc4af5ad-ad08-425a-92ea-2284f9d03b1a
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=8,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0,2,0])
	
	yrvec = year.(time)
	movec = month.(time)
	dtvec = Date(yrvec[1],movec[1]) : Month(1) : Date(yrvec[end],movec[end])
	yrmov = Date.(yrvec,movec)

	abin .= 0; bbin .= 0; cbin .= 0; dbin .= 0
	anum .= 0; bnum .= 0; cnum .= 0; dnum .= 0

	lvls = -75:5:0
	for dtii in dtvec
		ind = (yrmov .== dtii)
		ids = read(e5ds,evar,ereg,Date(dtii))
		for istn = 1 : 1
			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpii = prcps[ind,istn]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			δ18Oii = δ18Oμ[ind,istn]
			δ18Oii = mean(δ18Oii[.!isnan.(δ18Oii)])
			pwwgti = ids["p_wwgt"][ilon,ilat,month(dtii)] / 100

			if ismissing(pwwgti); pwwgti = NaN end

			if !isnan(prcpii) && !isnan(pwwgti)
				rind = argmin(abs.(prcpii.-rpnt))
				pind = argmin(abs.(pwwgti.-ppnt))
				abin[rind,pind] += δ18Oii; anum[rind,pind] += 1
				bbin[rind,pind] += δ18Oii; bnum[rind,pind] += 1
			end
			
			# a2[1].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)			
			# a2[2].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)
		end
		for istn = [3,4]
			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpii = prcps[ind,istn]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			δ18Oii = δ18Oμ[ind,istn]
			δ18Oii = mean(δ18Oii[.!isnan.(δ18Oii)])
			pwwgti = ids["p_wwgt"][ilon,ilat,month(dtii)] / 100

			if ismissing(pwwgti); pwwgti = NaN end

			if !isnan(prcpii) && !isnan(pwwgti)
				rind = argmin(abs.(prcpii.-rpnt))
				pind = argmin(abs.(pwwgti.-ppnt))
				abin[rind,pind] += δ18Oii; anum[rind,pind] += 1
				cbin[rind,pind] += δ18Oii; cnum[rind,pind] += 1
			end
			
			# a2[1].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)
			# a2[3].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)
		end
		for istn = 2 : 2
			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpii = prcps[ind,istn]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			δ18Oii = δ18Oμ[ind,istn]
			δ18Oii = mean(δ18Oii[.!isnan.(δ18Oii)])
			pwwgti = ids["p_wwgt"][ilon,ilat,month(dtii)] / 100

			if ismissing(pwwgti); pwwgti = NaN end

			if !isnan(prcpii) && !isnan(pwwgti)
				rind = argmin(abs.(prcpii.-rpnt))
				pind = argmin(abs.(pwwgti.-ppnt))
				abin[rind,pind] += δ18Oii; anum[rind,pind] += 1
				dbin[rind,pind] += δ18Oii; dnum[rind,pind] += 1
			end
			
			# a2[1].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)
			# a2[4].scatter(prcpii,pwwgti,c=δ18Oii,cmap="viridis",levels=lvls,extend="both",s=10)
		end
		close(ids)
	end

	a2[1].pcolormesh(rbin,pbin,(abin./anum)',cmap="viridis",levels=lvls,extend="both")
	a2[3].pcolormesh(rbin,pbin,(bbin./bnum)',cmap="viridis",levels=lvls,extend="both")
	a2[5].pcolormesh(rbin,pbin,(cbin./cnum)',cmap="viridis",levels=lvls,extend="both")
	a2[7].pcolormesh(rbin,pbin,(dbin./dnum)',cmap="viridis",levels=lvls,extend="both")

	anum[iszero.(anum)] .= NaN
	bnum[iszero.(bnum)] .= NaN
	cnum[iszero.(cnum)] .= NaN
	dnum[iszero.(dnum)] .= NaN
	
	cbnum = a2[2].pcolormesh(rbin,pbin,anum',levels=0:10,extend="both")
	a2[4].pcolormesh(rbin,pbin,bnum',levels=0:10,extend="both")
	a2[6].pcolormesh(rbin,pbin,cnum',levels=0:10,extend="both")
	a2[8].pcolormesh(rbin,pbin,dnum',levels=0:10,extend="both")

	a2[2].format(urtitle="(a) All Stations",suptitle="30-Day Moving Average")
	a2[4].format(urtitle="(b) San Andres")
	a2[6].format(urtitle="(c) Buenaventura")
	a2[6].text(-23,295,"Bahia Solano",fontsize=10)
	a2[8].format(urtitle="(d) Quibdo")
	for ax in a2
		ax.format(
			ylabel=L"$p_\omega$ / hPa",
			xlabel=L"Rainfall Rate / mm day$^{-1}$",
			ylim=(950,100),xlim=(0,100),ylocator=200:200:800,xlocator=0:50:100
		)
	end

	f2.colorbar(cbnum,label="Number of Observations",length=0.75)
	f2.savefig(plotsdir("03e-moprcpvpwwgt-dailystn.png"),transparent=false,dpi=400)
	load(plotsdir("03e-moprcpvpwwgt-dailystn.png"))
end

# ╔═╡ 492b98ea-e5d8-40e3-ae89-df68d183c100
md"Now, we do a 7-day moving average ..."

# ╔═╡ bbcf49bb-4be0-463d-87ea-83f674a5729d
begin
	pplt.close()
	f3,a3 = pplt.subplots(ncols=8,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0,2,0])

	ii_7dy = 4 : 7 : length(time); n7dy = length(ii_7dy)
	time_7dy  = time[ii_7dy]
	prcps_7dy = reshape(prcps[1:n7dy*7,:],7,:,4)
	δ18Oμ_7dy = reshape(δ18Oμ[1:n7dy*7,:],7,:,4)

	abin .= 0; bbin .= 0; cbin .= 0; dbin .= 0
	anum .= 0; bnum .= 0; cnum .= 0; dnum .= 0

	for idy = 1 : n7dy

		dt  = time_7dy[idy]; dayii = day(dt)
		ids = read(e5dy,evar,ereg,dt,smooth=true,smoothtime=7)

		prcps_ii = prcps_7dy[:,idy,:]
		δ18Oμ_ii = δ18Oμ_7dy[:,idy,:]

		for istn in [1]

			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; prcpsii = mean(prcpsii[.!isnan.(prcpsii)])
			δ18Oμii = δ18Oμ_ii[:,istn]; δ18Oμii = mean(δ18Oμii[.!isnan.(δ18Oμii)])
			pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
			if ismissing(pwwgtii); pwwgtii = NaN end

			if !isnan(prcpsii) && !isnan(pwwgtii)
				rind = argmin(abs.(prcpsii.-rpnt))
				pind = argmin(abs.(pwwgtii.-ppnt))
				abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
				bbin[rind,pind] += δ18Oμii; bnum[rind,pind] += 1
			end

			# a3[1].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)
			# a3[2].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)

		end

		for istn in [3,4]

			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; prcpsii = mean(prcpsii[.!isnan.(prcpsii)])
			δ18Oμii = δ18Oμ_ii[:,istn]; δ18Oμii = mean(δ18Oμii[.!isnan.(δ18Oμii)])
			pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
			if ismissing(pwwgtii); pwwgtii = NaN end

			if !isnan(prcpsii) && !isnan(pwwgtii)
				rind = argmin(abs.(prcpsii.-rpnt))
				pind = argmin(abs.(pwwgtii.-ppnt))
				abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
				cbin[rind,pind] += δ18Oμii; cnum[rind,pind] += 1
			end

			# a3[1].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)
			# a3[3].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)

		end

		for istn in [2]

			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpsii = prcps_ii[:,istn]; prcpsii = mean(prcpsii[.!isnan.(prcpsii)])
			δ18Oμii = δ18Oμ_ii[:,istn]; δ18Oμii = mean(δ18Oμii[.!isnan.(δ18Oμii)])
			pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
			if ismissing(pwwgtii); pwwgtii = NaN end

			if !isnan(prcpsii) && !isnan(pwwgtii)
				rind = argmin(abs.(prcpsii.-rpnt))
				pind = argmin(abs.(pwwgtii.-ppnt))
				abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
				dbin[rind,pind] += δ18Oμii; dnum[rind,pind] += 1
			end

			# a3[1].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)
			# a3[4].scatter(prcpsii,pwwgtii,c=δ18Oμii,cmap="viridis",levels=lvls,extend="both",s=10)

		end

		close(ids)

	end

	c = a3[1].pcolormesh(rbin,pbin,(abin./anum)',cmap="viridis",levels=lvls,extend="both")
	a3[3].pcolormesh(rbin,pbin,(bbin./bnum)',cmap="viridis",levels=lvls,extend="both")
	a3[5].pcolormesh(rbin,pbin,(cbin./cnum)',cmap="viridis",levels=lvls,extend="both")
	a3[7].pcolormesh(rbin,pbin,(dbin./dnum)',cmap="viridis",levels=lvls,extend="both")

	anum[iszero.(anum)] .= NaN
	bnum[iszero.(bnum)] .= NaN
	cnum[iszero.(cnum)] .= NaN
	dnum[iszero.(dnum)] .= NaN
	
	a3[2].pcolormesh(rbin,pbin,anum',levels=0:10,extend="both")
	a3[4].pcolormesh(rbin,pbin,bnum',levels=0:10,extend="both")
	a3[6].pcolormesh(rbin,pbin,cnum',levels=0:10,extend="both")
	a3[8].pcolormesh(rbin,pbin,dnum',levels=0:10,extend="both")

	a3[2].format(urtitle="(a) All Stations",suptitle="7-Day Moving Average")
	a3[4].format(urtitle="(b) San Andres")
	a3[6].format(urtitle="(c) Buenaventura")
	a3[6].text(-23,295,"Bahia Solano",fontsize=10)
	a3[8].format(urtitle="(d) Quibdo")
	for ax in a3
		ax.format(
			ylabel=L"$p_\omega$ / hPa",
			xlabel=L"Rainfall Rate / mm day$^{-1}$",
			ylim=(950,100),xlim=(0,100),ylocator=200:200:800,xlocator=0:50:100
		)
	end
	f3.colorbar(c,label=L"$\delta^{2}$H / $\perthousand$",length=0.75,locator=-75:15:0)

	f3.savefig(plotsdir("03e-7dyprcpvpwwgt-dailystn.png"),transparent=false,dpi=400)
	load(plotsdir("03e-7dyprcpvpwwgt-dailystn.png"))
end

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─4bcf89d6-d9a2-4777-9e2f-1cf984b88712
# ╟─c8221168-ddb1-40f9-bc3c-c27393ce7568
# ╟─07c9b363-5455-4a6d-bd60-c53899b4d929
# ╟─25498e3c-a6f0-4686-8e5c-980894f70e38
# ╟─bd5c89ef-d3ba-4843-909a-1a7768316364
# ╟─ad40ce5d-cf0f-40d8-b2b4-05672d1d1f95
# ╟─849f8f15-a3c7-48a0-a1db-376cbb943a1b
# ╟─487c65b1-17f7-46fe-8c07-20c6e5f1c3f0
# ╟─73ab7928-b428-4792-af0e-1f44b42beee9
# ╟─667fe7de-db5d-4649-b14a-0493f81c29ea
# ╟─8927283a-c4e2-4c3a-b32f-850391eef9c7
# ╟─f4597370-4bc8-4763-8030-3a4dc89533b6
# ╠═2c4f6363-7b5a-40a0-aa47-251c53ae7367
# ╟─f0b67cf7-7159-44d2-bc9b-77306bddacf8
# ╟─b2dc480b-b6fc-471a-9724-1034415f4cfc
# ╟─cc4af5ad-ad08-425a-92ea-2284f9d03b1a
# ╟─492b98ea-e5d8-40e3-ae89-df68d183c100
# ╟─bbcf49bb-4be0-463d-87ea-83f674a5729d
