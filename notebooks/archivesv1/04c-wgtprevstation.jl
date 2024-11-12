### A Pluto.jl notebook ###
# v0.18.1

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
	mdf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Monthly")...)
	mdf = mdf[:,fields]; mgdf = groupby(mdf,"Station"); nmo = length(mgdf)
	md"Loading Isotope and Precipitation Data information ..."
end

# ╔═╡ c8221168-ddb1-40f9-bc3c-c27393ce7568
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 07c9b363-5455-4a6d-bd60-c53899b4d929
e5ds = ERA5Monthly(dtbeg=Date(2013),dtend=Date(2021),eroot="/n/kuangdss01/lab")

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

# ╔═╡ 64c1fc09-2947-40fe-ac48-5a415916f5ac
infomo

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
	
	fig.savefig(plotsdir("04c-moprcpvpwwgt-monthlystn.png"),transparent=false,dpi=150)
	load(plotsdir("04c-moprcpvpwwgt-monthlystn.png"))
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

# ╔═╡ b7b82298-c7b7-4cae-9de4-400870f5edbd
infody

# ╔═╡ 2c4f6363-7b5a-40a0-aa47-251c53ae7367
begin
	ds = NCDataset(datadir("processed.nc"))
	prcps = ds["prcps"][:]
	δ18Oμ = ds["δ2Hμ"][:]
	time  = ds["time"][:]
	close(ds)
end

# ╔═╡ cc4af5ad-ad08-425a-92ea-2284f9d03b1a
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=1,sharey=0,axwidth=2)
	
	yrvec = year.(time)
	movec = month.(time)
	dtvec = Date(yrvec[1],movec[1]) : Month(1) : Date(yrvec[end],movec[end])
	yrmov = Date.(yrvec,movec)
	for dtii in dtvec
		ind = (yrmov .== dtii)
		ids = read(e5ds,evar,ereg,Date(dtii))
		for istn = 1 : ndy
			lon = infody[istn,2] + 360; ilon = argmin(abs.(lds.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lds.lat.-lat))
			prcpii = prcps[ind,istn]
			prcpii = mean(prcpii[.!isnan.(prcpii)])
			δ18Oii = δ18Oμ[ind,istn]
			δ18Oii = mean(δ18Oii[.!isnan.(δ18Oii)])
			pwwgti = ids["p_wwgt"][ilon,ilat,month(dtii)] / 100

			if ismissing(pwwgti); pwwgti = NaN end
			
			c = a2[1].scatter(prcpii,pwwgti,c=δ18Oii,levels=-100:20)
		end
	end
	c = a2[1].scatter(NaN,NaN,c=NaN,levels=-100:20)
	a2[1].format(
		ylabel="Weighted Column Pressure / hPa",
		xlabel=L"Rainfall Rate / mm day$^{-1}$",
		ylim=(800,300)
	)
	a2[1].colorbar(c,locator=-100:20:20,label=L"$\delta^2$H / $\perthousand$")

	f2.savefig(plotsdir("04c-moprcpvpwwgt-dailystn.png"),transparent=false,dpi=150)
	load(plotsdir("04c-moprcpvpwwgt-dailystn.png"))
end

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─4bcf89d6-d9a2-4777-9e2f-1cf984b88712
# ╟─c8221168-ddb1-40f9-bc3c-c27393ce7568
# ╟─07c9b363-5455-4a6d-bd60-c53899b4d929
# ╟─bd5c89ef-d3ba-4843-909a-1a7768316364
# ╟─ad40ce5d-cf0f-40d8-b2b4-05672d1d1f95
# ╠═849f8f15-a3c7-48a0-a1db-376cbb943a1b
# ╟─487c65b1-17f7-46fe-8c07-20c6e5f1c3f0
# ╟─73ab7928-b428-4792-af0e-1f44b42beee9
# ╠═64c1fc09-2947-40fe-ac48-5a415916f5ac
# ╟─667fe7de-db5d-4649-b14a-0493f81c29ea
# ╠═8927283a-c4e2-4c3a-b32f-850391eef9c7
# ╠═f4597370-4bc8-4763-8030-3a4dc89533b6
# ╠═b7b82298-c7b7-4cae-9de4-400870f5edbd
# ╠═2c4f6363-7b5a-40a0-aa47-251c53ae7367
# ╟─cc4af5ad-ad08-425a-92ea-2284f9d03b1a
