### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 906ab755-843f-4647-9ddc-1366d3e6199c
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 4208d36e-eec6-44b5-9316-98045ee88e51
begin
	@quickactivate "ColombiaIsotope"
	using DataFrames
	using Dates
	using DelimitedFiles
	using GLM
	using NCDatasets
	using Printf
	using Statistics
	using XLSX
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ db91546a-6db9-11ec-3d37-0f5cdb63ec1a
md"
# 02d. Rainfall vs Isotopes

In this notebook, we then proceed to compare station rainfall directly against the corresponding isotopic data.  We do this first for monthly data, and then for station data.
"

# ╔═╡ f737f273-65a6-4c55-983c-8d486e63f96d
begin
	fields = [
		:"Station",:"Collection Date/Time (Start)",
		:"Collection Date/Time (end)",
		:"δ2H, in ‰",:"δ18O, in ‰",:"Totalizer precipitation (mm)",
	]
	mdf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Monthly")...)
	mdf = mdf[:,fields]; mgdf = groupby(mdf,"Station"); nmo = length(mgdf)
	ddf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Daily")...)
	ddf = ddf[:,fields]; dgdf = groupby(ddf,"Station"); ndy = length(dgdf)
	md"Loading Isotope and Precipitation Data information ..."
end

# ╔═╡ 0ffac087-dbf1-49cc-9c16-f056966cfb9d
function precipitationvalue(prcp)
	for ip = 1 : length(prcp)
		if typeof(prcp[ip]) <: String
			prcp[ip] = parse(Float64,prcp[ip])
		elseif ismissing(prcp[ip])
			prcp[ip] = NaN
		end
	end

	return prcp
end

# ╔═╡ 3a7e5f50-a8c1-48e1-a456-45be985ebe72
md"
### A. Plotting Raw Daily and Monthly Values

In this section, we simply plot the raw data of station rainfall against the measured isotopic concentrations.
"

# ╔═╡ f5d67b19-6756-414e-a586-3c90c5c60b8f
begin
	pplt.close(); f1,a1 = pplt.subplots(nrows=3,ncols=5,axwidth=1.5)

	for istn = 1 : 15

		prcp_mo = precipitationvalue(mgdf[istn][:,:"Totalizer precipitation (mm)"])
		δ18O_mo = mgdf[istn][:,:"δ18O, in ‰"]
		for idata = 1 : length(δ18O_mo)
			if ismissing(δ18O_mo[idata]); δ18O_mo[idata] = NaN end
		end
		stn = mgdf[istn][1,:"Station"]
		a1[istn].scatter(prcp_mo/30,δ18O_mo,s=5)
		a1[istn].format(urtitle="$(stn)")
	
	end

	for ax in a1
		ax.format(
			xlim=(-1,31),xlocator=0:5:40,xlabel=L"Precipitation Rate / mm day$^{-1}$",
			ylim=(-21,6),ylocator=-25:5:10,ylabel=L"$\delta^{18}$O"
		)
	end
	
	f1.savefig(plotsdir("02d-prcpvO-mo.png"),transparent=false,dpi=150)
	load(plotsdir("02d-prcpvO-mo.png"))
end

# ╔═╡ dae18a65-da28-4328-a5e4-5d78926d7831
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=3,ncols=5,axwidth=1.5)

	for istn = 1 : 15

		prcp_mo = precipitationvalue(mgdf[istn][:,:"Totalizer precipitation (mm)"])
		δ2H_mo  = mgdf[istn][:,:"δ2H, in ‰"]
		for idata = 1 : length(δ2H_mo)
			if ismissing(δ2H_mo[idata]); δ2H_mo[idata] = NaN end
		end
		stn = mgdf[istn][1,:"Station"]
		a2[istn].scatter(prcp_mo/30,δ2H_mo,s=5)
		a2[istn].format(urtitle="$(stn)")
	
	end

	for ax in a2
		ax.format(
			xlim=(-1,31),xlocator=0:5:40,xlabel=L"Precipitation Rate / mm day$^{-1}$",
			ylim=(-155,35),ylocator=-150:25:25,ylabel=L"$\delta^{2}$H"
		)
	end
	
	f2.savefig(plotsdir("02d-prcpvH-mo.png.png"),transparent=false,dpi=150)
	load(plotsdir("02d-prcpvH-mo.png.png"))
end

# ╔═╡ 68c9ba0a-c329-43fb-9758-97a00bd78cac
begin
	pplt.close(); f3,a3 = pplt.subplots(ncols=4,axwidth=1.5)

	for istn = 1 : 4

		prcp = precipitationvalue(dgdf[istn][:,:"Totalizer precipitation (mm)"])
		dbeg = dgdf[istn][:,:"Collection Date/Time (Start)"]
		dend = dgdf[istn][:,:"Collection Date/Time (end)"]
		day  = Dates.value.(dend-dbeg)
		δ18O = dgdf[istn][:,:"δ18O, in ‰"]
		for idata = 1 : length(δ18O)
			if ismissing(δ18O[idata]); δ18O[idata] = NaN end
		end
		stn = dgdf[istn][1,:"Station"]
		a3[istn].scatter(prcp./day,δ18O,s=5)
		a3[istn].format(urtitle="$(stn)")
	
	end

	for ax in a3
		ax.format(
			xlim=(-5,125),xlocator=0:30:150,
			xlabel=L"Precipitation Rate / mm day$^{-1}$",
			ylim=(-17,2),ylocator=-25:5:10,ylabel=L"$\delta^{18}$O"
		)
	end
	
	f3.savefig(plotsdir("02d-prcpvO-dy.png"),transparent=false,dpi=150)
	load(plotsdir("02d-prcpvO-dy.png"))
end

# ╔═╡ 5cb82472-8f7b-40ed-9ee1-55bb01775f7f
begin
	pplt.close(); f4,a4 = pplt.subplots(ncols=4,axwidth=1.5)

	for istn = 1 : 4

		prcp = precipitationvalue(dgdf[istn][:,:"Totalizer precipitation (mm)"])
		dbeg = dgdf[istn][:,:"Collection Date/Time (Start)"]
		dend = dgdf[istn][:,:"Collection Date/Time (end)"]
		day  = Dates.value.(dend-dbeg)
		δ2H  = dgdf[istn][:,:"δ2H, in ‰"]
		for idata = 1 : length(δ2H)
			if ismissing(δ2H[idata]); δ2H[idata] = NaN end
		end
		stn = dgdf[istn][1,:"Station"]
		a4[istn].scatter(prcp./day,δ2H,s=5)
		a4[istn].format(urtitle="$(stn)")
	
	end

	for ax in a4
		ax.format(
			xlim=(-5,125),xlocator=0:30:150,
			xlabel=L"Precipitation Rate / mm day$^{-1}$",
			ylim=(-130,30),ylocator=-150:25:25,ylabel=L"$\delta^{2}$H"
		)
	end
	
	f4.savefig(plotsdir("02d-prcpvH-dy.png"),transparent=false,dpi=150)
	load(plotsdir("02d-prcpvH-dy.png"))
end

# ╔═╡ 19e224b1-4078-461b-964b-d776b0590ee0
md"
### B. Converting from Daily to Monthly data

The stations where daily data are available are (with the exception of San Andres) located along the Pacific coastline.  However, stations where monthly data is directly available are not located in the same regions.  Therefore, we perform monthly averages for the stations with daily data in order to see if there are any notable differences between the precipitation-isotope relationship along the Pacific coast compared to inland Colombia and the Caribbean coastline.
"

# ╔═╡ ba890902-e8ad-44c9-879d-c85999cbdcb9
begin
	pds = NCDataset(datadir("processed.nc"))
	dtvec = pds["time"][:]
	prcps = pds["prcps"][:]
	δ18Oμ = pds["δ18Oμ"][:]
	δ2Hμ  = pds["δ2Hμ"][:]
	close(pds)
	md"Loading data that was processed in notebook `02b` ..."
end

# ╔═╡ 2e0a9bc5-7036-4727-abc3-cf65474672f5
function averagenday(data::Array;n=7)

	ndt = size(data,1); nstn = size(data,2)
	ndn = Int(floor(ndt/n))
	ndata = zeros(ndn,nstn)
	
	for istn = 1 : nstn, idt = 1 : ndn
		idata = @view data[(idt-1)*n .+ (1:n),istn]
		if sum(isnan.(idata)) < (0.6*n)
			ndata[idt,istn] = mean(idata[.!isnan.(idata)])
		else
			ndata[idt,istn] = NaN
		end
	end
	
	return ndata
        
end

# ╔═╡ d65e32c7-0b45-4c44-9de9-e26eadc95f5a
function isotopenday(prcp::Array,iso::Array;n=7)

	ndt = size(prcp,1); nstn = size(prcp,2)
	ndn = Int(floor(ndt/n))
	niso = zeros(ndn,nstn)
	
	for istn = 1 : nstn, idt = 1 : ndn
		iprcp = @view prcp[(idt-1)*n .+ (1:n),istn]
		iiso  = @view iso[(idt-1)*n .+ (1:n),istn]
		iiso  = iprcp .* iiso
		inan  = .!isnan.(iiso)
		if sum(inan) > (0.6*n)
			niso[idt,istn] = sum(iiso[inan]) / sum(iprcp[inan])
		else
			niso[idt,istn] = NaN
		end
	end
	
	return niso
        
end

# ╔═╡ deb8d418-a6ea-47e9-8e04-1c2698f1c4a8
begin
	prcps_7d = averagenday(prcps)
	prcps_mo = averagenday(prcps,n=30)
	md"Performing weekly and monthly averages on precipitation data ..."
end

# ╔═╡ 01f7847d-0d08-4cd9-a4af-303aacc0382c
begin
	δ18Oμ_7d = isotopenday(prcps,δ18Oμ)
	δ2Hμ_7d  = isotopenday(prcps,δ2Hμ)
	δ18Oμ_mo = isotopenday(prcps,δ18Oμ,n=30)
	δ2Hμ_mo  = isotopenday(prcps,δ2Hμ,n=30)
	md"Performing weekly and monthly averages on isotope data ..."
end

# ╔═╡ 4730804e-b183-461f-80fb-a95c62680d0f
begin
	pplt.close(); f5,a5 = pplt.subplots(ncols=3,nrows=2,axwidth=1.5)
	
	a5[1].scatter(prcps[:,1],δ18Oμ[:,1],s=7)
	a5[1].scatter(prcps[:,2],δ18Oμ[:,2],s=7)
	a5[1].scatter(prcps[:,3],δ18Oμ[:,3],s=7)
	a5[1].scatter(prcps[:,4],δ18Oμ[:,4],s=7)
	
	a5[2].scatter(prcps_7d[:,1],δ18Oμ_7d[:,1],s=7)
	a5[2].scatter(prcps_7d[:,2],δ18Oμ_7d[:,2],s=7)
	a5[2].scatter(prcps_7d[:,3],δ18Oμ_7d[:,3],s=7)
	a5[2].scatter(prcps_7d[:,4],δ18Oμ_7d[:,4],s=7)
	
	a5[3].scatter(prcps_mo[:,1],δ18Oμ_mo[:,1])
	a5[3].scatter(prcps_mo[:,2],δ18Oμ_mo[:,2])
	a5[3].scatter(prcps_mo[:,3],δ18Oμ_mo[:,3])
	a5[3].scatter(prcps_mo[:,4],δ18Oμ_mo[:,4])
	
	a5[4].scatter(prcps[:,1],δ2Hμ[:,1],s=7)
	a5[4].scatter(prcps[:,2],δ2Hμ[:,2],s=7)
	a5[4].scatter(prcps[:,3],δ2Hμ[:,3],s=7)
	a5[4].scatter(prcps[:,4],δ2Hμ[:,4],s=7)
	
	a5[5].scatter(prcps_7d[:,1],δ2Hμ_7d[:,1],s=7,label="San Andres",legend="b")
	a5[5].scatter(prcps_7d[:,2],δ2Hμ_7d[:,2],s=7,label="Quibdó",legend="b")
	a5[5].scatter(prcps_7d[:,3],δ2Hμ_7d[:,3],s=7,label="Bahía Solano",legend="b")
	a5[5].scatter(
		prcps_7d[:,4],δ2Hμ_7d[:,4],s=7,
		label="Buenoventura",legend="b",legend_kw=Dict(
		"ncol"=>4,"frame"=>false
	))
	
	a5[6].scatter(prcps_mo[:,1],δ2Hμ_mo[:,1])
	a5[6].scatter(prcps_mo[:,2],δ2Hμ_mo[:,2])
	a5[6].scatter(prcps_mo[:,3],δ2Hμ_mo[:,3])
	a5[6].scatter(prcps_mo[:,4],δ2Hμ_mo[:,4])

	a5[1].format(
		leftlabels=[
			L"$\delta^{18}O$ / $\perthousand$",
			L"$\delta^{2}H$ / $\perthousand$"
		],xlabel=L"Precipitation Rate / mm day$^{-1}$"
	)
	
	f5.savefig(plotsdir("02d-dy2mo.png"),transparent=false,dpi=150)
	load(plotsdir("02d-dy2mo.png"))
end

# ╔═╡ 35b86280-642d-45ef-8270-d226dd9ad26a
md"
### C. Comparing Daily vs. Monthly San Andres Data
"

# ╔═╡ 77944ed0-5af0-4c80-8cf3-9e5ee76cc9d0
begin
	pplt.close(); f6,a6 = pplt.subplots(ncols=2,axwidth=1.5,sharey=0)
	
	prcp_mo = precipitationvalue(mgdf[12][:,:"Totalizer precipitation (mm)"])
	δ18O_mo = mgdf[12][:,:"δ18O, in ‰"]
	δ2H_mo = mgdf[12][:,:"δ2H, in ‰"]
	for idata = 1 : length(δ18O_mo)
		if ismissing(δ18O_mo[idata]); δ18O_mo[idata] = NaN end
		if ismissing(δ2H_mo[idata]); δ2H_mo[idata] = NaN end
	end
	
	a6[1].scatter(prcps_mo[:,1],δ18Oμ_mo[:,1],s=15)
	a6[1].scatter(prcp_mo/30,δ18O_mo,s=15)

	a6[2].scatter(
		prcps_mo[:,1],δ2Hμ_mo[:,1],s=15,
		label="30-Day Average of Daily Data",legend="r"
	)
	a6[2].scatter(prcp_mo/30,δ2H_mo,s=15
		,label="Monthly Data",legend="r",legend_kw=Dict(
		"ncol"=>1,"frame"=>false
	))

	a6[1].format(
		ylabel=L"$\delta^{18}O$ / $\perthousand$",
		xlabel=L"Precipitation Rate / mm day$^{-1}$",
		suptitle="San Andres Monthly Data"
	)
	a6[2].format(ylabel=L"$\delta^{2}H$ / $\perthousand$")
	
	f6.savefig(plotsdir("02d-SanAndres.png"),transparent=false,dpi=150)
	load(plotsdir("02d-SanAndres.png"))
end

# ╔═╡ Cell order:
# ╟─db91546a-6db9-11ec-3d37-0f5cdb63ec1a
# ╟─906ab755-843f-4647-9ddc-1366d3e6199c
# ╟─4208d36e-eec6-44b5-9316-98045ee88e51
# ╟─f737f273-65a6-4c55-983c-8d486e63f96d
# ╟─0ffac087-dbf1-49cc-9c16-f056966cfb9d
# ╟─3a7e5f50-a8c1-48e1-a456-45be985ebe72
# ╠═f5d67b19-6756-414e-a586-3c90c5c60b8f
# ╟─dae18a65-da28-4328-a5e4-5d78926d7831
# ╟─68c9ba0a-c329-43fb-9758-97a00bd78cac
# ╟─5cb82472-8f7b-40ed-9ee1-55bb01775f7f
# ╟─19e224b1-4078-461b-964b-d776b0590ee0
# ╟─ba890902-e8ad-44c9-879d-c85999cbdcb9
# ╟─2e0a9bc5-7036-4727-abc3-cf65474672f5
# ╟─d65e32c7-0b45-4c44-9de9-e26eadc95f5a
# ╟─deb8d418-a6ea-47e9-8e04-1c2698f1c4a8
# ╟─01f7847d-0d08-4cd9-a4af-303aacc0382c
# ╟─4730804e-b183-461f-80fb-a95c62680d0f
# ╟─35b86280-642d-45ef-8270-d226dd9ad26a
# ╟─77944ed0-5af0-4c80-8cf3-9e5ee76cc9d0
