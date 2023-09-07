### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ b3bc56eb-43a7-4736-bd66-704529911d60
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 9071763e-f6ad-4468-ae48-369307a85263
begin
	@quickactivate "ConvectionIsotopes"
	using Dates
	using DataFrames
	using DelimitedFiles
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	using XLSX
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ a8431d46-46fe-11ec-2b8d-e39caffdabec
md"
# 03b - Corroborating Station Rainfall w/ GPM

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ 5b392d8a-d68b-44df-ab57-d0e81aedd668
md"
### A. Defining the NASAPrecipitation Dataset
"

# ╔═╡ 0c9b54e3-78a4-400a-9323-b716d31d327c
npd = IMERGFinalHH(start=Date(2019,7,1),stop=Date(2021,6,30),path=datadir())

# ╔═╡ 307f237d-9ba7-49a8-9107-fc000d1aa147
geo = GeoRegion("OTREC")

# ╔═╡ 6bfec0ba-ac32-4874-be90-f91815c9f188
begin
	fnc  = npd.ID * "-" * geo.ID * "-20190601.nc"
	ds   = NCDataset(joinpath(npd.datapath,geo.ID,"2019","06",fnc))
	lon  = ds["longitude"][:]
	lat  = ds["latitude"][:]
	close(ds)
end

# ╔═╡ 1d31fb5c-6af3-4dec-a2f3-fdc91a5e3856
md"
### B. Loading Station Information
"

# ╔═╡ 7cd13617-40e2-4bc5-ad0a-47b47581ec86
begin
	df = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Daily"))
	df = df[:,[
		:"Station",:"Longitude",:"Latitude",
		:"δ2H, in ‰",:"δ2H Stdev",:"δ18O, in ‰",:"δ18O Stdev",
		:"Collection Date/Time (Start)",:"Collection Date/Time (end)",
		:"Totalizer precipitation (mm)"
	]]
	gdf = groupby(df,"Station")
	md"Loading Isotope Data information ..."
end

# ╔═╡ 66b70aef-f079-4e37-9214-bb1f31ca227b
begin
	infody = stninfody()[1:4,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ bbfc6d8c-e49e-4444-beb2-7d10b907aa6b
begin
	ilon = zeros(Int,nstn)
	ilat = zeros(Int,nstn)
	for istn = 1 : nstn
		ilon[istn] = argmin(abs.(lon.-mod.(infody[istn,2],360)))
		ilat[istn] = argmin(abs.(lat.-infody[istn,3]))
	end
	md"Finding nearest longitude/latitude coordinate points"
end

# ╔═╡ 60d01a60-f52c-4266-816a-959d73be8e4c
begin
	slon = zeros(4); arglon = zeros(Int16,4)
	slat = zeros(4); arglat = zeros(Int16,4)
	stn  = Vector{String}(undef,4)
	for ii = 1 : length(gdf)
	
		idf = gdf[ii]
		stn[ii]  = idf[:,:"Station"][1]
		slon[ii] = idf[:,:"Longitude"][1]
		slat[ii] = idf[:,:"Latitude"][1]
	
	end
	
	for iarg = 1 : 4
		arglon[iarg] = argmin(abs.(lon.-mod.(slon[iarg],360)))
		arglat[iarg] = argmin(abs.(lat.-slat[iarg]))
	end
	md"Loading positions of the stations in the grid"
end

# ╔═╡ 2d5898e9-93a5-4712-9da4-4b39eccc1435
md"
### C. Extracting the station-point data
"

# ╔═╡ f5bf67c1-745c-41a4-8309-06b0c4cb73cf
begin
	dtvec = npd.start : Day(1) : npd.stop
	dtgpm = npd.start : Day(1) : (npd.stop+Day(1))
	dtnc  = npd.start : Month(1) : npd.stop
	prcp  = zeros(48,length(dtgpm),nstn)
	prcpg = zeros(length(dtvec),4)
	prcps = zeros(length(dtvec),4)
	d18Oμ = zeros(length(dtvec),4)
	d18Oσ = zeros(length(dtvec),4)
	d2Hμ  = zeros(length(dtvec),4)
	d2Hσ  = zeros(length(dtvec),4)
	md"Preallocating data arrays for station GPM data ..."
end

# ╔═╡ 7a4f0cb7-ad33-4c0c-8fb3-8e02e1609021
begin
	for it = 1 : length(dtgpm)
		idt = dtgpm[it]
		ifo = joinpath(npd.datapath,geo.ID,Dates.format(idt,dateformat"yyyy/mm"))
		inc = npd.ID * "-" * geo.ID * "-" * Dates.format(idt,dateformat"yyyymmdd") * ".nc"
		ids = NCDataset(joinpath(ifo,inc))
		for istn = 1 : 4
			prcp[:,it,istn] .= ids["precipitation"][arglon[istn],arglat[istn],:] * 86400
		end
	end
	npdata = reshape(prcp,:,nstn)
	npdata = npdata[11:(end-38),:]
	npdata = dropdims(mean(reshape(npdata,48,:,nstn),dims=1),dims=1)
	md"Loading precipitation data ..."
end

# ╔═╡ 65a634c9-70bd-4166-bc3e-842575dfc85e
begin
	for iarg in 1 : 4
		ii = 0
		dtbeg = gdf[iarg][:,:"Collection Date/Time (Start)"]
		dtend = gdf[iarg][:,:"Collection Date/Time (end)"]
		iprcp = gdf[iarg][:,:"Totalizer precipitation (mm)"]
		
		for ip = 1 : length(iprcp)
			if typeof(iprcp[ip]) <: String
				iprcp[ip] = parse(Float64,iprcp[ip])
			elseif ismissing(iprcp[ip])
				iprcp[ip] = NaN
			end
		end
		
		iprcp = iprcp ./ Dates.value.(dtend .- dtbeg)
		
		for idt in dtvec
			ii += 1
			ind = (idt .>= dtbeg) .& (idt .<= (dtend .- Day(1)))
			if !iszero(sum(ind)) && isone(sum(ind))
				ind2 = findfirst((idt .>= dtbeg) .& (idt .<= (dtend .- Day(1))))
				prcps[ii,iarg] = iprcp[ind][1]
				ibeg = findfirst(dtbeg[ind] .== dtvec)
				iend = findfirst((dtend[ind] .- Day(1)) .== dtvec)
				if isnothing(ibeg); ibeg = 1 end
				if !isnothing(iend)
					prcpg[ii,iarg] = sum(npdata[ibeg:iend,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				else
					prcpg[ii,iarg] = sum(npdata[ibeg:end,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				end
				d18Oμ[ii,iarg] = gdf[iarg][:,:"δ18O, in ‰"][ind][1]
				d18Oσ[ii,iarg] = gdf[iarg][:,:"δ18O Stdev"][ind][1]
				d2Hμ[ii,iarg]  = gdf[iarg][:,:"δ2H, in ‰"][ind][1]
				d2Hσ[ii,iarg]  = gdf[iarg][:,:"δ2H Stdev"][ind][1]
			else
				prcps[ii,iarg] = NaN
				prcpg[ii,iarg] = NaN
				d18Oμ[ii,iarg] = NaN
				d18Oσ[ii,iarg] = NaN
				d2Hμ[ii,iarg]  = NaN
				d2Hσ[ii,iarg]  = NaN
			end
		end
	end
	md"Loading Station data into preallocated arrays ..."
end

# ╔═╡ 031de062-8269-457e-8e07-a83589db7e5a
function smooth(data)
	
	data = (data[1:(end-4)] + data[2:(end-3)] + data[3:(end-2)] + data[4:(end-1)] + data[5:end]) / 5
	
end

# ╔═╡ 415ed9a6-b5bd-4423-b8aa-00cf28952fa6
begin
	arr1 = [[1,1,1,2],[3,3,3,4],[5,5,5,6],[7,7,7,8]]
	pplt.close(); fig,axs = pplt.subplots(
		arr1,aspect=3.2,axwidth=3.5,
		sharex=0,sharey=0,wspace=1.5,hspace=1.5
	)
	
	axs[1].scatter(
		dtvec,prcpg[:,1],s=5,
		label="GPM",legend="ur",legend_kw=Dict("frame"=>false,"ncol"=>1)
	)
	axs[1].scatter(dtvec,prcps[:,1],s=5,label="Station",legend="ur")
	axs[2].scatter(prcps[:,1],prcpg[:,1],s=5)
	axs[3].scatter(dtvec,prcpg[:,2],s=5)
	axs[3].scatter(dtvec,prcps[:,2],s=5)
	axs[4].scatter(prcps[:,2],prcpg[:,2],s=5)
	axs[5].scatter(dtvec,prcpg[:,3],s=5)
	axs[5].scatter(dtvec,prcps[:,3],s=5)
	axs[6].scatter(prcps[:,3],prcpg[:,3],s=5)
	axs[7].scatter(dtvec,prcpg[:,4],s=5)
	axs[7].scatter(dtvec,prcps[:,4],s=5)
	axs[8].scatter(prcps[:,4],prcpg[:,4],s=5)
	
	for ii = [2,4,6,8]
		axs[ii].plot([0,100],[0,100],lw=0.5,c="k",linestyle=":")
		axs[ii].format(
			xlim=(0,100),ylim=(0,100),xlocator=0:20:100,
			ytickloc="right",ylabel=L"GPM Rainfall / mm $^{-1}$"
		)
	end

	for ii = [1,3,5,7]
		axs[ii].format(
			xlim=(npd.start-Month(1),npd.stop+Month(1)),
			ylabel=L"Rainfall Rate / mm $^{-1}$"
		)
	end
	
	for ii = 1 : 6
		axs[ii].format(xticklabels=[])
	end
	
	for ax in axs
		ax.format(ylim=(0,100),ylocator=0:20:100)
	end
	
	axs[1].format(ultitle="(a) $(stn[1])")
	axs[3].format(ultitle="(b) $(stn[2])")
	axs[5].format(ultitle="(c) $(stn[3])")
	axs[7].format(ultitle="(d) $(stn[4])")
	axs[8].format(xlabel=L"Station Rainfall / mm $^{-1}$")
	
	fig.savefig(plotsdir("03b-stationvsgpm.png"),transparent=false,dpi=400)
	load(plotsdir("03b-stationvsgpm.png"))
end

# ╔═╡ d5f4aa86-894d-4bbb-83e2-e053ed08c034
begin
	if isfile(datadir("processed.nc"))
		rm(datadir("processed.nc"),force=true)
	end
	
	nds = NCDataset(datadir("processed.nc"),"c")
	nds.dim["station"] = 4
	nds.dim["time"] = size(npdata,1)
	
	nclon = defVar(nds,"longitude",Float32,("station",),attrib = Dict(
		"units"     => "degrees_east",
		"long_name" => "longitude",
	))
	
	nclat = defVar(nds,"latitude",Float32,("station",),attrib = Dict(
		"units"     => "degrees_north",
		"long_name" => "latitude",
	))
	
	ncprcp = defVar(nds,"prcp",Float32,("time","station",),attrib = Dict(
		"units"     => "mm day**-1",
		"long_name" => "Raw GPM Precipitation Rate",
	))
	
	ncgpm = defVar(nds,"prcpg",Float32,("time","station",),attrib = Dict(
		"units"     => "mm day**-1",
		"long_name" => "Processed GPM Precipitation Rate",
	))
	
	ncstn = defVar(nds,"prcps",Float32,("time","station",),attrib = Dict(
		"units"     => "mm day**-1",
		"long_name" => "Station Precipitation Rate",
	))
	
	ncδ18Oμ = defVar(nds,"δ18Oμ",Float32,("time","station",),attrib = Dict(
		"units"     => "‰",
		"long_name" => "Mean δ18O",
	))
	
	ncδ18Oσ = defVar(nds,"δ18Oσ",Float32,("time","station",),attrib = Dict(
		"units"     => "‰",
		"long_name" => "Standard Deviation δ18O",
	))
	
	ncδ2Hμ = defVar(nds,"δ2Hμ",Float32,("time","station",),attrib = Dict(
		"units"     => "‰",
		"long_name" => "Mean δ2H",
	))
	
	ncδ2Hσ = defVar(nds,"δ2Hσ",Float32,("time","station",),attrib = Dict(
		"units"     => "‰",
		"long_name" => "Standard Deviation δ2H",
	))
	
	nctime = defVar(nds,"time",Int32,("time",),attrib = Dict(
		"units"     => "Days since 2019-07-01 00:00:00.0",
		"long_name" => "time",
		"calendar"  => "gregorian",
	))
	
	nclon[:]   = slon
	nclat[:]   = slat
	ncprcp[:]  = npdata
	ncgpm[:]   = prcpg
	ncstn[:]   = prcps
	ncδ18Oμ[:] = d18Oμ
	ncδ18Oσ[:] = d18Oσ
	ncδ2Hμ[:]  = d2Hμ
	ncδ2Hσ[:]  = d2Hσ
	nctime[:]  = 0 : (length(dtvec) - 1)
	
	close(nds)
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─5b392d8a-d68b-44df-ab57-d0e81aedd668
# ╟─0c9b54e3-78a4-400a-9323-b716d31d327c
# ╟─307f237d-9ba7-49a8-9107-fc000d1aa147
# ╟─6bfec0ba-ac32-4874-be90-f91815c9f188
# ╟─bbfc6d8c-e49e-4444-beb2-7d10b907aa6b
# ╟─1d31fb5c-6af3-4dec-a2f3-fdc91a5e3856
# ╟─7cd13617-40e2-4bc5-ad0a-47b47581ec86
# ╟─66b70aef-f079-4e37-9214-bb1f31ca227b
# ╟─60d01a60-f52c-4266-816a-959d73be8e4c
# ╟─2d5898e9-93a5-4712-9da4-4b39eccc1435
# ╠═f5bf67c1-745c-41a4-8309-06b0c4cb73cf
# ╟─7a4f0cb7-ad33-4c0c-8fb3-8e02e1609021
# ╟─65a634c9-70bd-4166-bc3e-842575dfc85e
# ╟─031de062-8269-457e-8e07-a83589db7e5a
# ╟─415ed9a6-b5bd-4423-b8aa-00cf28952fa6
# ╟─d5f4aa86-894d-4bbb-83e2-e053ed08c034
