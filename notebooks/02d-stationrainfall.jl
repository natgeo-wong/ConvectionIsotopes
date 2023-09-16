### A Pluto.jl notebook ###
# v0.19.27

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
# 02d. Corroborating Station Rainfall w/ GPM

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ 5b392d8a-d68b-44df-ab57-d0e81aedd668
md"
### A. Defining the NASAPrecipitation Dataset
"

# ╔═╡ 0c9b54e3-78a4-400a-9323-b716d31d327c
npd = IMERGFinalHH(start=Date(2019,1,1),stop=Date(2021,6,30),path=datadir())

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
	df_co = DataFrame(XLSX.readtable(datadir("IsotopeSummary.xlsx"),"ColombiaDaily"))
	df_co = df_co[:,[
		:"Station",
		:"δ2H, in ‰",:"δ2H, stdev",:"δ18O, in ‰",:"δ18O, stdev",
		:"Collection Date/Time (Start)",:"Collection Date/Time (end)",
		:"Totalizer precipitation (mm)"
	]]
	gdf_co = groupby(df_co,"Station")
	md"Loading Colombia Daily Isotope Data information ..."
end

# ╔═╡ e524c41a-681f-47ff-91ab-22a3be2e93cb
begin
	df_cr = DataFrame(XLSX.readtable(datadir("IsotopeSummary.xlsx"),"CostaRica"))
	df_cr = df_cr[:,[
		:"Station",
		:"δ2H, in ‰",:"δ2H, stdev",:"δ18O, in ‰",:"δ18O, stdev",
		:"Collection Date/Time (Start)",:"Collection Date/Time (end)",
		:"Totalizer precipitation (mm)"
	]]
	gdf_cr = groupby(df_cr,"Station")
	md"Loading Costa Rica Daily Isotope Data information ..."
end

# ╔═╡ 66b70aef-f079-4e37-9214-bb1f31ca227b
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ 60d01a60-f52c-4266-816a-959d73be8e4c
begin
	slon = zeros(12); arglon = zeros(Int16,12)
	slat = zeros(12); arglat = zeros(Int16,12)
	stn  = Vector{String}(undef,12)
	for ii = 1 : length(gdf_co)
	
		idf = gdf_co[ii]
		stn[ii]  = idf[:,:"Station"][1]
		slon[ii] = infody[infody[:,1] .== stn[ii],2][1]
		slat[ii] = infody[infody[:,1] .== stn[ii],3][1]
	
	end
	for ii = 1 : length(gdf_cr)
	
		idf = gdf_cr[ii]
		stn[ii+4]  = idf[:,:"Station"][1]
		slon[ii+4] = infody[infody[:,1] .== stn[ii+4],2][1]
		slat[ii+4] = infody[infody[:,1] .== stn[ii+4],3][1]
	
	end
	
	for iarg = 1 : 12
		arglon[iarg] = argmin(abs.(lon.-slon[iarg]))
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
	prcpg = zeros(length(dtvec),12)
	prcps = zeros(length(dtvec),12)
	d18Oμ = zeros(length(dtvec),12)
	d18Oσ = zeros(length(dtvec),12)
	d2Hμ  = zeros(length(dtvec),12)
	d2Hσ  = zeros(length(dtvec),12)
	md"Preallocating data arrays for station GPM data ..."
end

# ╔═╡ 7a4f0cb7-ad33-4c0c-8fb3-8e02e1609021
begin
	for it = 1 : length(dtgpm)
		idt = dtgpm[it]
		ifo = joinpath(npd.datapath,geo.ID,Dates.format(idt,dateformat"yyyy/mm"))
		inc = npd.ID * "-" * geo.ID * "-" * Dates.format(idt,dateformat"yyyymmdd") * ".nc"
		ids = NCDataset(joinpath(ifo,inc))
		for istn = 1 : 12
			prcp[:,it,istn] .= ids["precipitation"][arglon[istn],arglat[istn],:] * 86400
		end
	end
	npdata = reshape(prcp,:,nstn)
	npdata_co = npdata[11:(end-38),:]
	npdata_co = dropdims(mean(reshape(npdata_co,48,:,nstn),dims=1),dims=1)
	npdata_cr = npdata[13:(end-36),:]
	npdata_cr = dropdims(mean(reshape(npdata_cr,48,:,nstn),dims=1),dims=1)
	npdata = cat(npdata_co[:,1:4],npdata_cr[:,5:end],dims=2)
	md"Loading precipitation data ..."
end

# ╔═╡ 65a634c9-70bd-4166-bc3e-842575dfc85e
begin
	for iarg in 1 : 4
		ii = 0
		dtbeg = gdf_co[iarg][:,:"Collection Date/Time (Start)"]
		dtend = gdf_co[iarg][:,:"Collection Date/Time (end)"]
		iprcp = gdf_co[iarg][:,:"Totalizer precipitation (mm)"]
		
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
					prcpg[ii,iarg] = sum(npdata_co[ibeg:iend,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				else
					prcpg[ii,iarg] = sum(npdata_co[ibeg:end,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				end
				d18Oμ[ii,iarg] = gdf_co[iarg][:,:"δ18O, in ‰"][ind][1]
				d18Oσ[ii,iarg] = gdf_co[iarg][:,:"δ18O, stdev"][ind][1]
				d2Hμ[ii,iarg]  = gdf_co[iarg][:,:"δ2H, in ‰"][ind][1]
				d2Hσ[ii,iarg]  = gdf_co[iarg][:,:"δ2H, stdev"][ind][1]
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
	md"Loading Colombia Daily Station data into preallocated arrays ..."
end

# ╔═╡ 93e79d88-e8c3-4eb6-ab2a-e5559b089003
begin
	for iarg in 1 : 8
		ii = 0
		dtbeg = gdf_cr[iarg][:,:"Collection Date/Time (Start)"]
		dtend = gdf_cr[iarg][:,:"Collection Date/Time (end)"]
		iprcp = gdf_cr[iarg][:,:"Totalizer precipitation (mm)"]
		
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
				prcps[ii,iarg+4] = iprcp[ind][1]
				ibeg = findfirst(dtbeg[ind] .== dtvec)
				iend = findfirst((dtend[ind] .- Day(1)) .== dtvec)
				if isnothing(ibeg); ibeg = 1 end
				if !isnothing(iend)
					prcpg[ii,iarg+4] = sum(npdata_cr[ibeg:iend,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				else
					prcpg[ii,iarg+4] = sum(npdata_cr[ibeg:end,iarg]) / Dates.value(dtend[ind2] - dtbeg[ind2])
				end
				d18Oμ[ii,iarg+4] = gdf_cr[iarg][:,:"δ18O, in ‰"][ind][1]
				d18Oσ[ii,iarg+4] = NaN
				d2Hμ[ii,iarg+4]  = gdf_cr[iarg][:,:"δ2H, in ‰"][ind][1]
				d2Hσ[ii,iarg+4]  = NaN
			else
				prcps[ii,iarg+4] = NaN
				prcpg[ii,iarg+4] = NaN
				d18Oμ[ii,iarg+4] = NaN
				d18Oσ[ii,iarg+4] = NaN
				d2Hμ[ii,iarg+4]  = NaN
				d2Hσ[ii,iarg+4]  = NaN
			end
		end
	end
	md"Loading Costa Rica Daily Station data into preallocated arrays ..."
end

# ╔═╡ 031de062-8269-457e-8e07-a83589db7e5a
function smooth(data)
	
	data = (data[1:(end-4)] + data[2:(end-3)] + data[3:(end-2)] + data[4:(end-1)] + data[5:end]) / 5
	
end

# ╔═╡ d5f4aa86-894d-4bbb-83e2-e053ed08c034
begin
	if isfile(datadir("processed.nc"))
		rm(datadir("processed.nc"),force=true)
	end
	
	nds = NCDataset(datadir("processed.nc"),"c")
	nds.dim["station"] = 12
	nds.dim["time"] = size(prcpg,1)
	
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
		"units"     => "Days since $(dtvec[1]) 00:00:00.0",
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

# ╔═╡ 415ed9a6-b5bd-4423-b8aa-00cf28952fa6
begin
	arr1 = [[1,1,1,2],[3,3,3,4],[5,5,5,6],[7,7,7,8]]
	pplt.close(); fig,axs = pplt.subplots(
		arr1,aspect=3.2,axwidth=3.5,
		sharex=0,sharey=0,wspace=1.5,hspace=1.5
	)
	
	axs[1].scatter(
		dtvec,prcpg[:,5],s=5,
		label="GPM",legend="ur",legend_kw=Dict("frame"=>false,"ncol"=>1)
	)
	axs[1].scatter(dtvec,prcps[:,5],s=5,label="Station",legend="ur")
	axs[2].scatter(prcps[:,5],prcpg[:,5],s=5)
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

# ╔═╡ de139a33-a91d-4a6e-bef0-7833671c2a3a
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,axwidth=2,sharey=0,sharex=0)
	
	a2[1].scatter(prcpg,prcps,s=2)
	a2[1].format(xlim=(-10,210),ylim=(-10,210))
	
	a2[2].scatter(d18Oμ,d2Hμ,s=2)
	
	f2.savefig("test.png",transparent=false,dpi=300)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─5b392d8a-d68b-44df-ab57-d0e81aedd668
# ╟─0c9b54e3-78a4-400a-9323-b716d31d327c
# ╟─307f237d-9ba7-49a8-9107-fc000d1aa147
# ╟─6bfec0ba-ac32-4874-be90-f91815c9f188
# ╟─1d31fb5c-6af3-4dec-a2f3-fdc91a5e3856
# ╟─7cd13617-40e2-4bc5-ad0a-47b47581ec86
# ╟─e524c41a-681f-47ff-91ab-22a3be2e93cb
# ╟─66b70aef-f079-4e37-9214-bb1f31ca227b
# ╟─60d01a60-f52c-4266-816a-959d73be8e4c
# ╟─2d5898e9-93a5-4712-9da4-4b39eccc1435
# ╟─f5bf67c1-745c-41a4-8309-06b0c4cb73cf
# ╟─7a4f0cb7-ad33-4c0c-8fb3-8e02e1609021
# ╟─65a634c9-70bd-4166-bc3e-842575dfc85e
# ╟─93e79d88-e8c3-4eb6-ab2a-e5559b089003
# ╟─031de062-8269-457e-8e07-a83589db7e5a
# ╟─d5f4aa86-894d-4bbb-83e2-e053ed08c034
# ╟─415ed9a6-b5bd-4423-b8aa-00cf28952fa6
# ╠═de139a33-a91d-4a6e-bef0-7833671c2a3a
