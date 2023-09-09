### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ fba33adb-d8a9-495e-b927-b9b62aaa74ef
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9149a006-add6-467e-9ad1-be497a53fff7
begin
	@quickactivate "ConvectionIsotopes"
	using DataFrames
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using XLSX

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 8e744f98-21b4-11ed-017a-8571e8807e61
md"
# 05a. GNIP Stations in Southeast Asia
"

# ╔═╡ 3b799f5e-0789-49cb-a1c7-895edbf9c7e5
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	clon = coast[:,1]; clat = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 65059102-8d7c-4481-8931-2603480ba231
md"
### A. Loading GNIP Station Locations
"

# ╔═╡ 9f498a11-07d4-4c13-b932-4679f1056978
begin
	xf = XLSX.readxlsx(datadir("GNIP-SEA.xlsx"))["Data"]["A1:I973"]
	md"Loading GNIP data ..."
end

# ╔═╡ ef55c91a-a9cb-4749-89b1-fe9c7f2fd13c
begin
	ds  = NCDataset(datadir("flsm","flsm-TRP.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	lsm = ds["flsm"][:]
	close(ds)
	md"Loading filtered ERA5 land-sea mask ..."
end

# ╔═╡ cd06ed14-7bfc-4508-a0ef-75a842e65498
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=1,axwidth=3)

	lvls = vcat(10. .^(-4:-1),0.2,0.5,0.9,0.95,0.97,0.98,0.99)
	
	axs[1].plot(clon,clat,c="k",lw=0.5)
	c = axs[1].pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both")
	axs[1].scatter(xf[2:end,4],xf[2:end,3],c="r",s=10)
	axs[1].plot([99,99,105,105,99],[1,7,7,1,1],c="w",linestyle="--")
	axs[1].plot([99,112],[7,19],c="w",linestyle="--")
	axs[1].plot([105,124],[1,7.3],c="w",linestyle="--")
	axs[1].plot([115,119,119,115,115],[3,3,7,7,3],c="w",linestyle="--")
	axs[1].plot([115,116.2],[3,-9],c="w",linestyle="--")
	axs[1].plot([119,124.2],[7,-1.2],c="w",linestyle="--")
	axs[1].format(
		xlim=(95,125),ylim=(-10,20),xlocator=95:5:125,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
	)
	axs[1].colorbar(c,length=0.8)

	ix1 = fig.add_axes([0.51,0.615,0.27,0.33])
	ix1.pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both")
	ix1.scatter(xf[2:end,4],xf[2:end,3],c="r",s=10)
	ix1.plot(clon,clat,c="k",lw=0.5)
	ix1.format(xlim=(99,105),ylim=(1,7),xlocator=[],ylocator=[])

	ix2 = fig.add_axes([0.60,0.16,0.18,0.22])
	ix2.pcolormesh(lon,lat,lsm',levels=lvls,cmap="delta",extend="both")
	ix2.scatter(xf[2:end,4],xf[2:end,3],c="r",s=10)
	ix2.plot(clon,clat,c="k",lw=0.5)
	ix2.format(xlim=(115,119),ylim=(3,7),xlocator=[],ylocator=[])
	
	fig.savefig(plotsdir("05a-GNIPstations.png"),transparent=false,dpi=200)
	load(plotsdir("05a-GNIPstations.png"))
end

# ╔═╡ b597fa66-3849-49fb-bd76-d23e8e46984f
begin
	stlist = xf[2:end,1]
	sites  = unique(stlist)
	nsites = length(sites)
	stinfo = Array{Any,2}(undef,nsites,4)
	
	for isite = 1 : nsites
		ii = findfirst(stlist.==sites[isite])
		stinfo[isite,1] = xf[ii,1]
		stinfo[isite,2] = xf[ii,4]
		stinfo[isite,3] = xf[ii,3]
		stinfo[isite,4] = xf[ii,5]
	end
	open("gnipstninfo.csv","w") do io
        writedlm(io,stinfo,',')
	end
	md"Saving station information into a .csv file ..."
end

# ╔═╡ 3419ccfa-848d-49a5-886e-a283114feb72
md"The cluster of stations in northeast Borneo (on the eastern coast) is a collection of stations at TAWAU.  We assume, for the purposes of this project, that essentially they are the same station."

# ╔═╡ 8dca0258-9a69-4770-bf00-e8919c18e0df
md"
### B. Exploring the vertical-velocity-weighted column pressure in the region
"

# ╔═╡ 2e783000-8fb9-4e48-9471-16ecf4780596
e5ds = ERA5Monthly(start=Date(2013,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 0ffb4e9f-41e4-4714-8b3d-8795f0df9f5c
evar = SingleVariable("p_wwgt")

# ╔═╡ 42381a85-bf1b-4fcf-b50c-2c590e94f99c
geo = GeoRegion("TRP")

# ╔═╡ e160b2ed-5c60-4156-876c-b2a5bfde3631
ereg = ERA5Region(geo)

# ╔═╡ 0b2d990d-5d24-47fb-8110-8e1f81cde7b6
begin
	eds = read(e5ds,evar,ereg,Date(2001))
	var = nomissing(eds[evar.varID][:],0) / 100
	cnt = Int.(.!iszero.(var))
	close(eds)
	for dt in Date(2002) : Year(1) : Date(2020)
		ids = read(e5ds,evar,ereg,dt)
		var[:,:,:] += nomissing(ids[evar.varID][:],0) / 100
		cnt[:,:,:] += Int.(.!iszero.(nomissing(ids[evar.varID][:],0)))
		close(ids)
	end
end

# ╔═╡ 8438af8a-2a54-40cc-aac7-6a3c79178e37
begin
	wwgt_pre = var ./ cnt

	lds = NCDataset(datadir("flsm","flsm-$(ereg.geoID).nc"))
	nlon = length(lon); nlat = length(lat)
	close(lds)

	for ilat = 1 : nlat, ilon = 1 : nlon
		if lsm[ilon,ilat] > 0.9
			wwgt_pre[ilon,ilat,:] .= NaN
		end
	end
	md"Filtering out data based on the filtered land-sea mask ..."
end

# ╔═╡ 9bb240ae-097b-406b-b0e9-486bc1783b7b
begin
	pplt.close(); f2,a2 = pplt.subplots(
		nrows=3,ncols=4,aspect=1,axwidth=1,hspace=1,wspace=1
	)
	
	c2 = a2[1].pcolormesh(
		lon,lat,wwgt_pre[:,:,1]',
		levels=250:50:850,extend="both",cmap="delta"
	)
	for ii = 2 : 12
		a2[ii].pcolormesh(
			lon,lat,wwgt_pre[:,:,ii]',
			levels=250:50:850,extend="both",cmap="delta"
		)
	end
	
	for ax in a2
		ax.scatter(xf[2:end,4],xf[2:end,3],c="r",s=10)
		ax.plot(clon,clat,c="k",lw=0.5)
		ax.format(
			xlim=(95,125),ylim=(-10,20),
			xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
			suptitle="W-weighted Mean Pressure (2001-2020)"
		)
	end

	f2.colorbar(c2,length=0.7,label=L"$p_w$ / hPa")
	f2.savefig(plotsdir("05a-wwgtpre_GNIP.png"),transparent=false,dpi=400)
	load(plotsdir("05a-wwgtpre_GNIP.png"))
end

# ╔═╡ d2001450-fe8b-4412-b2d3-6f659a737f9d
md"
### C. Precipitation Data Comparison
"

# ╔═╡ 95d5d27f-69c2-405b-8376-1d02ae52a5d1
npd = IMERGMonthly(start=Date(2001),stop=Date(2020),path=datadir())

# ╔═╡ b671a99c-9824-4fb8-aca4-1ac35fbeb3eb
npd_lsd = getIMERGlsd(geo,path=datadir("imergmask"))

# ╔═╡ b3f3d984-d9b4-4869-946c-1a683f2e051e
begin
	gnip_prcp = xf[2:end,end]
	gnip_prcp[ismissing.(gnip_prcp)] .= NaN
	gnip_prcp = convert(Vector{Float64},gnip_prcp)
	npd_prcp  = zeros(length(gnip_prcp)) * NaN
end

# ╔═╡ 7025fc11-43b2-46c1-9a30-e4bbda22563b
for idata = 1 : length(npd_prcp)
	idt = Date(xf[idata+1,6])
	if idt >= Date(2001)
		imo = month(idt)
		ilon = argmin(abs.(npd_lsd.lon .- xf[idata+1,4]))
		ilat = argmin(abs.(npd_lsd.lat .- xf[idata+1,3]))
		npdds = read(npd,geo,idt)
		npd_prcp[idata] = npdds["precipitation"][ilon,ilat,imo] * 86400 * daysinmonth(idt)
		close(npdds)
	end
end

# ╔═╡ 83209cf1-3cf0-4dc8-9056-bc0240dc3311
begin
	pplt.close(); f3,a3 = pplt.subplots()
	
	a3[1].scatter(gnip_prcp,npd_prcp,s=5)
	a3[1].plot([0,5000],[0,5000],lw=1,c="k",linestyle="--")
	a3[1].format(
		xlim=(1,5000),xscale="log",xlabel=L"GNIP Precipitation / mm day$^{-1}$",
		ylim=(1,5000),yscale="log",ylabel=L"GPM IMERG Precipitation / mm day$^{-1}$"
	)
	
	f3.savefig(plotsdir("05a-gnipvsgpm.png"),transparent=false,dpi=150)
	load(plotsdir("05a-gnipvsgpm.png"))
end

# ╔═╡ Cell order:
# ╟─8e744f98-21b4-11ed-017a-8571e8807e61
# ╟─fba33adb-d8a9-495e-b927-b9b62aaa74ef
# ╟─9149a006-add6-467e-9ad1-be497a53fff7
# ╟─3b799f5e-0789-49cb-a1c7-895edbf9c7e5
# ╟─65059102-8d7c-4481-8931-2603480ba231
# ╟─9f498a11-07d4-4c13-b932-4679f1056978
# ╟─ef55c91a-a9cb-4749-89b1-fe9c7f2fd13c
# ╟─cd06ed14-7bfc-4508-a0ef-75a842e65498
# ╟─b597fa66-3849-49fb-bd76-d23e8e46984f
# ╟─3419ccfa-848d-49a5-886e-a283114feb72
# ╟─8dca0258-9a69-4770-bf00-e8919c18e0df
# ╟─2e783000-8fb9-4e48-9471-16ecf4780596
# ╟─0ffb4e9f-41e4-4714-8b3d-8795f0df9f5c
# ╠═42381a85-bf1b-4fcf-b50c-2c590e94f99c
# ╟─e160b2ed-5c60-4156-876c-b2a5bfde3631
# ╟─0b2d990d-5d24-47fb-8110-8e1f81cde7b6
# ╟─8438af8a-2a54-40cc-aac7-6a3c79178e37
# ╟─9bb240ae-097b-406b-b0e9-486bc1783b7b
# ╟─d2001450-fe8b-4412-b2d3-6f659a737f9d
# ╟─95d5d27f-69c2-405b-8376-1d02ae52a5d1
# ╟─b671a99c-9824-4fb8-aca4-1ac35fbeb3eb
# ╟─b3f3d984-d9b4-4869-946c-1a683f2e051e
# ╟─7025fc11-43b2-46c1-9a30-e4bbda22563b
# ╟─83209cf1-3cf0-4dc8-9056-bc0240dc3311
