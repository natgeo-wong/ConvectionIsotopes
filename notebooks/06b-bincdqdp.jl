### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 1d7446ba-f464-44df-89e2-ae2a5726e849
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ab294fae-101f-4587-a2f4-7d72254dd421
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01d. Creating GeoRegions for Stations

In this notebook, we define additional GeoRegions of interest for plotting and for analysis based on WRF modelling output and as necessary for figures.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ 582476a4-e280-458b-ac4c-7c681ff96a74
function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

# ╔═╡ 6aff97ec-0bd3-4d84-9d5c-93393941ca4e
function smooth(data::AbstractVector,days)

	buffer,weights = calculatebufferweights(days)

	ndt = length(data)
	ndata = fill(NaN,ndt)
	smth  = zeros(1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		for ismth = 0 : (buffer*2)
			smth[ismth+1] = data[ii+ismth-buffer] * weights[ismth+1]
		end
		ndata[ii] = sum(smth)

	end

	return ndata

end

# ╔═╡ 10d1c691-00a7-47de-a8ca-8debcd3346c1
function extractbudget(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-∇decompose-20190801_20201231.nc"
	))
	advc = smooth(dropdims(mean(reshape(ds["ADV"][:],24,:),dims=1),dims=1),days) * 86400
	divg = smooth(dropdims(mean(reshape(ds["DIV"][:],24,:),dims=1),dims=1),days) * 86400
	close(ds)

	dsp = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-p_wwgt-$dystr.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	return pwgt,prcp,evap,advc,divg
	
end

# ╔═╡ d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf3","processed","$geoname-cdhodq-$(dystr).nc"))
	c1HDO = ds["cdHDOdH2O"][1,:]
	c2HDO = ds["cdHDOdH2O"][2,:] * 1e4
	c1O18 = ds["cdO18dH2O"][1,:]
	c2O18 = ds["cdO18dH2O"][2,:] * 1e4
	close(ds)

	return c1HDO,c2HDO,c1O18,c2O18
	
end

# ╔═╡ c793412d-71b6-4f2c-a9f4-15da6ec039e4
function plotcdqdp(
	axes,ii;
	nID, days=0, prfx = "", cinfo = false
)

	c1HDOedge = 0.25 : 0.01 : 1.1
	c2HDOedge = -0.02 : 0.001 : 0.1
	c1O18edge = 0.9 : 0.001 : 1.02
	c2O18edge = -0.005 : 0.0002 : 0.015
	binc1HDO = zeros(length(c1HDOedge)-1,nID)
	binc2HDO = zeros(length(c2HDOedge)-1,nID)
	binc1O18 = zeros(length(c1O18edge)-1,nID)
	binc2O18 = zeros(length(c2O18edge)-1,nID)
	c1HDOplt = (c1HDOedge[1:(end-1)] .+ c1HDOedge[2:end])/2
	c2HDOplt = (c2HDOedge[1:(end-1)] .+ c2HDOedge[2:end])/2
	c1O18plt = (c1O18edge[1:(end-1)] .+ c1O18edge[2:end])/2
	c2O18plt = (c2O18edge[1:(end-1)] .+ c2O18edge[2:end])/2
	μc1HDO = zeros(nID); σc1HDO = zeros(nID)
	μc2HDO = zeros(nID); σc2HDO = zeros(nID)
	μc1O18 = zeros(nID); σc1O18 = zeros(nID)
	μc2O18 = zeros(nID); σc2O18 = zeros(nID)
	IDplt = 1 : nID

	for stn in 1 : nID
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_$(prfx)$stnstr"
		c1HDO,c2HDO,c1O18,c2O18 = extract(geoname,days)
		pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt))
		binc1HDO[:,stn] += fit(Histogram,c1HDO[it],c1HDOedge).weights
		binc2HDO[:,stn] += fit(Histogram,c2HDO[it],c2HDOedge).weights
		binc1O18[:,stn] += fit(Histogram,c1O18[it],c1O18edge).weights
		binc2O18[:,stn] += fit(Histogram,c2O18[it],c2O18edge).weights
		μc1HDO[stn] = mean(c1HDO[it]); σc1HDO[stn] = std(c1HDO[it])
		μc2HDO[stn] = mean(c2HDO[it]); σc2HDO[stn] = std(c2HDO[it])
		μc1O18[stn] = mean(c1O18[it]); σc1O18[stn] = std(c1O18[it])
		μc2O18[stn] = mean(c2O18[it]); σc2O18[stn] = std(c2O18[it])
	end

	lvls = 2:2:20
	c1 = 
	axes[4*ii-3].pcolormesh(IDplt,c1HDOplt,binc1HDO,extend="both",levels=lvls)
	c2 = 
	axes[4*ii-2].pcolormesh(IDplt,c2HDOplt,binc2HDO,extend="both",levels=lvls)
	c3 = 
	axes[4*ii-1].pcolormesh(IDplt,c1O18plt,binc1O18,extend="both",levels=lvls)
	c4 = 
	axes[4*ii-0].pcolormesh(IDplt,c2O18plt,binc2O18,extend="both",levels=lvls)

	
	axes[4*ii-3].plot(IDplt,μc1HDO)
	axes[4*ii-2].plot(IDplt,μc2HDO)
	axes[4*ii-1].plot(IDplt,μc1O18)
	axes[4*ii-0].plot(IDplt,μc2O18)
	# axes[4*ii-3].errorbar(IDplt,μc1HDO,σc1HDO)
	# axes[4*ii-2].errorbar(IDplt,μc2HDO,σc2HDO)
	# axes[4*ii-1].errorbar(IDplt,μc1O18,σc1O18)
	# axes[4*ii-0].errorbar(IDplt,μc2O18,σc2O18)

	if cinfo
		return c1,c2,c3,c4
	else
		return
	end

end

# ╔═╡ 64ed26fa-fd4e-4ac7-a8d1-741dc40d91a7
function axesformat!(axes)

	

	return

end

# ╔═╡ 30424aa0-cc38-4f50-8eb6-efd4f6c4c9c4
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=4,aspect=2,axwidth=2,)

	c1,_,_,_ = 
	plotcdqdp(axs,1,nID=25,prfx="ITCZ",days=30,cinfo=true)
	
	axesformat!(axs)

	fig.colorbar(c1,loc="b",label="Number of Observations")
	fig.savefig(plotsdir("05g-pvdhqdp-smooth07.png"),transparent=false,dpi=400)
	load(plotsdir("05g-pvdhqdp-smooth07.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─582476a4-e280-458b-ac4c-7c681ff96a74
# ╟─6aff97ec-0bd3-4d84-9d5c-93393941ca4e
# ╟─10d1c691-00a7-47de-a8ca-8debcd3346c1
# ╟─d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
# ╠═c793412d-71b6-4f2c-a9f4-15da6ec039e4
# ╠═64ed26fa-fd4e-4ac7-a8d1-741dc40d91a7
# ╠═30424aa0-cc38-4f50-8eb6-efd4f6c4c9c4
