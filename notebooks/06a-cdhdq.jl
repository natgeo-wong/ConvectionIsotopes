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

# ╔═╡ ecdc9af4-9995-4803-be99-ee21bcef32da
geo = GeoRegion("OTREC_wrf_ITCZ05",path=srcdir())

# ╔═╡ af4914a2-597f-4f26-96d0-94f32434b733
begin
	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-cdhodq-daily-20190801_20201231-smooth_07days.nc"))
	time = ds["time"][:]; nt = length(time)
	cHDO = ds["cdHDOdH2O"][:,:]
	cO18 = ds["cdO18dH2O"][:,:]
	close(ds)

	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-dhodq-daily-20190801_20201231-smooth_07days.nc"))
	p = ds["P"][:,:] ./ 100
	dhdq = ds["dHDOdH2O"][:,:]
	dodq = ds["dO18dH2O"][:,:]
	close(ds)
end

# ╔═╡ 633b254b-482c-440f-8920-ffad264619c8
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,sharey=0)
	
	axs[1].scatter(cHDO[1,:],cHDO[2,:]*1e4,s=2)
	axs[2].scatter(cO18[1,:],cO18[2,:]*1e4,s=2)
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ f18abedb-62f5-4bd0-b3c7-403a1342018a
begin
	pvec = 200 : 20 : 1000; psmp = (pvec[1:(end-1)] .+ pvec[2:end])./2;
	np = length(psmp)
	pmat = ones(np,nt) .* psmp
	dHDOdH2O = zeros(np,nt)
	dO18dH2O = zeros(np,nt)
	dhdqbin = 0 : 0.01 : 1.2
	dodqbin = 0.9 : 0.001 : 1.02
	for it = 1 : nt
		dHDOdH2O[:,it] = cHDO[1,it] .+ psmp * cHDO[2,it] * 1e2
		dO18dH2O[:,it] = cO18[1,it] .+ psmp * cO18[2,it] * 1e2
	end
end

# ╔═╡ 74b2fb7e-100f-4b4c-a5bf-8e28d0331c30
begin
	hHDO = Float64.(fit(Histogram,(dHDOdH2O[:],pmat[:],),(dhdqbin,pvec)).weights)
	hHDO = hHDO ./ sum(hHDO) * 100
	md"Binning HDO"
end

# ╔═╡ 6c6c89c9-2ce0-4ebc-a2c5-97bd389d1fb5
begin
	hO18 = Float64.(fit(Histogram,(dO18dH2O[:],pmat[:]),(dodqbin,pvec)).weights)
	hO18 = hO18 ./ sum(hO18) * 100
	md"Binning O18"
end

# ╔═╡ 95aff4f0-690f-4c21-8e68-62103d6e9fa1
begin
	hdhdq = Float64.(fit(Histogram,(dhdq[:],p[:],),(dhdqbin,200:50:1000)).weights)
	hdhdq = hdhdq ./ sum(hdhdq) * 100
	md"Binning wrf HDO"
end

# ╔═╡ a88bd578-f863-4852-9917-9ed559f97dea
begin
	hdodq = Float64.(fit(Histogram,(dodq[:],p[:],),(dodqbin,200:50:1000)).weights)
	hdodq = hdodq ./ sum(hdodq) * 100
	md"Binning wrf O18"
end

# ╔═╡ 99ce6a37-1f32-4d5a-bc11-5b645777d41a
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,aspect=0.5,axwidth=0.5,wspace=0)

	hHDO[iszero.(hHDO)] .= NaN
	hdhdq[iszero.(hdhdq)] .= NaN
	hO18[iszero.(hO18)] .= NaN
	hdodq[iszero.(hdodq)] .= NaN
	
	hHDO[hHDO.<0.1] .= NaN
	hdhdq[hdhdq.<0.1] .= NaN
	hO18[hHDO.<0.1] .= NaN
	hdodq[hdodq.<0.1] .= NaN
	
	
	c = 
	a2[1].pcolormesh(dhdqbin,200:50:1000,hdhdq',levels=(1:15)./10,extend="both")
	a2[2].pcolormesh(dodqbin,200:50:1000,hdodq',levels=(1:15)./10,extend="both")
	c = 
	a2[1].pcolormesh(dhdqbin,pvec,hHDO',levels=(1:15)./10,extend="both",cmap="greys",alpha=0.8)
	a2[2].pcolormesh(dodqbin,pvec,hO18',levels=(1:15)./10,extend="both",cmap="greys",alpha=0.8)

	a2[1].format(xlim=(0.5,1.1))
	a2[2].format(xlim=(0.92,1.01))
	for ax in a2
		ax.format(ylim=(1000,100))
	end

	f2.colorbar(c)
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╠═ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╠═ecdc9af4-9995-4803-be99-ee21bcef32da
# ╠═af4914a2-597f-4f26-96d0-94f32434b733
# ╠═633b254b-482c-440f-8920-ffad264619c8
# ╟─f18abedb-62f5-4bd0-b3c7-403a1342018a
# ╟─74b2fb7e-100f-4b4c-a5bf-8e28d0331c30
# ╟─6c6c89c9-2ce0-4ebc-a2c5-97bd389d1fb5
# ╟─95aff4f0-690f-4c21-8e68-62103d6e9fa1
# ╟─a88bd578-f863-4852-9917-9ed559f97dea
# ╠═99ce6a37-1f32-4d5a-bc11-5b645777d41a
