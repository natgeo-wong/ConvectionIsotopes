### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e32a00ee-5f32-47a1-a983-91fb77bc5d18
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using ERA5Reanalysis
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 03b. Station W-Weighted Pressure, $\sigma$
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 75925166-d1e8-450e-bae6-29affd25d635
md"
### A. Defining Datasets, Variables and Regions
"

# ╔═╡ d3f4a6a9-0915-470b-83af-e05e00dff5d2
npd = IMERGMonthly(start=Date(2013),stop=Date(2020),path=datadir())

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Daily(start=Date(2019,2,1),stop=Date(2021,6,30),path=datadir())

# ╔═╡ 9473b30b-787a-495b-bdfe-b386c9995761
evar_p = SingleVariable("p_wwgt")

# ╔═╡ c7c3ee50-92eb-4848-8b1b-1422559a44dd
evar_σ = SingleVariable("σ_wwgt")

# ╔═╡ 3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
egeo = ERA5Region("OTREC")

# ╔═╡ cb7c6118-e25b-462a-84a5-751ea0682b52
elsd = getLandSea(e5ds,egeo,smooth=true,σlon=2,σlat=2)

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ 59c930cd-5b7f-4047-8660-615148d1bd9f
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ f5f20a42-f325-4c57-8366-3f42cc87f23c
istn = 2

# ╔═╡ c01a9836-fea7-4251-9b8a-ff1e3f5f464c
begin
	ilon = argmin(abs.(elsd.lon.-infody[istn,2]))
	ilat = argmin(abs.(elsd.lat.-infody[istn,3]))
end

# ╔═╡ 1986bfac-7020-4a6e-b480-e797e239f5a3
begin
	ndt = length(e5ds.start : Day(1) : e5ds.stop)
	wp7 = zeros(ndt); wp30 = zeros(ndt)
	wσ7 = zeros(ndt); wσ30 = zeros(ndt)
	for idt = e5ds.start : Month(1) : e5ds.stop
		ii = Dates.value(idt-e5ds.start) .+ (1:daysinmonth(idt))
		ds = read(e5ds,evar_p,egeo,idt,smooth=true,smoothtime=7,quiet=true)
		wp7[ii] .= nomissing(ds[evar_p.ID][ilon,ilat,:],NaN) ./ 100
		close(ds)
		ds = read(e5ds,evar_σ,egeo,idt,smooth=true,smoothtime=7,quiet=true)
		wσ7[ii] .= nomissing(ds[evar_σ.ID][ilon,ilat,:],NaN)
		close(ds)
		ds = read(e5ds,evar_p,egeo,idt,smooth=true,smoothtime=30,quiet=true)
		wp30[ii] .= nomissing(ds[evar_p.ID][ilon,ilat,:],NaN) ./ 100
		close(ds)
		ds = read(e5ds,evar_σ,egeo,idt,smooth=true,smoothtime=30,quiet=true)
		wσ30[ii] .= nomissing(ds[evar_σ.ID][ilon,ilat,:],NaN)
		close(ds)
	end
end

# ╔═╡ d8843fd3-31b9-4014-abf3-663ca330270f
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=2,aspect=3,axwidth=4)
	
	axs[1].plot(e5ds.start : Day(1) : e5ds.stop,wp7)
	axs[1].plot(e5ds.start : Day(1) : e5ds.stop,wp30)
	axs[1].format(ylim=(1000,0))
	
	axs[2].plot(e5ds.start : Day(1) : e5ds.stop,wσ7)
	axs[2].plot(e5ds.start : Day(1) : e5ds.stop,wσ30)
	axs[2].format(ylim=(1,0))
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─d3f4a6a9-0915-470b-83af-e05e00dff5d2
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─9473b30b-787a-495b-bdfe-b386c9995761
# ╟─c7c3ee50-92eb-4848-8b1b-1422559a44dd
# ╟─3d36cdbe-6ef1-4350-bfc8-9e27b1654bff
# ╟─cb7c6118-e25b-462a-84a5-751ea0682b52
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╠═f5f20a42-f325-4c57-8366-3f42cc87f23c
# ╟─c01a9836-fea7-4251-9b8a-ff1e3f5f464c
# ╟─1986bfac-7020-4a6e-b480-e797e239f5a3
# ╠═d8843fd3-31b9-4014-abf3-663ca330270f
