### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 90f9e8ea-579b-434e-99eb-04dfe210f1b6
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 6809a6f8-b758-4464-9ea2-5d212109e0cf
begin
	@quickactivate "ColumbiaIsotope"
	using Dates
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ 33e7e5cc-54b8-11ec-396b-c74ec0257d4b
md"
# 03a. Seasonal Isotopic Variability

In this notebook, we investigate if there is noticeable seasonal variability in the isotopic data that we have previously processed and explored in notebook `02c-prcpvsisotopes`.
"

# ╔═╡ 3db056fa-6473-4c84-90ae-c23fefd1c5d0
md"
### A. Loading the Daily Isotope and GPM data
"

# ╔═╡ 5cc36809-413b-46b9-9d58-68777a054e64
pthres = 2.4

# ╔═╡ 7ca0b4fb-6539-4f0c-8de1-8f388cd6ba87
begin
	ds = NCDataset(datadir("processed.nc"))
	prcp  = ds["prcp"][:]
	prcps = ds["prcps"][:]
	prcpg = ds["prcpg"][:]
	δ18Oμ = ds["δ18Oμ"][:]
	δ18Oσ = ds["δ18Oσ"][:]
	δ2Hμ  = ds["δ2Hμ"][:]
	δ2Hσ  = ds["δ2Hσ"][:]
	dtvec = ds["time"][:]; ndt = length(dtvec)
	close(ds)
	
	isless24 = prcps .< pthres
	prcps[isless24,:] .= NaN
	prcpg[isless24,:] .= NaN
	δ18Oμ[isless24,:] .= NaN
	δ18Oσ[isless24,:] .= NaN
	δ2Hμ[isless24,:]  .= NaN
	δ2Hσ[isless24,:]  .= NaN
	md"Loading processed data, and if precipitation is less than $pthres mm/day, set to NaN"
end

# ╔═╡ Cell order:
# ╟─33e7e5cc-54b8-11ec-396b-c74ec0257d4b
# ╟─90f9e8ea-579b-434e-99eb-04dfe210f1b6
# ╟─6809a6f8-b758-4464-9ea2-5d212109e0cf
# ╟─3db056fa-6473-4c84-90ae-c23fefd1c5d0
# ╠═5cc36809-413b-46b9-9d58-68777a054e64
# ╠═7ca0b4fb-6539-4f0c-8de1-8f388cd6ba87
