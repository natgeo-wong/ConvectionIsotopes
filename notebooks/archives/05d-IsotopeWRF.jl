### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 3c4b037f-c854-4148-ae6e-fff50bbb51cf
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using NCDatasets
	using NumericalIntegration
	using Statistics
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
md"
# 05d. Isotope vs W-Weighted Pressure
"

# ╔═╡ c853d30e-5a51-47f4-a6be-88ab857253c1
md"
### A. Loading some Sample Data
"

# ╔═╡ 14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
begin
	rain = []
	rHDO = []
	pwgt = []
	for imo = 8 : 12
		ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-RAINNC"))
		irin = ds["RAINNC"][:]
		close(ds)
		ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-HDO_RAINNC"))
		iHDO = ds["HDO_RAINNC"][:]
		iHDO = (iHDO./irin .-1) *1000
		close(ds)
		ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-p_wwgt"))
		ipwg = ds["p_wwgt"][:]
		close(ds)
		rain = vcat(rain,irin[:])
		rHDO = vcat(rHDO,iHDO[:])
		pwgt = vcat(pwgt,ipwg[:])
	end
end

# ╔═╡ 42f8a973-89c4-4a65-b752-a9235d637b21
begin
	rainl = rain[:]
	pwgtl = pwgt[:] ./ 100
	rHDOl = rHDO[:]

	rainb = 5  : 5 : 95; nr = length(rainb)
	rainu = 10 : 5 : 100
	rainm = 5  : 5 : 100
	pwgtb = 100 : 25 : 925; np = length(pwgtb)
	pwgtu = 125 : 25 : 950
	pwgtm = 100 : 25 : 950

	binmat = zeros(nr,np)

	for ip = 1 : np, ir = 1 : nr
		try
			binmat[ir,ip] = mean(rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .&
									(rainl.>rainb[ir]) .& (rainl.<rainu[ir])])
		catch
			binmat[ir,ip] = NaN
		end

	end
end

# ╔═╡ e884be27-b9ec-45ba-b800-afccbfb54a78
begin
	pplt.close(); f4,a4 = pplt.subplots(axwidth=2,aspect=2/3)
	
	c4 = a4[1].pcolormesh(
		rainm,pwgtm,binmat',
		cmap="viridis",levels=-120:5:-50,extend="both"
	)
	a4[1].format(
		ylim=(maximum(pwgtm),minimum(pwgtm)),ylabel=L"p$_w$ / hPa",
		xlim=(minimum(rainm),maximum(rainm)),xlabel=L"Rain Rate / mm day$^{-1}$"
	)

	f4.colorbar(c4,length=0.75,label=L"$\delta^2$H / $\perthousand$")
	f4.savefig(plotsdir("05d-IsotopeWRF.png"),transparent=false,dpi=200)
	load(plotsdir("05d-IsotopeWRF.png"))
end

# ╔═╡ Cell order:
# ╠═c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
# ╟─3c4b037f-c854-4148-ae6e-fff50bbb51cf
# ╟─4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
# ╟─c853d30e-5a51-47f4-a6be-88ab857253c1
# ╟─14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
# ╟─42f8a973-89c4-4a65-b752-a9235d637b21
# ╟─e884be27-b9ec-45ba-b800-afccbfb54a78
