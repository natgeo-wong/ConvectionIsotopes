### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 02c9bf80-07c6-411a-884c-d7b275583fec
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ d8b1349a-e91a-48f5-9c58-8e34a6a4db56
begin
	@quickactivate "ColombiaIsotope"
	using DataFrames
	using DelimitedFiles
	using GLM
	using Printf
	using Statistics
	using XLSX
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ b4e60800-6d92-11ec-0046-fb8e0610278d
md"
# 03c. Looking at Station Isotopic Data

In this notebook, we look at the monthly and daily isotopic data for each station and compare them against each other.
"

# ╔═╡ 7bf91579-8bf2-41ef-aab6-7f9da1d270b7
begin
	mdf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Monthly"))
	mdf = mdf[:,[:"Station",:"δ2H, in ‰",:"δ2H Stdev",:"δ18O, in ‰",:"δ18O Stdev",]]
	ddf = DataFrame(XLSX.readtable(datadir("IsotopeDataSummary.xlsx"),"Daily"))
	ddf = ddf[:,[:"Station",:"δ2H, in ‰",:"δ2H Stdev",:"δ18O, in ‰",:"δ18O Stdev",]]
	md"Loading Isotope Data information ..."
end

# ╔═╡ bac5aa56-bb3a-408f-9c29-f215a37fa93c
begin
	δ2Hμ_mo  = mdf[:,:"δ2H, in ‰"]
	δ18Oμ_mo = mdf[:,:"δ18O, in ‰"]
	ind_mo = .!ismissing.(δ18Oμ_mo) .& .!ismissing.(δ2Hμ_mo)
	δ2Hμ_mo  = δ2Hμ_mo[ind_mo]
	δ18Oμ_mo = δ18Oμ_mo[ind_mo]
	md"Removing missing data from monthly data ..."
end

# ╔═╡ f7ccbe2d-a4cc-4066-9a1d-09c4491827b7
begin
	δ2Hμ_dy  = ddf[:,:"δ2H, in ‰"]
	δ18Oμ_dy = ddf[:,:"δ18O, in ‰"]
	ind_dy = .!ismissing.(δ18Oμ_dy) .& .!ismissing.(δ2Hμ_dy)
	δ2Hμ_dy  = δ2Hμ_dy[ind_dy]
	δ18Oμ_dy = δ18Oμ_dy[ind_dy]
	md"Removing missing data from daily data ..."
end

# ╔═╡ 514b0cdf-1776-41ec-8d09-61b2b8d080fa
begin
	ddf2 = DataFrame(H=Float64.(δ2Hμ_dy),O=Float64.(δ18Oμ_dy))
	mdy  = lm(@formula(H ~ O), ddf2)
	c_dy,m_dy = coef(mdy.model)
	mdf2 = DataFrame(H=Float64.(δ2Hμ_mo),O=Float64.(δ18Oμ_mo))
	mmo  = lm(@formula(H ~ O), mdf2)
	c_mo,m_mo = coef(mmo.model)
	md"Performing linear regression on monthly and daily data ..."
end

# ╔═╡ f5b6f440-2c9f-4815-9f0a-a8bd71f32806
begin
	mgdf = groupby(mdf,"Station"); nmostn = length(mgdf)
	dgdf = groupby(ddf,"Station"); ndystn = length(dgdf)
	clr_mo = pplt.Colors("Phase",nmostn)
	clr_dy = pplt.Colors("Twilight",ndystn)
	md"Loading colors and splitting dataframes by station ..."
end

# ╔═╡ 0caadf1d-d68d-41dd-a589-7c51e640ad96
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=2,axwidth=2)

	for istn = 1 : nmostn
		δ2Hμ  = mgdf[istn][:,:"δ2H, in ‰"]
		δ18Oμ = mgdf[istn][:,:"δ18O, in ‰"]
		ind = .!ismissing.(δ18Oμ) .& .!ismissing.(δ2Hμ)
		δ2Hμ  = δ2Hμ[ind]
		δ18Oμ = δ18Oμ[ind]
		a1[1].scatter(δ18Oμ,δ2Hμ,s=1,c=clr_mo[istn])
	end

	for istn = 1 : ndystn
		δ2Hμ  = dgdf[istn][:,:"δ2H, in ‰"]
		δ18Oμ = dgdf[istn][:,:"δ18O, in ‰"]
		ind = .!ismissing.(δ18Oμ) .& .!ismissing.(δ2Hμ)
		δ2Hμ  = δ2Hμ[ind]
		δ18Oμ = δ18Oμ[ind]
		a1[2].scatter(δ18Oμ,δ2Hμ,s=1,c=clr_dy[istn])
	end
	
	a1[1].plot([-25,8],m_mo*[-25,8].+c_mo,c="k",linestyle="--")
	a1[2].plot([-25,8],m_dy*[-25,8].+c_dy,c="k",linestyle="--")

	a1[1].text(-15,-150,"y = " * @sprintf("%4.2f",m_mo) * " x + " * @sprintf("%4.2f",c_mo))
	a1[2].text(-15,-150,"y = " * @sprintf("%4.2f",m_dy) * " x + " * @sprintf("%4.2f",c_dy))

	a1[1].format(ultitle="(a) Monthly Data")
	a1[2].format(ultitle="(b) Daily Data")
	for ax in a1
		ax.format(
			xlim=(-25,8),xlabel=L"$\delta^{18}$O",
			ylim=(-175,50),ylabel=L"$\delta^{2}$H"
		)
	end
	
	f1.savefig(plotsdir("03c-stationisotopes.png"),transparent=false,dpi=400)
	load(plotsdir("03c-stationisotopes.png"))
end

# ╔═╡ 7bd359fe-7d92-4304-b340-72b159560978
md"
Here, we see and confirm that the daily data has been adjusted such that the $\delta^2$H and $\delta^{18}$O columns are correct.  For both monthly and daily data, we see that indeed the relationship between the two isotopic concentrations in rainwater have a consistent linear relationship that spans across different timescales.  By using different stations for different colors, we see that this relationship is also constant across different station locations.
"

# ╔═╡ Cell order:
# ╟─b4e60800-6d92-11ec-0046-fb8e0610278d
# ╟─02c9bf80-07c6-411a-884c-d7b275583fec
# ╟─d8b1349a-e91a-48f5-9c58-8e34a6a4db56
# ╟─7bf91579-8bf2-41ef-aab6-7f9da1d270b7
# ╟─bac5aa56-bb3a-408f-9c29-f215a37fa93c
# ╟─f7ccbe2d-a4cc-4066-9a1d-09c4491827b7
# ╟─514b0cdf-1776-41ec-8d09-61b2b8d080fa
# ╟─f5b6f440-2c9f-4815-9f0a-a8bd71f32806
# ╟─0caadf1d-d68d-41dd-a589-7c51e640ad96
# ╟─7bd359fe-7d92-4304-b340-72b159560978
