### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ConvectionIsotopes"
	using ERA5Reanalysis
	using DelimitedFiles
	using NCDatasets
	using PlutoUI
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 02b. Vertical Wind Profiles for Weighted Column-Mean $\sigma$

Here, we plot the profiles of vertical velocity for different $\sigma$ in different regions within the defined OTREC region.
"

# ╔═╡ 30d2be4a-9bdd-4939-8785-413cd7a59f78
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 1a69f197-a5a8-476b-86e7-c472061ca1d5
geolist = ["TRP","CIS_WPC","CIS_CPC","OTREC","OTREC_PAC","OTREC_ATL"]

# ╔═╡ d413d3cb-cf76-4aa5-80ff-c929d396929d
load(plotsdir("01b-georegions.png"))

# ╔═╡ 92621f09-ed8c-46aa-afc0-4cab187dde68
clr = pplt.Colors("oslo",12)

# ╔═╡ 4227bb02-d323-4653-a071-279eff3626d9
md"### B. Averaged vertical velocity profiles for each $\sigma_w$ bin"

# ╔═╡ 059b955b-1fa1-4a60-a136-730146cce921
begin
	
	pplt.close()
	arr = [1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7]
	wsp = [0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0]
	fig,axs = pplt.subplots(arr,aspect=1/6,axwidth=0.6,sharex=0)

	for ii = 1 : 6

		gID = geolist[ii]

		ds = NCDataset(datadir("wprofile","wwgtpre-profile-$(gID)x0.25.nc"))
		σ  = ds["σ"][:]
		σw = ds["σ_wgt"][:]; σp = (σw[1:(end-1)] .+ σw[2:end]) / 2
		pf = ds["p_freq"][:]
		wa = ds["w_profile"][:] ./ pf'
		close(ds)

		axs[1].plot(
			pf./sum(pf)*100,σp,label=gID,legend="l",
			legend_kw=Dict("ncol"=>1,"frame"=>false)
		)
		
		axs[ii+1].plot(wa[:,40],σ,c=clr[11])
		axs[ii+1].plot(wa[:,45],σ,c=clr[10])
		axs[ii+1].plot(wa[:,50],σ,c=clr[9])
		axs[ii+1].plot(wa[:,55],σ,c=clr[8])
		axs[ii+1].plot(wa[:,60],σ,c=clr[7])
		axs[ii+1].plot(wa[:,65],σ,c=clr[6])
		axs[ii+1].plot(wa[:,70],σ,c=clr[5])
		axs[ii+1].plot(wa[:,75],σ,c=clr[4])
		axs[ii+1].plot(wa[:,80],σ,c=clr[3])
		axs[ii+1].format(
			ylim=(1,0.05),#yscale="log",
			ylabel=L"$\sigma$",xlabel=L"Pa s$^{-1}$",xlim=(0.025,-0.2),
			urtitle=gID
		)

		if ii == 6

			axs[ii+1].plot(
				wa[:,40],σ,label=L"\sigma=0.40",legend="r",c=clr[11],
				legend_kw=Dict("ncol"=>1,"frame"=>false)
			)
			axs[ii+1].plot(wa[:,45],σ,label=L"\sigma=0.45",legend="r",c=clr[10])
			axs[ii+1].plot(wa[:,50],σ,label=L"\sigma=0.50",legend="r",c=clr[9])
			axs[ii+1].plot(wa[:,55],σ,label=L"\sigma=0.55",legend="r",c=clr[8])
			axs[ii+1].plot(wa[:,60],σ,label=L"\sigma=0.60",legend="r",c=clr[7])
			axs[ii+1].plot(wa[:,65],σ,label=L"\sigma=0.65",legend="r",c=clr[6])
			axs[ii+1].plot(wa[:,70],σ,label=L"\sigma=0.70",legend="r",c=clr[5])
			axs[ii+1].plot(wa[:,75],σ,label=L"\sigma=0.75",legend="r",c=clr[4])
			axs[ii+1].plot(wa[:,80],σ,label=L"\sigma=0.80",legend="r",c=clr[3])

		end

		axs[1].format(ultitle="(a)",xlim=(0,7.5))
		axs[2].format(ultitle="(b)")
		axs[3].format(ultitle="(c)")
		axs[4].format(ultitle="(d)")
		axs[5].format(ultitle="(e)")
		axs[6].format(ultitle="(f)")
		axs[7].format(ultitle="(g)")

	end
	
	fig.savefig(plotsdir("02b-wvertprofile.png"),transparent=false,dpi=400)
	load(plotsdir("02b-wvertprofile.png"))
	
end

# ╔═╡ 53295e98-0833-4396-baee-58be12dfda35
md"
Here, we plot similar curves to Torri et al. (2017) for different regions.  However, there is one difference, in that we did not adjust for the C/E = 2 ratio like they did in their paper.  In this manner, we therefore see some differences in the vertical profiles here compared to Torri et al. (2017), although they used ERA-40, which is a relatively old dataset.  Nonetheless, we see that the peak of $\sigma_w$ is not as sharply centered around 0.5, which was what we would have expected from Torri et al. (2017).

Again, likely related to C/E = 2 ratio not being taken into account here, and also because we didn't use as long a dataset (Torri et al. (2017) used datasets from 1957 to 2002, much longer time period), so there will likely also be some noise.
"

# ╔═╡ 40d5ffc7-b59c-4715-ba65-525e49a36ccf
md"
### C. Sample Vertical Profiles
"

# ╔═╡ 7c663b1d-2043-4990-8b7f-cb7f757abbf9
load(plotsdir("02a-wwgt_pre-OTREC.png"))

# ╔═╡ f9552789-9b6a-4a36-b9e2-2bf5441fd963
md"
We wish to zoom in on the region centered about 276-278.5E and 1-3.5N to check the vertical profile of the region without the distraction of binning the profile by $\sigma$.  We define this region as \"OTREC_BTM\".
" 

# ╔═╡ 077fc9b0-4ab8-4604-a175-ac4b5a3dae9f
e5ds = ERA5Monthly(dtbeg=Date(2013,1,1),dtend=Date(2021,12,31),eroot=datadir())

# ╔═╡ c20ea9f0-27be-46ed-96aa-91d0c379b445
geo = RectRegion(
	"OTREC_BTM","OTREC","OTREC Bottom Heavy",
	[3.5,1,278.5,276],savegeo=false
)

# ╔═╡ a6f989cd-552d-41f0-b394-891dda857cad
md"Extract for GeoRegion? $(@bind doextract PlutoUI.Slider(0:1))"

# ╔═╡ a5aeebb7-c304-4eae-b590-a14ebc8befd0
begin
	if isone(doextract)
		evar_wp = SingleVariable("p_wwgt")
		extract(e5ds,evar_wp,ERA5Region(geo))
		md"Extracting W-weighted Mean Pressure data ..."
	else
		md"Not extracting W-weighted Mean Pressure data."
	end
end

# ╔═╡ bc4cd897-2d5f-4b5b-8306-e54a3b290393
begin
	if isone(doextract)
		p = era5Pressures(); p = p[p.>=10]
		for ip in p
		    evar_w = PressureVariable("w",hPa=ip)
		    extract(e5ds,evar_w,ERA5Region(geo))
		end
		md"Extracting vertical velocity data ..."
	else
		md"Not extracting vertical velocity data."
	end
end

# ╔═╡ 0e6c8ba8-7678-47d9-93ab-13a70ab21508
lsd = getLandSea(ERA5Region(geo),path=datadir("emask"))

# ╔═╡ 7d10f096-67e3-46b0-beb9-7307178036e3
begin
	
	nlon = length(lsd.lon)
	nlat = length(lsd.lat)
	wprf = zeros(32)
	
	var_wp = SingleVariable("p_wwgt")
	var_w  = Vector{PressureVariable}(undef,32)

	wp  = zeros(nlon,nlat)
	w   = zeros(nlon,nlat,32)
	cnt = 0
	
	pp  = era5Pressures(); pp = pp[pp.>=10]
	for ip = 1 : 32
		var_w[ip] = PressureVariable("w",hPa=pp[ip])
	end

	for dt = Date(2013) : Year(1) : Date(2021)

		ids = read(e5ds,var_wp,ERA5Region(geo),dt)
		wp .= nomissing(ids[var_wp.varID][:,:,1],NaN)
		close(ids)

		for ip = 1 : 32

			ids = read(e5ds,var_w[ip],ERA5Region(geo),dt)
			w[:,:,ip] .= ids[var_w[ip].varID][:,:,1]
			close(ids)

		end

		for ilon = 1 : nlon, ilat = 1 : nlat

			if !isnan(wp[ilon,ilat])
				global cnt += 1
				wprf[:] += w[ilon,ilat,:]
			end
			
		end
		
	end

	wprf = wprf / cnt
	md"Loading average vertical profile for the region ..."
	
end

# ╔═╡ c8c54bb6-fd87-465d-a25d-793d90d83e82
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=1/3,axwidth=1)
	
	a2[1].plot(wprf,pp)
	a2[1].format(
		ylim=(1000,10),xlim=(0.05,-0.15),
		xlabel=L"$\omega$ / Pa s$^{-1}$",
		ylabel="Pressure / hPa"
	)
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─30d2be4a-9bdd-4939-8785-413cd7a59f78
# ╠═1a69f197-a5a8-476b-86e7-c472061ca1d5
# ╟─d413d3cb-cf76-4aa5-80ff-c929d396929d
# ╟─92621f09-ed8c-46aa-afc0-4cab187dde68
# ╟─4227bb02-d323-4653-a071-279eff3626d9
# ╟─059b955b-1fa1-4a60-a136-730146cce921
# ╟─53295e98-0833-4396-baee-58be12dfda35
# ╟─40d5ffc7-b59c-4715-ba65-525e49a36ccf
# ╟─7c663b1d-2043-4990-8b7f-cb7f757abbf9
# ╟─f9552789-9b6a-4a36-b9e2-2bf5441fd963
# ╟─077fc9b0-4ab8-4604-a175-ac4b5a3dae9f
# ╟─c20ea9f0-27be-46ed-96aa-91d0c379b445
# ╟─a6f989cd-552d-41f0-b394-891dda857cad
# ╟─a5aeebb7-c304-4eae-b590-a14ebc8befd0
# ╟─bc4cd897-2d5f-4b5b-8306-e54a3b290393
# ╟─0e6c8ba8-7678-47d9-93ab-13a70ab21508
# ╟─7d10f096-67e3-46b0-beb9-7307178036e3
# ╟─c8c54bb6-fd87-465d-a25d-793d90d83e82
