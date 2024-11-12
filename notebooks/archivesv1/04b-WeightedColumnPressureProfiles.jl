### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 04b. Vertical Wind Profiles for Weighted Column-Mean $\sigma$

Here, we plot the profiles of vertical velocity for different $\sigma$ in different regions within the defined OTREC region.
"

# ╔═╡ 30d2be4a-9bdd-4939-8785-413cd7a59f78
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 1a69f197-a5a8-476b-86e7-c472061ca1d5
geolist = ["OTREC","CIS_SAN","CIS_CLB","CIS_MNT","CIS_INL","CIS_PAC","CIS_CRB"]

# ╔═╡ d413d3cb-cf76-4aa5-80ff-c929d396929d
load(plotsdir("01a-ExploreDomain.png"))

# ╔═╡ 92621f09-ed8c-46aa-afc0-4cab187dde68
clr = pplt.Colors("oslo",12)

# ╔═╡ 4227bb02-d323-4653-a071-279eff3626d9
md"### B. Averaged vertical velocity profiles for each $\sigma_w$ bin"

# ╔═╡ 059b955b-1fa1-4a60-a136-730146cce921
begin
	
	pplt.close()
	arr = [1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8]
	wsp = [0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0,0.1,0,0,0]
	fig,axs = pplt.subplots(arr,aspect=1/6,axwidth=0.5,sharex=0,wspace=wsp)

	for ii = 1 : 7

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
		
		axs[ii+1].plot(wa[:,40],σ,c=clr[10])
		axs[ii+1].plot(wa[:,45],σ,c=clr[9])
		axs[ii+1].plot(wa[:,50],σ,c=clr[8])
		axs[ii+1].plot(wa[:,55],σ,c=clr[7])
		axs[ii+1].plot(wa[:,60],σ,c=clr[6])
		axs[ii+1].plot(wa[:,65],σ,c=clr[5])
		axs[ii+1].plot(wa[:,70],σ,c=clr[4])
		axs[ii+1].format(
			ylim=(1,0.05),#yscale="log",
			ylabel=L"$\sigma$",xlabel=L"Pa s$^{-1}$",xlim=(0.1,-0.5),
			urtitle=gID
		)

		if ii == 7

			axs[ii+1].plot(
				wa[:,40],σ,label=L"\sigma=0.40",legend="r",c=clr[10],
				legend_kw=Dict("ncol"=>1,"frame"=>false)
			)
			axs[ii+1].plot(wa[:,45],σ,label=L"\sigma=0.45",legend="r",c=clr[9])
			axs[ii+1].plot(wa[:,50],σ,label=L"\sigma=0.50",legend="r",c=clr[8])
			axs[ii+1].plot(wa[:,55],σ,label=L"\sigma=0.55",legend="r",c=clr[7])
			axs[ii+1].plot(wa[:,60],σ,label=L"\sigma=0.60",legend="r",c=clr[6])
			axs[ii+1].plot(wa[:,65],σ,label=L"\sigma=0.65",legend="r",c=clr[5])
			axs[ii+1].plot(wa[:,70],σ,label=L"\sigma=0.70",legend="r",c=clr[4])
			axs[ii+1].format(
				ylim=(1,0.05),xlim=(0.05,-0.25),#yscale="log",
				ylabel=L"$\sigma$",xlabel=L"Pa s$^{-1}$"
			)

		end

		axs[1].format(ultitle="(a)")
		axs[2].format(ultitle="(b)")
		axs[3].format(ultitle="(c)")
		axs[4].format(ultitle="(d)")
		axs[5].format(ultitle="(e)")
		axs[6].format(ultitle="(f)")
		axs[7].format(ultitle="(g)",xlim=(0.1,-0.5))
		axs[8].format(ultitle="(h)")

	end
	
	fig.savefig(plotsdir("04a-wvertprofile.png"),transparent=false,dpi=150)
	load(plotsdir("04a-wvertprofile.png"))
	
end

# ╔═╡ 53295e98-0833-4396-baee-58be12dfda35
md"
Here, we plot similar curves to Torri et al. (2017) for different regions.  However, there is one difference, in that we did not adjust for the C/E = 2 ratio like they did in their paper.  In this manner, we therefore see some differences in the vertical profiles here compared to Torri et al. (2017), although they used ERA-40, which is a relatively old dataset.  Nonetheless, we see that the peak of $\sigma_w$ is not as sharply centered around 0.5, which was what we would have expected from Torri et al. (2017).

Again, likely related to C/E = 2 ratio not being taken into account here, and also because we didn't use as long a dataset (Torri et al. (2017) used datasets from 1957 to 2002, much longer time period), so there will likely also be some noise.
"

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
