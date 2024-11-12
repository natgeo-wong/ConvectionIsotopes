### A Pluto.jl notebook ###
# v0.19.36

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
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ a8431d46-46fe-11ec-2b8d-e39caffdabec
md"
# Figure 2. Station Rainfall and Isotopes

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ 5b392d8a-d68b-44df-ab57-d0e81aedd668
md"
### A. Loading the Processed Data
"

# ╔═╡ cceac73b-b7d9-4418-b145-d45d384b61b8
begin
	ds = NCDataset(datadir("processed.nc"))
	dtvec = ds["time"][:]
	lon   = ds["longitude"][:]
	lat   = ds["latitude"][:]
	prcp  = ds["prcp"][:] 
	prcpg = ds["prcpg"][:] 
	prcps = ds["prcps"][:] 
	δ18Oμ = ds["δ18Oμ"][:] 
	δ18Oσ = ds["δ18Oσ"][:] 
	δ2Hμ  = ds["δ2Hμ"][:] 
	δ2Hσ  = ds["δ2Hσ"][:] 
	close(ds)
	md"Loading processed station and GPM rainfall and isotopic data ..."
end

# ╔═╡ 0374f823-90a0-4dc6-be0e-9449a7fa7202
begin
	arr = [[1,1,1,2,3],[4,4,4,5,6],[7,7,7,8,9]]
	pplt.close(); fig,axs = pplt.subplots(arr,aspect=3,axwidth=3.5,sharex=0,sharey=0,wspace=[0,0,1,8],hspace=[1,1])

	for istn in 1 
		axs[1].scatter(dtvec,prcpg[:,istn],s=1,c="b")
		axs[1].scatter(dtvec,prcps[:,istn],s=1,c="r")
		axs[2].scatter(prcps[:,istn],prcpg[:,istn],s=1,c="k")
		axs[3].scatter(δ18Oμ[:,istn],δ2Hμ[:,istn],s=1,c="k")
	end

	for istn in 2 : 4
		axs[4].scatter(dtvec,prcpg[:,istn],s=1,c="b")
		axs[4].scatter(dtvec,prcps[:,istn],s=1,c="r")
		axs[5].scatter(prcps[:,istn],prcpg[:,istn],s=1,c="k")
		axs[6].scatter(δ18Oμ[:,istn],δ2Hμ[:,istn],s=1,c="k")
	end

	for istn in 5 : 12
		axs[7].scatter(dtvec,prcpg[:,istn],s=1,c="b")
		axs[7].scatter(dtvec,prcps[:,istn],s=1,c="r")
		axs[8].scatter(prcps[:,istn],prcpg[:,istn],s=1,c="k")
		axs[9].scatter(δ18Oμ[:,istn],δ2Hμ[:,istn],s=1,c="k")
	end

	for ii in [1,4,7]
		axs[ii].format(ylim=(-10,210),xlim=(Date(2019,1),Date(2021,12)))
	end

	for ii in [1,2,3,4,5,6]
		axs[ii].format(xticklabels=[])
	end

	for ii in [2,5,8]
		axs[ii].format(xlim=(-10,210),ylim=(-10,210),ytickloc="r")
	end

	for ii in [3,6,9]
		axs[ii].format(xlim=(-20,2),ylim=(-150,25),ytickloc="r")
	end

	for ax in axs
		ax.format(xminorlocator=[])
	end

	axs[1].format(ltitle="(a) Time-Series", ultitle="(i) San Andres",)
	axs[2].format(ltitle="(b) Station vs GPM")
	axs[3].format(ltitle=L"(c) $\delta^2$H vs $\delta^{18}$O")
	axs[4].format(ultitle="(ii) Colombia Daily Stations",ylabel=L"Rain Rate / mm day$^{-1}$")
	axs[5].format(ylabel=L"GPM Rainfall / mm day$^{-1}$")
	axs[6].format(ylabel=L"$\delta^2$H")
	axs[7].format(ultitle="(iii) Costa Rica Stations")
	axs[8].format(xlabel=L"Station Rainfall / mm day$^{-1}$")
	axs[9].format(xlabel=L"$\delta^{18}$O")
		
	
	fig.savefig(plotsdir("Figure2.png"),transparent=false,dpi=400)
	load(plotsdir("Figure2.png"))
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─5b392d8a-d68b-44df-ab57-d0e81aedd668
# ╟─cceac73b-b7d9-4418-b145-d45d384b61b8
# ╟─0374f823-90a0-4dc6-be0e-9449a7fa7202
