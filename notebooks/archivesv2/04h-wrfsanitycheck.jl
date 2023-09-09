### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 8488ffe7-7781-44aa-8bb0-742cf815f5cd
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b5301b23-1ef8-4cca-b92f-c16773c510a0
begin
	@quickactivate "ConvectionIsotopes"
	using Dates
	using DelimitedFiles
	using NCDatasets
	using PlutoUI
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ a6a3baf2-c81a-11ec-0d65-edc25a50a22d
md"
# 05a. Exploring WRF Data Output
"

# ╔═╡ 05c78719-5347-45c8-b26d-af98fb56a30a
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 0eb4efc5-a7cf-4b08-83a9-44e568752cc3
md"
### Sanity Check at Pressure Levels
"

# ╔═╡ 503c98f2-eb53-41a2-a46a-611a92fd2e79
begin
	ds  = NCDataset(datadir("wrf","3D","NOV-QVAPOR_ISO"))
	tln = ds["longitude"][:]
	tlt = ds["latitude"][:]
	qvp = ds["QVAPOR"][:]
	hdo = ds["HDO_QVAPOR"][:]
	O18 = ds["O18_QVAPOR"][:]
	δ2H = (hdo ./ qvp .- 1) * 1000
	δ18 = (O18 ./ qvp .- 1) * 1000
	close(ds)
	ds  = NCDataset(datadir("wrf","3D","NOV-W"))
	pre = ds["P"][:]
	close(ds)
end

# ╔═╡ 41f67bc0-df2e-4623-8c37-ee9800d8f7ac
pre_lvl = 300

# ╔═╡ c8942f27-46fd-4682-8895-64deee292786
begin
	ip = argmin(abs.(pre./100 .-pre_lvl),dims=3)
	md"Finding the right pressure levels ..."
end

# ╔═╡ 93517ec7-e38e-4b66-8fd9-eacb7946dd68
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,axwidth=2)
				
	c2_1 = a2[1].pcolormesh(
		tln,tlt,δ2H[ip][:,:,1],
		cmap="Blues",extend="both",#levels=vcat(5:5:50)
	)
	a2[1].plot(x,y,lw=1,c="k")
	a2[1].colorbar(c2_1,label=L"$\delta^{2}$H / $\perthousand$",loc="l")
	a2[1].format(xlim=(-90,-75),xlocator=-90:3:-75,ylim=(0,15),ylocator=0:3:15)
	
	c2_2 = a2[2].pcolormesh(
		tln,tlt,δ18[ip][:,:,1],
		cmap="Blues",extend="both",#levels=vcat(5:5:50)
	)
	a2[2].plot(x,y,lw=1,c="k")
	a2[2].colorbar(c2_2,label=L"$\delta^{18}$O / $\perthousand$")
	a2[2].format(xlim=(-90,-75),xlocator=-90:3:-75,ylim=(0,15),ylocator=0:3:15)

	for ax in a2
		ax.format(suptitle="Pressure Level: $pre_lvl hPa")
	end
	
	f2.savefig(plotsdir("04h-sanitycheck-$pre_lvl.png"),transparent=false,dpi=200)
	load(plotsdir("04h-sanitycheck-$pre_lvl.png"))
end

# ╔═╡ 30a9c234-089d-4915-b588-15365591195d
md"
Okay, so the accumulation starts from the beginning of each run and ends at the end of each run.  So we need to just sum all the accumulation for each month, and then find the delD and delO18 for each month, and their isotopes, etc.
"

# ╔═╡ ee2f7b09-5be0-4c56-9ee2-ee15e8e35759
md"
### C. Loading some 3D data and vertical profiles?
"

# ╔═╡ a86e485a-f2df-4bf8-8130-35819888d5d9
ilon = 50

# ╔═╡ a3eb0f0f-48ac-4492-80eb-121d53f8b19b
ilat = 100

# ╔═╡ 40e84bd1-16c7-4283-927a-eaad0b33290d
begin
	d3D = NCDataset(datadir("wrf","wrfout","AUGp02","wrfout3D_d02_2019-08-11_00:00:00"))
	lvl = dropdims(mean(d3D["PB"][ilon,ilat,:,:],dims=2),dims=2) / 100
	lvl = dropdims(mean(d3D["P"][ilon,ilat,:,:],dims=2),dims=2)  / 100 .+ lvl
	geo = dropdims(mean(d3D["PHB"][ilon,ilat,:,:],dims=2),dims=2) / 100
	geo = dropdims(mean(d3D["PH"][ilon,ilat,:,:],dims=2),dims=2)  / 100 .+ geo
	geo = (geo[1:(end-1)].+geo[2:end]) / 2
	w   = dropdims(mean(d3D["W"][ilon,ilat,:,:],dims=2),dims=2)
	w   = (w[1:(end-1)].+w[2:end]) / 2
	w   = w[2:(end-1)] .* (lvl[3:end].-lvl[1:(end-2)]) ./ (geo[3:end].-geo[1:(end-2)]) * 100
	close(d3D)
end

# ╔═╡ 584e8235-69b0-46d3-b5fc-05b1a06e7bbf
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=1/3,axwidth=1)
	
	a3[1].plot(w,lvl[2:(end-1)])
	a3[1].format(ylim=(1000,100),xlim=(2.5,-2.5))
	
	f3.savefig("test3D.png",transparent=false,dpi=150)
	load("test3D.png")
end

# ╔═╡ Cell order:
# ╟─a6a3baf2-c81a-11ec-0d65-edc25a50a22d
# ╟─8488ffe7-7781-44aa-8bb0-742cf815f5cd
# ╟─b5301b23-1ef8-4cca-b92f-c16773c510a0
# ╟─05c78719-5347-45c8-b26d-af98fb56a30a
# ╟─0eb4efc5-a7cf-4b08-83a9-44e568752cc3
# ╟─503c98f2-eb53-41a2-a46a-611a92fd2e79
# ╠═41f67bc0-df2e-4623-8c37-ee9800d8f7ac
# ╟─c8942f27-46fd-4682-8895-64deee292786
# ╠═93517ec7-e38e-4b66-8fd9-eacb7946dd68
# ╟─30a9c234-089d-4915-b588-15365591195d
# ╟─ee2f7b09-5be0-4c56-9ee2-ee15e8e35759
# ╠═a86e485a-f2df-4bf8-8130-35819888d5d9
# ╠═a3eb0f0f-48ac-4492-80eb-121d53f8b19b
# ╠═40e84bd1-16c7-4283-927a-eaad0b33290d
# ╠═584e8235-69b0-46d3-b5fc-05b1a06e7bbf
