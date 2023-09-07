### A Pluto.jl notebook ###
# v0.19.0

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

# ╔═╡ 8488ffe7-7781-44aa-8bb0-742cf815f5cd
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b5301b23-1ef8-4cca-b92f-c16773c510a0
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using NCDatasets
	using NumericalIntegration
	using PlutoUI
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColumbiaIsotope project..."
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

# ╔═╡ 0aba4321-ae74-4c3b-859d-646dfa346480
md"
### A. Just some OLR Animation Testing ...
"

# ╔═╡ 6d7d8839-11f5-4734-bbb6-1ba04175c569
md"Create Animation? $(@bind createanim PlutoUI.Slider(0:1))"

# ╔═╡ d4a3296a-2eb9-462f-a555-2f558cb46101
begin
	if isone(createanim)
	tt = 0
		for it = 1 : 10, ihr = 1 : 8
	
			dt = Date(2019,8,1) + Day(it-1)
			fol_2D = datadir("wrf","2D")
			fnc_2D = joinpath(fol_2D,"wrfout2D_d02_$(dt)_00:00:00")
			d2D = NCDataset(fnc_2D)
			lon = d2D["XLONG"][:,1,1]
			lat = d2D["XLAT"][1,:,1]
			var = d2D["OLR"][:]
			close(d2D)
			
			global tt += 1
			pplt.close(); f1,a1 = pplt.subplots(axwidth=2)
			
			c1 = a1[1].pcolormesh(
				lon,lat,var[:,:,ihr]',
				extend="both",levels=100:20:300
			)
			a1[1].plot(x,y,lw=1,c="k")
			a1[1].colorbar(c1)
			a1[1].format(xlim=(-90,-75),ylim=(0,15))
		
			f1.savefig(
				joinpath("test$(@sprintf("%03d",tt)).png"),
				transparent=false,dpi=200
			)
		end
		md"Creating sample animation from OLR data ..."
	else
		md"Not creating animation, skipping to next cell ..."
	end

end

# ╔═╡ 0eb4efc5-a7cf-4b08-83a9-44e568752cc3
md"
### B. Precipitation Accumulation?
"

# ╔═╡ 503c98f2-eb53-41a2-a46a-611a92fd2e79
begin
	ds  = NCDataset(datadir("wrf","AUGp02","wrfout_d02_2019-08-11_00:00:00"))
	tln = ds["XLONG"][:,1,1]
	tlt = ds["XLAT"][1,:,1]
	rn1 = ds["RAINNC"][:]
	d2H1 = ds["HDO_RAINNC"][:]
	close(ds)
	ds  = NCDataset(datadir("wrf","AUGp02","wrfout_d02_2019-08-21_00:00:00"))
	rn2 = ds["RAINNC"][:]
	d2H2 = ds["HDO_RAINNC"][:]
	close(ds)
end

# ╔═╡ 93517ec7-e38e-4b66-8fd9-eacb7946dd68
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,axwidth=2)
				
	c2_1 = a2[1].pcolormesh(
		tln,tlt,(rn2[:,:,1].-rn1[:,:,1])'./9,
		cmap="Blues",extend="both",levels=vcat(5:5:50)
	)
	a2[1].plot(x,y,lw=1,c="k")
	a2[1].colorbar(c2_1,label=L"Rain / mm day$^{-1}$")
	a2[1].format(xlim=(-90,-75),xlocator=-90:3:-75,ylim=(0,15),ylocator=0:3:15)
	
	c2_2 = a2[2].pcolormesh(
		tln,tlt,((d2H2[:,:,1].-d2H1[:,:,1])./(rn2[:,:,1].-rn1[:,:,1]) .- 1)' *1000,
		cmap="Blues",extend="both"
	)
	a2[2].plot(x,y,lw=1,c="k")
	a2[2].colorbar(c2_2,label=L"$\delta^2$H / $\perthousand$")
	a2[2].format(xlim=(-90,-75),xlocator=-90:3:-75,ylim=(0,15),ylocator=0:3:15)
	
	f2.savefig("testaccum.png",transparent=false,dpi=200)
	load("testaccum.png")
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
	d3D = NCDataset(datadir("wrf","AUGp02","wrfout3D_d02_2019-08-11_00:00:00"))
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

# ╔═╡ e5057cbb-5658-419a-a65c-7fe9c74d3aa5
integrate(lvl[2:(end-1)],lvl[2:(end-1)].*w) / integrate(lvl[2:(end-1)],w)

# ╔═╡ Cell order:
# ╟─a6a3baf2-c81a-11ec-0d65-edc25a50a22d
# ╟─8488ffe7-7781-44aa-8bb0-742cf815f5cd
# ╟─b5301b23-1ef8-4cca-b92f-c16773c510a0
# ╟─05c78719-5347-45c8-b26d-af98fb56a30a
# ╟─0aba4321-ae74-4c3b-859d-646dfa346480
# ╟─6d7d8839-11f5-4734-bbb6-1ba04175c569
# ╟─d4a3296a-2eb9-462f-a555-2f558cb46101
# ╟─0eb4efc5-a7cf-4b08-83a9-44e568752cc3
# ╟─503c98f2-eb53-41a2-a46a-611a92fd2e79
# ╟─93517ec7-e38e-4b66-8fd9-eacb7946dd68
# ╟─30a9c234-089d-4915-b588-15365591195d
# ╟─ee2f7b09-5be0-4c56-9ee2-ee15e8e35759
# ╠═a86e485a-f2df-4bf8-8130-35819888d5d9
# ╠═a3eb0f0f-48ac-4492-80eb-121d53f8b19b
# ╠═40e84bd1-16c7-4283-927a-eaad0b33290d
# ╟─584e8235-69b0-46d3-b5fc-05b1a06e7bbf
# ╠═e5057cbb-5658-419a-a65c-7fe9c74d3aa5
