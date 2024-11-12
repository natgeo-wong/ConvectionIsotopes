### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ e4b60273-bf9b-4d79-b460-97a942a333f8
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ a57a1f6a-6a7b-44ed-baaf-fe9980adb72e
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using Glob
	using NCDatasets
	using NumericalIntegration
	using PlutoUI
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ 497776f8-c837-11ec-2301-5170cf3d601f
md"
# 05b. Monthly WRF Aggregates and Means

In this notebook, we compile and monthly aggregates and means for the following variables:
* Vertical wind (Pa s$^{-1}$)
* Rain accumlation / daily rate, $\delta^{18}$O and $\delta^{2}$H fraction
"

# ╔═╡ f4d62c89-060e-44f4-9def-573d733c88ac
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 30b8b534-4dc1-4038-8ce9-b3fa01d2dfa7
imo = 12

# ╔═╡ f67ffe0e-7cb0-4c26-925c-244106fd58e1
md"
### A. Extracting and Saving Monthly 2D Data
"

# ╔═╡ 3e06da77-0fe0-4422-923b-f71738d1731f
@bind ncID Select([
	"RAINNC" => "Total Accumulated Rainfall",
	"HDO_RAINNC" => "Total Accumulated HDO Rainfall",
	"O18_RAINNC" => "Total Accumulated 18O Rainfall",
])

# ╔═╡ 670c4316-1f62-46f8-8c21-ae71d57cd39d
begin
	gds  = NCDataset(datadir("wrf","AUGp02","wrfout_d02_2019-08-11_00:00:00"))
	tln = gds["XLONG"][:,1,1]
	tlt = gds["XLAT"][1,:,1]
	lon = gds["XLONG"][:,:,1]
	lat = gds["XLAT"][:,:,1]
	close(gds)
end

# ╔═╡ 2532b35e-03bf-403c-b005-e0fdbf42ef74
begin
	v1  = zeros(Float32,543,543)
	v2  = zeros(Float32,543,543)
	md"Preallocation of arrays ..."
end

# ╔═╡ e077e997-aed6-4c7a-a1d2-97d754517dd8
begin
	dvar = zeros(Float32,543,543)
	for ifol = 1 : 3
		dsvec = glob(
			"wrfout_d02*",
			datadir("wrf","$(uppercase(monthabbr(imo)))p0$ifol")
		)
		fnc1 = dsvec[2]
		fnc2 = dsvec[end]

		ds = NCDataset(fnc1)
		NCDatasets.load!(ds[ncID].var,v1,:,:,1)
		close(ds)

		ds = NCDataset(fnc2)
		NCDatasets.load!(ds[ncID].var,v2,:,:,1)
		close(ds)

		for iy = 1 : 543, ix = 1 : 543
			dvar[ix,iy] += (v2[ix,iy] - v1[ix,iy])
		end
	end
	dvar = dvar / daysinmonth(Date(2019,imo))
	md"Loading monthly accumulations"
end

# ╔═╡ 224f19e7-e62c-47b2-bf16-c281764b5ba2
begin
	pplt.close(); f1,a1 = pplt.subplots(axwidth=2)
	
	c = a1[1].contourf(tln,tlt,dvar',cmap="blues",extend="both",levels=5:2.5:35)
	a1[1].plot(x,y,lw=1,c="k")
	a1[1].format(xlim=(-90,-75),xlocator=-90:3:-75,ylim=(0,15),ylocator=0:3:15)
	a1[1].colorbar(c)

	f1.savefig("testmean.png",transparent=false,dpi=150)
	load("testmean.png")
end

# ╔═╡ 4a1ae392-d9b5-4232-8e99-42c395c2364b
begin
	fnc = datadir("wrf","2D","$(uppercase(monthabbr(imo)))-$ncID")
	if isfile(fnc); rm(fnc,force=true) end
	ds  = NCDataset(fnc,"c")
	defDim(ds,"longitude",length(tln))
	defDim(ds,"latitude",length(tlt))

	nclon = defVar(ds,"longitude",Float32,("longitude","latitude",))
	nclat = defVar(ds,"latitude",Float32,("longitude","latitude",))
	ncvar = defVar(ds,ncID,Float32,("longitude","latitude",))

	nclon[:] = lon
	nclat[:] = lat
	ncvar[:] = dvar
	
	close(ds)
end

# ╔═╡ 90ad821d-5864-4e65-b1de-392c794677eb
md"
### B. Extracting and Saving Monthly 3D Data
"

# ╔═╡ f9b1e976-788d-4f09-883d-979658f419f0
begin
	v3D_1 = zeros(Float32,543,543,50,8)
	v3D_2 = zeros(Float32,543,543,51,8)
	md"Temporary Arrays ..."
end

# ╔═╡ d1c38d2c-49c8-4fe8-a893-9b5704b94a22
begin
	v3D_p = zeros(Float32,543,543,50)
	v3D_z = zeros(Float32,543,543,51)
	v3D_w = zeros(Float32,543,543,51)
	for ifol = 1 : 3
		dsvec = glob(
			"wrfout3D_*",
			datadir("wrf","$(uppercase(monthabbr(imo)))p0$ifol")
		)
		for inc = 2 : (length(dsvec) - 1)
			fnc = dsvec[inc]
			d3D = NCDataset(fnc)
			
			NCDatasets.load!(d3D["P"].var,v3D_1,:,:,:,:)
			for it = 1 : 8
				v3D_p[:,:,:] += view(v3D_1,:,:,:,it)
			end
			
			NCDatasets.load!(d3D["PB"].var,v3D_1,:,:,:,:)
			for it = 1 : 8
				v3D_p[:,:,:] += view(v3D_1,:,:,:,it)
			end
			
			NCDatasets.load!(d3D["PH"].var,v3D_2,:,:,:,:)
			for it = 1 : 8
				v3D_z[:,:,:] += view(v3D_2,:,:,:,it)
			end
			
			NCDatasets.load!(d3D["PHB"].var,v3D_2,:,:,:,:)
			for it = 1 : 8
				v3D_z[:,:,:] += view(v3D_2,:,:,:,it)
			end
			
			NCDatasets.load!(d3D["W"].var,v3D_2,:,:,:,:)
			for it = 1 : 8
				v3D_w[:,:,:] += view(v3D_2,:,:,:,it)
			end
			
			close(d3D)
		end
	end
	v3D_p = v3D_p / daysinmonth(Date(2019,imo)) / 8
	v3D_z = v3D_z / daysinmonth(Date(2019,imo)) / 8
	v3D_w = v3D_w / daysinmonth(Date(2019,imo)) / 8
	md"Loading monthly accumulations"
end

# ╔═╡ 4f31465f-7fb4-4ae7-bf52-3337852652e5
begin
	f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-W")
	if isfile(f3D); rm(f3D,force=true) end
	d3D = NCDataset(f3D,"c")
	
	defDim(d3D,"longitude",length(tln))
	defDim(d3D,"latitude",length(tlt))
	defDim(d3D,"level",50)
	defDim(d3D,"level_staggered",51)

	nclon_3D = defVar(d3D,"longitude",Float32,("longitude","latitude",))
	nclat_3D = defVar(d3D,"latitude",Float32,("longitude","latitude",))
	ncp_3D = defVar(d3D,"P",Float32,("longitude","latitude","level",))
	ncz_3D = defVar(d3D,"Z",Float32,("longitude","latitude","level_staggered",))
	ncw_3D = defVar(d3D,"W",Float32,("longitude","latitude","level_staggered",))

	nclon_3D[:] = lon
	nclat_3D[:] = lat
	ncp_3D[:] = v3D_p
	ncz_3D[:] = v3D_z
	ncw_3D[:] = v3D_w
	
	close(d3D)
end

# ╔═╡ Cell order:
# ╟─497776f8-c837-11ec-2301-5170cf3d601f
# ╟─e4b60273-bf9b-4d79-b460-97a942a333f8
# ╟─a57a1f6a-6a7b-44ed-baaf-fe9980adb72e
# ╟─f4d62c89-060e-44f4-9def-573d733c88ac
# ╠═30b8b534-4dc1-4038-8ce9-b3fa01d2dfa7
# ╟─f67ffe0e-7cb0-4c26-925c-244106fd58e1
# ╟─3e06da77-0fe0-4422-923b-f71738d1731f
# ╟─670c4316-1f62-46f8-8c21-ae71d57cd39d
# ╟─2532b35e-03bf-403c-b005-e0fdbf42ef74
# ╟─e077e997-aed6-4c7a-a1d2-97d754517dd8
# ╟─224f19e7-e62c-47b2-bf16-c281764b5ba2
# ╟─4a1ae392-d9b5-4232-8e99-42c395c2364b
# ╟─90ad821d-5864-4e65-b1de-392c794677eb
# ╟─f9b1e976-788d-4f09-883d-979658f419f0
# ╟─d1c38d2c-49c8-4fe8-a893-9b5704b94a22
# ╟─4f31465f-7fb4-4ae7-bf52-3337852652e5
