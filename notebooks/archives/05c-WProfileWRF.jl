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
# 05c. Testing for Isotope vs W-Weighted Pressure
"

# ╔═╡ dd06e4a3-b507-446a-ae9c-fe3e173c434e
begin
	cst  = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	clon = cst[:,1]
	clat = cst[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ cc24d989-0b3e-4213-bce0-efa420191ac8
imo = 8

# ╔═╡ c853d30e-5a51-47f4-a6be-88ab857253c1
md"
### A. Loading some Sample Data
"

# ╔═╡ 14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
begin
	ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-RAINNC"))
	lon  = ds["longitude"][:] .+ 360
	lat  = ds["latitude"][:]
	rain = ds["RAINNC"][:]
	close(ds)
	ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-HDO_RAINNC"))
	rHDO = ds["HDO_RAINNC"][:]
	rHDO = (rHDO./rain .-1) *1000
	close(ds)
	ds   = NCDataset(datadir("wrf","2D","$(uppercase(monthabbr(imo)))-p_wwgt"))
	pwgt = ds["p_wwgt"][:]
	close(ds)
	ds = NCDataset(datadir("wrf","3D","$(uppercase(monthabbr(imo)))-W"))
	p  = ds["P"][:]
	z  = ds["Z"][:]
	w  = ds["W"][:]
	close(ds)
end

# ╔═╡ 9f0922fe-c451-427c-8fef-b85cfa4229f3
begin
	pplt.close(); f1,a1 = pplt.subplots(axwidth=2.5)
	
	c1 = a1[1].pcolormesh(lon,lat,rHDO,extend="both",levels=-150:10:0)
	a1[1].plot(clon,clat,lw=0.5,c="k")
	a1[1].format(xlim=(270,285),ylim=(0,15))

	f1.colorbar(c1,length=0.75)
	f1.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 4c223141-fb32-4491-af88-026b0411ddde
begin
	pplt.close(); f2,a2 = pplt.subplots(axwidth=2.5)
	
	c2 = a2[1].pcolormesh(lon,lat,p[:,:,1])
	a2[1].plot(clon,clat,lw=0.5,c="k")
	a2[1].format(xlim=(270,285),ylim=(0,15))

	f2.colorbar(c2,length=0.75)
	f2.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 247a88a8-596f-4b6a-9c5c-aae91365228a
begin
	pplt.close(); f3,a3 = pplt.subplots(axwidth=2.5)
	
	c3 = a3[1].pcolormesh(lon,lat,pwgt./100,extend="both")
	a3[1].plot(clon,clat,lw=0.5,c="k")
	a3[1].format(xlim=(270,285),ylim=(0,15))

	f3.colorbar(c3,length=0.75)
	f3.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ e884be27-b9ec-45ba-b800-afccbfb54a78
begin
	pplt.close(); f4,a4 = pplt.subplots(axwidth=2.5)

	rainl = rain[:]
	pwgtl = pwgt[:] ./ 100
	rHDOl = rHDO[:]

	rainb = 0 : 2 : 58; nr = length(rainb)
	rainu = 2 : 2 : 60
	rainm = 0 : 2 : 60
	pwgtb = 100 : 25 : 925; np = length(pwgtb)
	pwgtu = 125 : 25 : 950
	pwgtm = 100 : 25 : 950

	binmat = zeros(nr,np)

	for ip = 1 : np, ir = 1 : nr

		binmat[ir,ip] = mean(rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .&
									(rainl.>rainb[ir]) .& (rainl.<rainu[ir])])

	end
	
	c4 = a4[1].pcolormesh(rainm,pwgtm,binmat',cmap="viridis",)#levels=-65:-50)
	a4[1].format(ylim=(maximum(pwgtm),minimum(pwgtm)))

	f4.colorbar(c4,length=0.75)
	f4.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
# ╟─3c4b037f-c854-4148-ae6e-fff50bbb51cf
# ╟─4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
# ╟─dd06e4a3-b507-446a-ae9c-fe3e173c434e
# ╠═cc24d989-0b3e-4213-bce0-efa420191ac8
# ╟─c853d30e-5a51-47f4-a6be-88ab857253c1
# ╟─14ec25a0-9bd8-4e61-a7c1-4eaee00f8f4e
# ╟─9f0922fe-c451-427c-8fef-b85cfa4229f3
# ╟─4c223141-fb32-4491-af88-026b0411ddde
# ╟─247a88a8-596f-4b6a-9c5c-aae91365228a
# ╟─e884be27-b9ec-45ba-b800-afccbfb54a78
