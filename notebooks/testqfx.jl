### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ e32a00ee-5f32-47a1-a983-91fb77bc5d18
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
begin
	@quickactivate "ConvectionIsotopes"
	using Dates
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 03b. Station W-Weighted Pressure, $\sigma$
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 59c930cd-5b7f-4047-8660-615148d1bd9f
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ a926c18c-9536-4be9-94dc-02aa9104f243
begin
	dsg = NCDataset(datadir("wrf","processed","OTREC_STN02-QBUDGET-daily.nc"));
	dt1 = dsg["time"][:]
	∇   = dsg["∇"][:] * 20000
	close(dsg)
	dsg = NCDataset(datadir("wrf","processed","OTREC_STN02-∇.nc"));
	dt2 = dsg["time"][:]
	∇c  = dsg["∇"][:] * 20000
	close(dsg)
end

# ╔═╡ fde75120-f6f0-49e7-82c5-c92524fa8841
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=4,axwidth=4)
	
	axs[1].plot(dt1,∇)
	axs[1].plot(dt2,∇c)

	# axs[1].format(ylim=(-50,50))
	
	fig.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 80d6921a-73e7-4d84-987e-6e08641b7935
begin
	pplt.close(); f2,a2 = pplt.subplots()
	
	a2[1].scatter(∇,∇c,s=10)
	a2[1].plot([-75,30],[-75,30],c="k",linestyle="--")
	# a2[1].plot([-45,15]/2,[-45,15],c="k",linestyle="--")
	a2[1].format(xlim=(-75,30),ylim=(-75,30))
	
	f2.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╠═a926c18c-9536-4be9-94dc-02aa9104f243
# ╠═fde75120-f6f0-49e7-82c5-c92524fa8841
# ╠═80d6921a-73e7-4d84-987e-6e08641b7935
