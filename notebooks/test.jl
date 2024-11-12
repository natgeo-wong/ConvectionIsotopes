### A Pluto.jl notebook ###
# v0.19.46

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
	using GeoRegions
	using DelimitedFiles
	using NASAPrecipitation
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

# ╔═╡ 6ebfd76d-9ee3-451f-a63b-64cdf936d13b
begin
	dslnlt = NCDataset(datadir("wrf3","grid.nc"))
	lon = dslnlt["longitude"][:,:]
	lat = dslnlt["latitude"][:,:]
	close(dslnlt)
end

# ╔═╡ 7bebe8a4-c780-4464-8b14-0a87a63cef2f
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 1afc643b-4ec3-41e5-bb60-7fde65abf198
npd = IMERGMonthly(start=Date(2020),stop=Date(2020),path=datadir())

# ╔═╡ 7c42a107-16a7-49d7-949b-f11193cedeeb
lsd = getLandSea(npd,geo)

# ╔═╡ 3abd9e8a-5d8e-44ff-8a64-5b98217ddee1
begin
	npdds = read(npd,geo,Date(2020))
	prcp  = dropdims(mean(npdds["precipitation"][:,:,1:5],dims=3),dims=3) * 86400
	close(npdds)
end

# ╔═╡ 87b04b4c-5ccf-4332-99cc-40e5edfb5c02
begin
	dsH2O = NCDataset(datadir("wrf3","2D","RAINNC-daily.nc"))
	ds7dy = NCDataset(datadir("wrf3","2D","RAINNC-daily-smooth_07days.nc"))
	dsO18 = NCDataset(datadir("wrf3","2D","O18_RAINNC-daily-smooth_07days.nc"))

	prcp_H2O = dropdims(mean(dsH2O["RAINNC"][:,:,:],dims=3),dims=3)
	prcp_7dy = dropdims(mean(ds7dy["RAINNC"][:,:,4:(end-3)],dims=3),dims=3)
	prcp_O18 = dropdims(mean(dsO18["O18_RAINNC"][:,:,4:(end-3)],dims=3),dims=3)
	
	close(dsH2O)
	close(ds7dy)
	close(dsO18)
end

# ╔═╡ 431c7e32-551e-44b1-8c8d-7940a87487f7
begin
	pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=3)

	axs[1].pcolormesh(lsd.lon,lsd.lat,prcp',levels=2:20,cmap="Blues",extend="both")
	axs[2].pcolormesh(lon,lat,prcp_H2O,levels=2:20,cmap="Blues",extend="both")
	axs[3].pcolormesh(lon,lat,(prcp_O18./prcp_7dy .- 1)*1000,levels=-15:0,cmap="viridis",extend="both")

	for ax in axs
		ax.plot(x,y,c="k",lw=0.4)
		ax.format(xlim=(-100,-60),ylim=(-15,25))
	end
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╠═bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─6ebfd76d-9ee3-451f-a63b-64cdf936d13b
# ╠═7bebe8a4-c780-4464-8b14-0a87a63cef2f
# ╠═1afc643b-4ec3-41e5-bb60-7fde65abf198
# ╠═7c42a107-16a7-49d7-949b-f11193cedeeb
# ╟─3abd9e8a-5d8e-44ff-8a64-5b98217ddee1
# ╟─87b04b4c-5ccf-4332-99cc-40e5edfb5c02
# ╠═431c7e32-551e-44b1-8c8d-7940a87487f7
