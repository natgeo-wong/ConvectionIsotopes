### A Pluto.jl notebook ###
# v0.20.4

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
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 04a. Moisture Budget from WRF
"

# ╔═╡ c7f79feb-b3f8-4154-a2f3-0c15d48bf61d
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ a15013a7-d990-40db-b608-a18035dcb437
begin
	infody  = stninfody()
	md"Loading station location information ..."
end

# ╔═╡ d447fdfd-b5c5-44f7-a2a2-161d5bb157b5
geo = GeoRegion("OTREC_wrf_stn02",path=srcdir())

# ╔═╡ 43bb0b2c-995a-4829-a30e-b953651305de
begin
	ds = NCDataset(datadir("wrf3","regridded","era-RAINNC-20190801_20201231.nc"))
	dt  = ds["time"][:]; nt = length(dt)
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	prc = ds["RAINNC"][:,:,:]
	close(ds)
end

# ╔═╡ 6bd91e94-ac11-4867-bbbb-b534707fce64
nt

# ╔═╡ d7eec81d-2bf8-4eb8-808d-d6fbf0fd3fd4
doanim = true

# ╔═╡ 92fd5032-baf7-41f3-bdb8-d0e273171923
it = 1006

# ╔═╡ 839bf2a5-0f04-45d1-b683-99f275619ff5
begin
	if doanim
		for it = 1 : 1000
			pplt.close(); fig,axs = pplt.subplots(ncols=1,axwidth=1.5)
			
			c = axs[1].pcolormesh(lon,lat,prc[:,:,it]',levels=0:0.2:2,extend="both",cmap="Blues")
		
			axs[1].format(ultitle="(a) 7-day average")
			for ax in axs
				ax.plot(x,y,lw=0.5,c="k")
				ax.scatter(infody[:,2].-360,infody[:,3],c="r",s=2)
				ax.format(xlim=(-80,-75),ylim=(3,8),suptitle=dt[it])
			end
		
			fig.colorbar(c)
			fig.savefig(joinpath("testrainera","$it.png"),transparent=false,dpi=120)
		end
	end
end

# ╔═╡ eaa40336-ab8c-489e-8b36-df9be729998b
load(joinpath("testrainera","5.png"))

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╠═bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╠═c7f79feb-b3f8-4154-a2f3-0c15d48bf61d
# ╠═a15013a7-d990-40db-b608-a18035dcb437
# ╠═d447fdfd-b5c5-44f7-a2a2-161d5bb157b5
# ╠═43bb0b2c-995a-4829-a30e-b953651305de
# ╠═6bd91e94-ac11-4867-bbbb-b534707fce64
# ╠═d7eec81d-2bf8-4eb8-808d-d6fbf0fd3fd4
# ╠═92fd5032-baf7-41f3-bdb8-d0e273171923
# ╠═839bf2a5-0f04-45d1-b683-99f275619ff5
# ╠═eaa40336-ab8c-489e-8b36-df9be729998b
