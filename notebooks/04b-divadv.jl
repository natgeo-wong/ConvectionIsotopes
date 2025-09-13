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
# 04b. Divergence vs Advection
"

# ╔═╡ 90e7029c-6756-4f80-9f63-7b05304dfe49
begin
	infody  = stninfody()
	md"Loading station location information ..."
end

# ╔═╡ a16e2a7f-2c8c-4b4a-9a9e-cb91fd80fdfb
function smooth(data::AbstractArray;days=0)

	if !iszero(days)
		hrs = days * 24
		nt = length(data) - hrs
		ndata = zeros(nt)
	
		for it = 1 : nt
			for ii = 0 : (hrs-1)
				ndata[it] += data[it+ii]  / hrs
			end
		end
	else
		ndata = deepcopy(data)
	end

	return ndata

end

# ╔═╡ 4d2d13e0-178c-4a96-8afa-988a84515373
begin
	xbin =-20:0.5:20
	ybin = -20:0.5:20
end

# ╔═╡ 5045b20b-925f-46a5-a6f6-d8cf20e2d79e
begin
	pplt.close()
	f1,a1 = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)
	
	for istn = 1 : 12

		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
	
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = ds["DIV"][:] * 3600
		adv = ds["ADV"][:] * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = ds["∇"][:] * 3600
		close(ds)

		h = fit(Histogram,(∇,div.+adv),(xbin,ybin))
		wgts[:,:] += h.weights

		for boxx = 1 : 4
			boxstr = @sprintf("%d",boxx)
			geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
			div = ds["DIV"][:] * 3600
			adv = ds["ADV"][:] * 3600
			close(ds)
			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
			∇   = ds["∇"][:] * 3600
			close(ds)

			h = fit(Histogram,(∇,div.+adv),(xbin,ybin))
			wgts[:,:] += h.weights
		end

		c1 = a1[istn].pcolormesh(xbin,ybin,wgts',levels=[0,1,2,5,10,20,50,100,200,500]*10,extend="both")
		a1[istn].plot([-20,20],[-20,20],c="k",lw=1,linestyle=":")
		a1[istn].format(
			ultitle="$(infody[istn,1])",suptitle=L"Calculated vs Given (0.5$\degree$)",
			xlim=(-12,12),xlocator=-20:5:20,
			xlabel=L"$\nabla\cdot(uq)$ / kg m$^{-2}$ hr$^{-1}$",
			ylim=(-12,12),ylocator=-20:5:20,
			ylabel=L"$u\cdot\nabla q + q\nabla\cdot u$ / kg m$^{-2}$ hr$^{-1}$"
		)

		if istn == 8
			f1.colorbar(c1,length=0.5,label="Number of Observations")
		end

	end
	
	f1.savefig(plotsdir("04b-∇decompose.png"),transparent=false,dpi=400)
	load(plotsdir("04b-∇decompose.png"))
end

# ╔═╡ 9836a8a9-1165-49d7-a84c-a08220ed8315
begin
	pplt.close()
	f2,a2 = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)

	for istn = 1 : 12

		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
	
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = ds["DIV"][:] * 3600
		adv = ds["ADV"][:] * 3600
		close(ds)

		h = fit(Histogram,(div,adv),(xbin,ybin))
		wgts[:,:] += h.weights

		for ibox = 1 : 4
			boxstr = @sprintf("%d",ibox)
			geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
			div = ds["DIV"][:] * 3600
			adv = ds["ADV"][:] * 3600
			close(ds)
	
			h = fit(Histogram,(div,adv),(xbin,ybin))
			wgts[:,:] += h.weights
		end
	
		c2 = a2[istn].pcolormesh(xbin,ybin,wgts',levels=[0,1,2,5,10,20,50,100,200,500]*10,extend="both")
		a2[istn].plot([-20,20],[-20,20],c="k",lw=1,linestyle=":")
		a2[istn].plot([-20,20],[20,-20],c="k",lw=1,linestyle=":")
		a2[istn].format(
			suptitle=L"Decomposing $\nabla$ (0.5$\degree$)",
			ultitle="$(infody[istn,1])",
			xlim=(-12,12),xlocator=-20:5:20,
			xlabel=L"$q\nabla\cdot u$ / kg m$^{-2}$ hr$^{-1}$",
			ylim=(-12,12),ylocator=-20:5:20,
			ylabel=L"$u\cdot\nabla q$ / kg m$^{-2}$ hr$^{-1}$"
		)

		if istn == 8
			f2.colorbar(c2,length=0.5,label="Number of Observations")
		end
		
	end
	
	f2.savefig(plotsdir("04b-divadv.png"),transparent=false,dpi=400)
	load(plotsdir("04b-divadv.png"))
end

# ╔═╡ 77f7d812-ca24-45ed-bbf0-9e80f11ba24f
begin
	xcbin =-2.5:0.1:2.5
	ycbin = -2.5:0.1:2.5
end

# ╔═╡ 40621124-2edf-4828-9145-de8cdc6488ed
begin
	pplt.close()
	f3,a3 = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)
	
	for istn = 1 : 12

		wgts = zeros(length(xcbin)-1,length(ycbin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
	
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = smooth(ds["DIV"][:],days=30) * 3600
		adv = smooth(ds["ADV"][:],days=30) * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = smooth(ds["∇"][:],days=30) * 3600
		close(ds)

		h = fit(Histogram,(∇,div.+adv),(xcbin,ycbin))
		wgts[:,:] += h.weights

		for boxx = 1 : 4
			boxstr = @sprintf("%d",boxx)
			geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
			div = smooth(ds["DIV"][:],days=30) * 3600
			adv = smooth(ds["ADV"][:],days=30) * 3600
			close(ds)
			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
			∇   = smooth(ds["∇"][:],days=30) * 3600
			close(ds)

			h = fit(Histogram,(∇,div.+adv),(xcbin,ycbin))
			wgts[:,:] += h.weights
		end

		c3 = 
		a3[istn].pcolormesh(xcbin,ycbin,wgts',levels=100:100:1000,extend="both")
		a3[istn].plot([-20,20],[-20,20],c="k",lw=1,linestyle=":")
		a3[istn].format(
			ultitle="$(infody[istn,1])",suptitle=L"Calculated vs Given (0.5$\degree$)",
			xlim=(-2.5,2.5),xlocator=-5:5,
			ylim=(-2.5,2.5),ylocator=-5:5,
			xlabel=L"$\nabla\cdot(uq)$ / kg m$^{-2}$ hr$^{-1}$",
			ylabel=L"$u\cdot\nabla q + q\nabla\cdot u$ / kg m$^{-2}$ hr$^{-1}$"
		)

		if istn == 8
			f3.colorbar(c3,length=0.5,label="Number of Observations")
		end

	end
	
	f3.savefig(plotsdir("04b-∇decompose-smooth30.png"),transparent=false,dpi=400)
	load(plotsdir("04b-∇decompose-smooth30.png"))
end

# ╔═╡ 79ace2fa-98eb-4908-8d17-50f244532aa1
begin
	Δ = zeros(12,5)
	for istn = 1 : 12
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = ds["DIV"][:] * 3600
		adv = ds["ADV"][:] * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = ds["∇"][:] * 3600
		close(ds)
	
		r = sqrt.((div.+adv.-∇).^2)/2
		Δ[istn,1] = mean(r[.!isnan.(r)])
		
	end
	for istn = 1 : 12, ibox = 1 : 4
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		boxstr = @sprintf("%d",ibox)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = ds["DIV"][:] * 3600
		adv = ds["ADV"][:] * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = ds["∇"][:] * 3600
		close(ds)
	
		r = sqrt.((div.+adv.-∇).^2)/2
		Δ[istn,ibox+1] = mean(r[.!isnan.(r)])
		
	end
	open(datadir("wrf3","wrfvscalc.txt"),"w") do io
		writedlm(io,Δ)
	end
end

# ╔═╡ 55de4e3e-f2d5-4d31-9679-e4285c93d8f3
begin
	Δ .= 0
	for istn = 1 : 12
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = smooth(ds["DIV"][:],days=30) * 3600
		adv = smooth(ds["ADV"][:],days=30) * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = smooth(ds["∇"][:],days=30) * 3600
		close(ds)
	
		r = sqrt.((div.+adv.-∇).^2)/2
		Δ[istn,1] = mean(r[.!isnan.(r)])
		
	end
	for istn = 1 : 12, ibox = 1 : 4
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		boxstr = @sprintf("%d",ibox)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		div = smooth(ds["DIV"][:],days=30) * 3600
		adv = smooth(ds["ADV"][:],days=30) * 3600
		close(ds)
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		∇   = smooth(ds["∇"][:],days=30) * 3600
		close(ds)
	
		r = sqrt.((div.+adv.-∇).^2)/2
		Δ[istn,ibox+1] = mean(r[.!isnan.(r)])
		
	end
	open(datadir("wrf3","wrfvscalc-smooth30.txt"),"w") do io
		writedlm(io,Δ)
	end
end

# ╔═╡ ad78d452-0bb4-4a80-b0a7-38270386cc65
Δ .> 0.15

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╠═bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─90e7029c-6756-4f80-9f63-7b05304dfe49
# ╟─a16e2a7f-2c8c-4b4a-9a9e-cb91fd80fdfb
# ╟─4d2d13e0-178c-4a96-8afa-988a84515373
# ╟─5045b20b-925f-46a5-a6f6-d8cf20e2d79e
# ╟─9836a8a9-1165-49d7-a84c-a08220ed8315
# ╠═77f7d812-ca24-45ed-bbf0-9e80f11ba24f
# ╟─40621124-2edf-4828-9145-de8cdc6488ed
# ╟─79ace2fa-98eb-4908-8d17-50f244532aa1
# ╟─55de4e3e-f2d5-4d31-9679-e4285c93d8f3
# ╠═ad78d452-0bb4-4a80-b0a7-38270386cc65
