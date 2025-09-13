### A Pluto.jl notebook ###
# v0.20.5

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
# 04a. Moisture Budget from WRF
"

# ╔═╡ 6ecf692f-a83b-4a97-abe0-5b9b0dd768c2
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

# ╔═╡ afdf0526-f170-4ede-a4de-30cd6a12aefe
begin
	infoall = stninfoall()
	infody  = stninfody(); ndystn = size(infody,1)
	infomo  = stninfomo()
	md"Loading station location information ..."
end

# ╔═╡ 5045b20b-925f-46a5-a6f6-d8cf20e2d79e
begin
	pplt.close()
	f1,a1 = pplt.subplots(nrows=3,ncols=4,axwidth=1.2,wspace=1.5,hspace=1.5)

	xbin = -5:0.2:10
	ybin = -5:0.2:10
	lvls = vcat(0.1,0.5,1,2:2:10,15,20) .* 100

	for istn = 1 : 12

		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = ds["P"][:]
		e  = ds["E"][:]
		∇  = ds["∇"][:] * 3600
		Δ  = ds["ΔWVP"][:]
		close(ds)

		h = fit(Histogram,(p.+Δ,-∇),(xbin,ybin))
		wgts[:,:] += h.weights

		for ibox = 1 : 4

			geo = GeoRegion("OTREC_wrf_stn$(stnstr)_box$ibox",path=srcdir())

			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
			p  = ds["P"][:]
			e  = ds["E"][:]
			∇  = ds["∇"][:] * 3600
			Δ  = ds["ΔWVP"][:]
			close(ds)

			h = fit(Histogram,(p.+Δ,-∇),(xbin,ybin))
			wgts[:,:] += h.weights

		end

		global c = a1[istn].pcolormesh(xbin,ybin,wgts',levels=lvls,extend="both")
		a1[istn].plot([-6,9],[-6,9],c="k",lw=1,linestyle="--")
		a1[istn].plot([-6,9],[-6,9].+0.5*sqrt(2),c="grey",lw=1,linestyle=":")
		a1[istn].plot([-6,9],[-6,9].-0.5*sqrt(2),c="grey",lw=1,linestyle=":")
		a1[istn].format(
			suptitle="(a) Moisture Budget Balancing",
			ultitle=infody[istn,1],
			xlim=(-6,6),#xlocator=-20:10:70,
			ylim=(-6,6),#ylocator=-20:10:70,
			xlabel=L"P - E + $\Delta$ / kg m$^{-2}$ hr$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ hr$^{-1}$"
		)

	end

	f1.colorbar(c,label="Number of Observations")
	f1.savefig(plotsdir("04a-wrfqbudget.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget.png"))
end

# ╔═╡ 7fdb55f7-2888-4937-8409-88331c8c53c4
begin
	pplt.close()
	f2,a2 = pplt.subplots(nrows=3,ncols=4,axwidth=1.2,wspace=1.5,hspace=1.5)

	xbin07 = -1:0.1:5
	ybin07 = -1:0.1:5

	for istn = 1 : 12

		wgts = zeros(length(xbin07)-1,length(ybin07)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=7)
		e  = smooth(ds["E"][:],days=7)
		∇  = smooth(ds["∇"][:],days=7) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=7)
		close(ds)

		h = fit(Histogram,(p,-∇),(xbin07,ybin07))
		wgts[:,:] += h.weights

		for ibox = 1 : 4

			geo = GeoRegion("OTREC_wrf_stn$(stnstr)_box$ibox",path=srcdir())

			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
			p  = smooth(ds["P"][:],days=7)
			e  = smooth(ds["E"][:],days=7)
			∇  = smooth(ds["∇"][:],days=7) * 3600
			Δ  = smooth(ds["ΔWVP"][:],days=7)
			close(ds)

			h = fit(Histogram,(p,-∇),(xbin07,ybin07))
			wgts[:,:] += h.weights

		end

		global c2 = a2[istn].pcolormesh(xbin07,ybin07,wgts',levels=lvls,extend="both")
		a2[istn].plot([-2,5],[-2,5],c="k",lw=1,linestyle="--")
		a2[istn].plot([-2,5],[-2,5].+0.2*sqrt(2),c="grey",lw=1,linestyle=":")
		a2[istn].plot([-2,5],[-2,5].-0.2*sqrt(2),c="grey",lw=1,linestyle=":")
		a2[istn].format(
			suptitle="(b) Moisture Budget Balancing (7 Days Smoothing)",
			ultitle=infody[istn,1],
			xlim=(-1,4),#xlocator=-20:10:70,
			ylim=(-1,4),#ylocator=-20:10:70,
			xlabel=L"P / kg m$^{-2}$ hr$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ hr$^{-1}$"
		)

	end

	f2.colorbar(c2,label="Number of Observations")
	f2.savefig(plotsdir("04a-wrfqbudget-smooth07.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget-smooth07.png"))
end

# ╔═╡ f4ce33b4-a0f5-4db1-8b0c-85b6929e7417
begin
	pplt.close()
	f3,a3 = pplt.subplots(nrows=3,ncols=4,axwidth=1.2,wspace=1.5,hspace=1.5)

	xbin30 = -1:0.05:5
	ybin30 = -1:0.05:5

	for istn = 1 : 12

		wgts = zeros(length(xbin30)-1,length(ybin30)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=30)
		e  = smooth(ds["E"][:],days=30)
		∇  = smooth(ds["∇"][:],days=30) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=30)
		close(ds)

		h = fit(Histogram,(p,-∇),(xbin30,ybin30))
		wgts[:,:] += h.weights

		for ibox = 1 : 4

			geo = GeoRegion("OTREC_wrf_stn$(stnstr)_box$ibox",path=srcdir())

			ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
			p  = smooth(ds["P"][:],days=30)
			e  = smooth(ds["E"][:],days=30)
			∇  = smooth(ds["∇"][:],days=30) * 3600
			Δ  = smooth(ds["ΔWVP"][:],days=30)
			close(ds)

			h = fit(Histogram,(p,-∇),(xbin30,ybin30))
			wgts[:,:] += h.weights

		end

		global c3 = a3[istn].pcolormesh(xbin30,ybin30,wgts',levels=lvls,extend="both")
		a3[istn].plot([-200,300],[-200,300],c="k",lw=1,linestyle="--")
		a3[istn].plot([-200,300],[-200,300].+0.1*sqrt(2),c="grey",lw=1,linestyle=":")
		a3[istn].plot([-200,300],[-200,300].-0.1*sqrt(2),c="grey",lw=1,linestyle=":")
		a3[istn].format(
			suptitle="Moisture Budget Balancing (30 Days Smoothing)",
			xlim=(-0.5,3),#xlocator=-20:10:70,
			ylim=(-0.5,3),#ylocator=-20:10:70,
			xlabel=L"P / kg m$^{-2}$ hr$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ hr$^{-1}$"
		)

	end

	f3.colorbar(c3)
	f3.savefig(plotsdir("04a-wrfqbudget-smooth30.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget-smooth30.png"))
end

# ╔═╡ 55121124-de02-4b81-86e0-08ba1b8db1dc
begin
	Δmat = zeros(12,5)
	for istn = 1 : 12
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=7)
		e  = smooth(ds["E"][:],days=7)
		∇  = smooth(ds["∇"][:],days=7) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=7)
		close(ds)
	
		r = sqrt.((p.-e.+Δ.+∇).^2)/2
		Δmat[istn,1] = mean(r[.!isnan.(r)])
		
	end
	for istn = 1 : 12, ibox = 1 : 4
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		boxstr = @sprintf("%d",ibox)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=7)
		e  = smooth(ds["E"][:],days=7)
		∇  = smooth(ds["∇"][:],days=7) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=7)
		close(ds)
	
		r = sqrt.((p.-e.+Δ.+∇).^2)/2
		Δmat[istn,ibox+1] = mean(r[.!isnan.(r)])
		
	end
	open(datadir("wrf3","qbudgetdiff-smooth07.txt"),"w") do io
		writedlm(io,Δmat)
	end
end

# ╔═╡ 617614bb-c457-418b-a31d-2a27ccb0199d
begin
	Δmat .= 0
	for istn = 1 : 12
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=30)
		e  = smooth(ds["E"][:],days=30)
		∇  = smooth(ds["∇"][:],days=30) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=30)
		close(ds)
	
		r = sqrt.((p.-e.+Δ.+∇).^2)/2
		Δmat[istn,1] = mean(r[.!isnan.(r)])
		
	end
	for istn = 1 : 12, ibox = 1 : 4
	
		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		boxstr = @sprintf("%d",ibox)
		geo    = GeoRegion("OTREC_wrf_stn$(stnstr)_box$(boxstr)",path=srcdir())
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=30)
		e  = smooth(ds["E"][:],days=30)
		∇  = smooth(ds["∇"][:],days=30) * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=30)
		close(ds)
	
		r = sqrt.((p.-e.+Δ.+∇).^2)/2
		Δmat[istn,ibox+1] = mean(r[.!isnan.(r)])
		
	end
	open(datadir("wrf3","qbudgetdiff-smooth30.txt"),"w") do io
		writedlm(io,Δmat)
	end
end

# ╔═╡ 5696050b-6f56-4d95-b315-ebaf1b4e4f55
Δmat .> 0.05

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─6ecf692f-a83b-4a97-abe0-5b9b0dd768c2
# ╠═afdf0526-f170-4ede-a4de-30cd6a12aefe
# ╟─5045b20b-925f-46a5-a6f6-d8cf20e2d79e
# ╟─7fdb55f7-2888-4937-8409-88331c8c53c4
# ╟─f4ce33b4-a0f5-4db1-8b0c-85b6929e7417
# ╟─55121124-de02-4b81-86e0-08ba1b8db1dc
# ╟─617614bb-c457-418b-a31d-2a27ccb0199d
# ╠═5696050b-6f56-4d95-b315-ebaf1b4e4f55
