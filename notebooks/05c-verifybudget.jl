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

# ╔═╡ aaa79ea7-d5ed-457e-9b99-57a81b6c3d95
prfx = "CrossITCZ"

# ╔═╡ 5045b20b-925f-46a5-a6f6-d8cf20e2d79e
begin
	pplt.close()
	f1,a1 = pplt.subplots(nrows=5,ncols=5,axwidth=0.8,wspace=1.5,hspace=1.5)

	xbin = -150:2.5:150
	ybin = -150:2.5:150
	lvls = vcat(0.5,1,2:2:10,15,20) .* 10

	for istn = 1 : 25

		wgts = zeros(length(xbin)-1,length(ybin)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_$(prfx)$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = ds["P"][:] * 24
		e  = ds["E"][:] * 24
		∇  = ds["∇"][:] * 24 * 3600
		Δ  = ds["ΔWVP"][:] * 24
		close(ds)

		h = fit(Histogram,(p.-e.+Δ,-∇),(xbin,ybin))
		wgts[:,:] += h.weights

		a1[istn].pcolormesh(xbin,ybin,wgts',levels=lvls,extend="both")
		a1[istn].plot([-100,150],[-100,150],c="k",lw=1,linestyle=":")
		a1[istn].format(
			suptitle="Moisture Budget Balancing | $(prfx)",ultitle="$istn",
			xlim=(-60,60),#xlocator=-20:10:70,
			ylim=(-60,60),#ylocator=-20:10:70,
			xlabel=L"P - E + $\Delta$ / kg m$^{-2}$ day$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ day$^{-1}$"
		)

	end

	c1 = a1[1].pcolormesh([0,1],[0,1],zeros(2,2)*NaN,levels=lvls,extend="both")

	f1.colorbar(c1,row=[2,4])
	f1.savefig(plotsdir("05c-wrfqbudget-$(prfx).png"),transparent=false,dpi=400)
	load(plotsdir("05c-wrfqbudget-$(prfx).png"))
end

# ╔═╡ 1f30ec64-344a-4545-81b9-05c52fd14bc4
begin
	pplt.close()
	f2,a2 = pplt.subplots(nrows=5,ncols=5,axwidth=0.8,wspace=1.5,hspace=1.5)

	xbin07 = -150:2.5:150
	ybin07 = -150:2.5:150

	for istn = 1 : 25

		wgts = zeros(length(xbin07)-1,length(ybin07)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_$(prfx)$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=7) * 24
		e  = smooth(ds["E"][:],days=7) * 24
		∇  = smooth(ds["∇"][:],days=7) * 24 * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=7) * 24
		close(ds)

		h = fit(Histogram,(p.-e.+Δ,-∇),(xbin07,ybin07))
		wgts[:,:] += h.weights

		a2[istn].pcolormesh(xbin07,ybin07,wgts',levels=lvls,extend="both")
		a2[istn].plot([-100,150],[-100,150],c="k",lw=1,linestyle=":")
		a2[istn].format(
			suptitle="Moisture Budget Balancing | $(prfx) | 7 Days Smoothing",
			ultitle="$istn",
			xlim=(-20,120),#xlocator=-20:10:70,
			ylim=(-20,120),#ylocator=-20:10:70,
			xlabel=L"P - E + $\Delta$ / kg m$^{-2}$ day$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ day$^{-1}$"
		)

	end

	c2 = a2[1].pcolormesh([0,1],[0,1],zeros(2,2)*NaN,levels=lvls,extend="both")

	f2.colorbar(c1,row=[2,4])
	f2.savefig(plotsdir("05c-wrfqbudget-$(prfx)-smooth07.png"),transparent=false,dpi=400)
	load(plotsdir("05c-wrfqbudget-$(prfx)-smooth07.png"))
end

# ╔═╡ e014c072-b43f-47fb-9b40-093ef8f8df7a
begin
	pplt.close()
	f3,a3 = pplt.subplots(nrows=5,ncols=5,axwidth=0.8,wspace=1.5,hspace=1.5)

	xbin30 = -20:2:100
	ybin30 = -20:2:100

	for istn = 1 : 25

		wgts = zeros(length(xbin30)-1,length(ybin30)-1)
		stnstr = @sprintf("%02d",istn)
		geo    = GeoRegion("OTREC_wrf_$(prfx)$(stnstr)",path=srcdir())
			
		ds  = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
		p  = smooth(ds["P"][:],days=30) * 24
		e  = smooth(ds["E"][:],days=30) * 24
		∇  = smooth(ds["∇"][:],days=30) * 24 * 3600
		Δ  = smooth(ds["ΔWVP"][:],days=30) * 24
		close(ds)

		h = fit(Histogram,(p.-e.+Δ,-∇),(xbin30,ybin30))
		wgts[:,:] += h.weights

		a3[istn].pcolormesh(xbin30,ybin30,wgts',levels=lvls,extend="both")
		a3[istn].plot([-100,150],[-100,150],c="k",lw=1,linestyle=":")
		a3[istn].format(
			suptitle="Moisture Budget Balancing | $(prfx) | 30 Days Smoothing",
			ultitle="$istn",
			xlim=(-10,60),#xlocator=-20:10:70,
			ylim=(-10,60),#ylocator=-20:10:70,
			xlabel=L"P - E + $\Delta$ / kg m$^{-2}$ day$^{-1}$",
			ylabel=L"$-\nabla$ / kg m$^{-2}$ day$^{-1}$"
		)

	end

	c3 = a3[1].pcolormesh([0,1],[0,1],zeros(2,2)*NaN,levels=lvls,extend="both")

	f3.colorbar(c3,row=[2,4])
	f3.savefig(plotsdir("05c-wrfqbudget-$(prfx)-smooth30.png"),transparent=false,dpi=400)
	load(plotsdir("05c-wrfqbudget-$(prfx)-smooth30.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╠═6ecf692f-a83b-4a97-abe0-5b9b0dd768c2
# ╠═aaa79ea7-d5ed-457e-9b99-57a81b6c3d95
# ╟─5045b20b-925f-46a5-a6f6-d8cf20e2d79e
# ╟─1f30ec64-344a-4545-81b9-05c52fd14bc4
# ╟─e014c072-b43f-47fb-9b40-093ef8f8df7a
