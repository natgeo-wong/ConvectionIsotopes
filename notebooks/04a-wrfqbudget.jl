### A Pluto.jl notebook ###
# v0.19.32

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
	using NCDatasets
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

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
		hrs = days * 8
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

# ╔═╡ 5045b20b-925f-46a5-a6f6-d8cf20e2d79e
begin
	pplt.close()
	f1,a1 = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)

	xbin = -20:1:50
	ybin = -20:1:50

	for istn = 1 : 12

		wgts = zeros(length(xbin)-1,length(ybin)-1)

		for ibox = 1 : 4
			stnstr = @sprintf("%02d",istn)
			boxstr = @sprintf("%02d",ibox)
			gname  = "OTREC_STN$(stnstr)_$(boxstr)"
			
			ds = NCDataset(datadir("wrf","processed","$(gname)-QBUDGET.nc"))
			p  = smooth(ds["P"][:],days=30)    / 3 * 24
			e  = smooth(ds["E"][:],days=30) * 3600 * 24
			∇  = smooth(ds["∇"][:],days=30) * 3600 * 24
			Δ  = smooth(ds["ΔWVP"][:],days=30) / 3 * 24
			close(ds)
			ds = NCDataset(datadir("wrf","processed","$(gname)-∇decompose.nc"))
			A  = smooth(ds["ADV"][:],days=30) * 3600 * 24
			D  = smooth(ds["DIV"][:],days=30) * 3600 * 24
			close(ds)

			h = fit(Histogram,(p,e.-∇.-Δ),(xbin,ybin))
			wgts[:,:] += h.weights
		end

		a1[istn].pcolormesh(xbin,ybin,wgts',levels=vcat(0:10)*10,extend="both")
		a1[istn].plot([-100,150],[-100,150],c="k",lw=1,linestyle=":")
		a1[istn].format(
			suptitle="Moisture Budget Balancing",
			xlim=(-25,25),xlocator=-20:10:50,
			xlabel=L"P / kg m$^{-2}$ day$^{-1}$",
			ylim=(-25,25),ylocator=-20:10:50,
			ylabel=L"E $-\nabla - \Delta$ / kg m$^{-2}$ day$^{-1}$"
		)

	end
	
	f1.savefig(plotsdir("04a-wrfqbudget.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget.png"))
end

# ╔═╡ 5074b541-5283-4b8d-9c2d-b93fe909fa93
begin
	pplt.close()
	fig,axs = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)

	for istn = 1 : 12

		stnstr = @sprintf("%02d",istn)
		ds = NCDataset(datadir("wrf","processed","OTREC_STN$stnstr-QBUDGET.nc"))
		p  = smooth(ds["P"][:],days=30)    / 3 * 24
		e  = smooth(ds["E"][:],days=30) * 3600 * 24
		∇  = smooth(ds["∇"][:],days=30) * 3600 * 24
		close(ds)

		h = fit(Histogram,(p,e.-∇),(xbin,ybin))

		axs[istn].pcolormesh(xbin,ybin,h.weights',levels=0:2:20,extend="both")
		axs[istn].plot([-100,100],[-100,100],c="k",lw=1,linestyle=":")
		axs[istn].format(
			suptitle="Moisture Budget Balancing",
			xlim=(-5,50),xlocator=0:10:50,
			xlabel=L"P / kg m$^{-2}$ day$^{-1}$",
			ylim=(-5,50),ylocator=0:10:50,
			ylabel=L"E $-\nabla - \Delta$ / kg m$^{-2}$ day$^{-1}$"
		)

	end
	
	fig.savefig(plotsdir("04a-wrfqbudget2.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget2.png"))
end

# ╔═╡ 37d37e56-528a-40f1-bd3f-e1a44e397b43
begin
	pplt.close()
	f2,a2 = pplt.subplots(aspect=6,axwidth=4,nrows=3)

	hbin = -50 : 0.2 : 50

	for istn = 1 : 1

		stnstr = @sprintf("%02d",istn)
		ds = NCDataset(datadir("wrf","processed","OTREC_STN$stnstr-QBUDGET.nc"))
		∇  = smooth(ds["∇"][:],days=7) * 86400
		p  = smooth(ds["P"][:],days=7) * 8
		close(ds)

		h = fit(Histogram,∇,hbin)

		a2[1].plot((hbin[1:(end-1)] .+ hbin[2:end])/2,h.weights)
		a2[1].plot(ones(2)*mean(∇),[0,250],c="k",lw=1)
		a2[1].plot(ones(2)*mean(p),[0,250],c="k",lw=1,linestyle="--")

	end

	for istn = 2 : 4

		stnstr = @sprintf("%02d",istn)
		ds = NCDataset(datadir("wrf","processed","OTREC_STN$stnstr-QBUDGET.nc"))
		∇  = smooth(ds["∇"][:],days=7) * 86400
		p  = smooth(ds["P"][:],days=7) * 8
		close(ds)

		h = fit(Histogram,∇,hbin)

		a2[2].plot((hbin[1:(end-1)] .+ hbin[2:end])/2,h.weights)
		a2[2].plot(ones(2)*mean(∇),[0,250],c="k",lw=1)
		a2[2].plot(ones(2)*mean(p),[0,250],c="k",lw=1,linestyle="--")

	end

	for istn = 5 : 12

		stnstr = @sprintf("%02d",istn)
		ds = NCDataset(datadir("wrf","processed","OTREC_STN$stnstr-QBUDGET.nc"))
		∇  = smooth(ds["∇"][:],days=7) * 86400
		p  = smooth(ds["P"][:],days=7) * 8
		close(ds)

		h = fit(Histogram,∇,hbin)

		a2[3].plot((hbin[1:(end-1)] .+ hbin[2:end])/2,h.weights)
		a2[3].plot(ones(2)*mean(∇),[0,250],c="k",lw=1)
		a2[3].plot(ones(2)*mean(p),[0,250],c="k",lw=1,linestyle="--")

	end

	for axs in a2
		axs.format(xlabel=L"$\nabla$ / kg m$^{-2}$ day$^{-1}$",ylim=(0,50),xlim=(-50,50))
	end
	
	f2.savefig(plotsdir("04a-∇histogram.png"),transparent=false,dpi=400)
	load(plotsdir("04a-∇histogram.png"))
end

# ╔═╡ 3ef863fa-d34a-48a5-be60-f4f4527146ed
begin
	pplt.close()
	f3,a3 = pplt.subplots(axwidth=2,wspace=1.5,hspace=1.5)

	for istn = 1 : 12

		stnstr = @sprintf("%02d",istn)

		if istn == 1
			clr = "b"
			boxvec = 1 : 4
		elseif istn >= 2 && istn <= 4
			clr = "cyan"
			boxvec = [1,3]
		elseif istn >= 5 && istn <= 7
			clr = "r"
			boxvec = []
		elseif istn >= 9 && istn <= 11
			clr = "yellow4"
			boxvec = [4]
		else
			clr = "orange"
			boxvec = [1]
		end

		xx = 0
		yy = 0

		for ibox = 1 : 4
			boxstr = @sprintf("%02d",ibox)
			gname  = "OTREC_STN$(stnstr)_$(boxstr)"
			ds = NCDataset(datadir("wrf","processed","$(gname)-QBUDGET.nc"))
			p  = ds["P"][:] * 8
			e  = ds["E"][:] * 86400
			Δ  = ds["ΔWVP"][:] * 8
			∇  = ds["∇"][:] * 86400
			close(ds)
			ds = NCDataset(datadir("wrf","processed","$(gname)-∇decompose.nc"))
			A  = ds["ADV"][:] * 86400
			D  = ds["DIV"][:] * 86400
			close(ds)

			xx += mean(e.-∇) / 4
			yy += mean(p) / 4
	
			a3[1].scatter(mean(-∇),mean(p),c=clr)
			a3[1].scatter(mean(-∇.+A),mean(p.+A),c=clr,s=10)
		end
		gname  = "OTREC_STN$(stnstr)"
		ds = NCDataset(datadir("wrf","processed","$(gname)-QBUDGET.nc"))
		p  = ds["P"][:] * 8
		e  = ds["E"][:] * 86400
		Δ  = ds["ΔWVP"][:] * 8
		∇  = ds["∇"][:] * 86400
		close(ds)
		ds = NCDataset(datadir("wrf","processed","$(gname)-∇decompose.nc"))
			A  = ds["ADV"][:] * 86400
			D  = ds["DIV"][:] * 86400
		close(ds)
		
		# a3[1].scatter(mean(-∇.+A),mean(p.+A),c=clr,s=10)
		# a3[1].scatter(mean(-∇),mean(p),c=clr)
		# a3[1].scatter(xx,yy,c=clr)
		

	end

	a3[1].plot([-100,100],[-100,100],c="k",lw=1,linestyle=":")
	a3[1].format(
		suptitle="Moisture Budget Balancing",
		xlim=(-75,100),xlocator=-100:25:100,
		xlabel=L"$\mu(E-\nabla)$ / kg m$^{-2}$ day$^{-1}$",
		ylim=(-75,100),ylocator=-100:25:100,
		ylabel=L"$\mu(P)$ / kg m$^{-2}$ day$^{-1}$"
	)
	
	f3.savefig(plotsdir("04a-wrfqbudget-mean.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget-mean.png"))
end

# ╔═╡ b7a241dc-4f92-4952-9c14-6311d7f36cad
begin
	pplt.close()
	f4,a4 = pplt.subplots(nrows=3,ncols=4,axwidth=1,wspace=1.5,hspace=1.5)

	xbin1 = -3:0.1:3
	ybin1 = -3:0.1:3

	for istn = 1 : 12

		stnstr = @sprintf("%02d",istn)
		ds = NCDataset(datadir("wrf","processed","OTREC_STN$stnstr-QBUDGET.nc"))
		p  = ds["P"][:] / 10800 * 1000
		e  = ds["E"][:] * 1000
		Δ  = ds["ΔWVP"][:] / 10800 * 1000
		∇  = ds["∇"][:] * 1000
		close(ds)

		h = fit(Histogram,(e.-∇,Δ.+p),(xbin1,ybin1))

		a4[istn].pcolormesh(xbin1,ybin1,h.weights',levels=vcat(0,1,2,3,5,7,10,15,20,30,50,70,100))
		a4[istn].plot([-20,20],[-20,20],c="k",lw=1,linestyle=":")
		a4[istn].format(
			suptitle="Moisture Budget Balancing",
			xlim=(-1.5,1.5),xlocator=-2:2,
			xlabel=L"E - $\nabla$ (Import) / g m$^{-2}$ s$^{-1}$",
			ylim=(-1.5,1.5),ylocator=-2:2,
			ylabel=L"P + $\Delta$ (Export) / g m$^{-2}$ s$^{-1}$"
		)

	end
	
	f4.savefig(plotsdir("04a-wrfqbudget-invsout.png"),transparent=false,dpi=400)
	load(plotsdir("04a-wrfqbudget-invsout.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╠═6ecf692f-a83b-4a97-abe0-5b9b0dd768c2
# ╟─5045b20b-925f-46a5-a6f6-d8cf20e2d79e
# ╟─5074b541-5283-4b8d-9c2d-b93fe909fa93
# ╟─37d37e56-528a-40f1-bd3f-e1a44e397b43
# ╟─3ef863fa-d34a-48a5-be60-f4f4527146ed
# ╠═b7a241dc-4f92-4952-9c14-6311d7f36cad
