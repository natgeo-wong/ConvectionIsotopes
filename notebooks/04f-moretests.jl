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
	using DelimitedFiles
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
# 03b. Station W-Weighted Pressure, $\sigma$
"

# ╔═╡ 59c930cd-5b7f-4047-8660-615148d1bd9f
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ 4319fd0e-fd9f-424e-9286-3b3b5a844b73
function extract(geoname,iso,days)

	dystr  = "-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-p_wwgt-daily$(dystr).nc"
	))
	p_wwgt = ds["p_wwgt"][:] / 100
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-$(iso)ωdqdp-daily$(dystr).nc"
	))
	ωdqdp = ds["⟨$(iso)ωdqdp⟩"][:]
	ω = ds["⟨ω⟩"][:]
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-$(iso)ωdqdpfixed-daily$(dystr).nc"
	))
	ωdqdp_f = ds["⟨$(iso)ωdqdp⟩"][:]
	ω_f = ds["⟨ω⟩"][:]
	close(ds)

	return p_wwgt,ωdqdp*1e7 ./ω,ωdqdp_f*1e7 ./ω_f
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function plotωdqdp(
	axes,ii,dpwwgtbin,dωdqdpbin;
	ID, days=0, box = false, bnum = 1, cinfo = false
)

	binHDO = zeros(length(dωdqdpbin)-1,length(dpwwgtbin)-1)
	binO18 = zeros(length(dωdqdpbin)-1,length(dpwwgtbin)-1)
	binHDO_f = zeros(length(dωdqdpbin)-1,length(dpwwgtbin)-1)
	binO18_f = zeros(length(dωdqdpbin)-1,length(dpwwgtbin)-1)
	dωdqdpbinplt = (dωdqdpbin[1:(end-1)] .+ dωdqdpbin[2:end])/2
	dpwwgtbinplt = (dpwwgtbin[1:(end-1)] .+ dpwwgtbin[2:end])/2

	for stn in ID
		stnstr = @sprintf("%02d",stn)
		if box
			for ibox = 1 : bnum
				boxstr = @sprintf("%02d",ibox)
				geoname = "OTREC_STN$(stnstr)_$(boxstr)"
				pwwgt,ωdqdp,ωdqdp_f = extract(geoname,"HDO_",days)
				binHDO += fit(Histogram,(ωdqdp,pwwgt),(dωdqdpbin/8,dpwwgtbin)).weights
				binHDO_f += fit(Histogram,(ωdqdp_f,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
				pwwgt,ωdqdp,ωdqdp_f = extract(geoname,"O18_",days)
				binO18 += fit(Histogram,(ωdqdp,pwwgt),(dωdqdpbin/8,dpwwgtbin)).weights
				binO18_f += fit(Histogram,(ωdqdp_f,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
			end
		else
			geoname = "OTREC_STN$stnstr"
			pwwgt,ωdqdp,ωdqdp_f = extract(geoname,"HDO_",days)
			binHDO += fit(Histogram,(ωdqdp,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
			binHDO_f += fit(Histogram,(ωdqdp_f,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
			pwwgt,ωdqdp,ωdqdp_f = extract(geoname,"O18_",days)
			binO18 += fit(Histogram,(ωdqdp,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
			binO18_f += fit(Histogram,(ωdqdp_f,pwwgt),(dωdqdpbin,dpwwgtbin)).weights
		end
	end
	c1 = axes[4*ii].pcolormesh(dωdqdpbinplt,dpwwgtbinplt,(binHDO_f./sum(binHDO_f))'*100,extend="both")
	c2 = axes[4*ii-2].pcolormesh(dωdqdpbinplt,dpwwgtbinplt,(binO18_f./sum(binO18_f))'*100,extend="both")
	c3 = axes[4*ii-1].pcolormesh(dωdqdpbinplt/8,dpwwgtbinplt,(binHDO./sum(binHDO))'*100,extend="both")
	c4 = axes[4*ii-3].pcolormesh(dωdqdpbinplt/8,dpwwgtbinplt,(binO18./sum(binO18))'*100,extend="both")

	if cinfo
		return c1,c2,c3,c4
	else
		return
	end

end

# ╔═╡ 1343fbae-0ebd-4237-8273-0ebab8325424
function axesformat!(axes)

	for ax in axes
		ax.format(
			ylim=(1000,0),ylabel=L"$p_\omega$ / hPa",
			# xlabel=L"$\partial_p(q_h/q)$ / 10$^{-6}$ Pa$^{-1}$",
			suptitle="7-Day Moving Average"
		)
	end

	# for ii = 1 : 8
	# 	axes[2*ii].format(xlim=(-2,10),lrtitle="HDO")
	# 	axes[2*ii-1].format(xlim=(-2,10)./10,lrtitle="O18")
	# end

	# axes[2].format(ultitle="(a) All Stations")
	# axes[2].format(ultitle="(a) Colombia")
	# axes[4].format(ultitle="(b) San Andres")
	# axes[6].format(ultitle="(c) Buenaventura")
	# axes[6].text(0.48,39.5,"Bahia Solano",fontsize=10)
	# axes[8].format(ultitle="(d) Quibdo")
	# axes[10].format(ultitle="(e) Costa Rica")
	# axes[12].format(ultitle="(f) EEFMB")
	# axes[12].text(0.41,39.5,"ADMQ",fontsize=10)
	# axes[12].text(0.41,34.5,"CGFI",fontsize=10)
	# axes[14].format(ultitle="(g) Cahuita")
	# axes[14].text(0.48,39.5,"Bataan",fontsize=10)
	# axes[14].text(0.48,34.5,"Limon",fontsize=10)
	# axes[16].format(ultitle="(h) Liberia")
	# axes[16].text(0.48,39.5,"OSA",fontsize=10)

	return

end

# ╔═╡ 8c211620-d632-4f23-85f5-a702faf82270
begin
	pplt.close(); fig,axs = pplt.subplots(
		[
			[04,03,08,07,12,11,16,15],
			[02,01,06,05,10,09,14,13],
			[20,19,24,23,28,27,32,31],
			[18,17,22,21,26,25,30,29],
		],aspect=1,axwidth=0.6,
		wspace=[0,1.5,0,1.5,0,1.5,0],
		hspace=[0,1.5,0,]
	)

	c1,_,_,_ = 
	plotωdqdp(axs,1,0:20:1000,0:0.2:8,ID=1:4,box=true,bnum=4,days=7,cinfo=true)
	plotωdqdp(axs,2,0:20:1000,0:0.2:8,ID=1,box=true,bnum=4,days=7)
	plotωdqdp(axs,3,0:20:1000,0:0.2:8,ID=3:4,box=true,bnum=4,days=7)
	plotωdqdp(axs,4,0:20:1000,0:0.2:8,ID=2,box=true,bnum=4,days=7)
	plotωdqdp(axs,5,0:20:1000,0:0.2:8,ID=5:12,box=true,bnum=4,days=7)
	plotωdqdp(axs,6,0:20:1000,0:0.2:8,ID=5:7,box=true,bnum=4,days=7)
	plotωdqdp(axs,7,0:20:1000,0:0.2:8,ID=9:11,box=true,bnum=4,days=7)
	plotωdqdp(axs,8,0:20:1000,0:0.2:8,ID=[8,12],box=true,bnum=4,days=7)
	
	axesformat!(axs)

	# fig.colorbar(c1,length=0.75)
	fig.savefig(plotsdir("04f-moretests.png"),transparent=false,dpi=400)
	load(plotsdir("04f-moretests.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╠═4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╠═5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╠═1343fbae-0ebd-4237-8273-0ebab8325424
# ╠═8c211620-d632-4f23-85f5-a702faf82270
