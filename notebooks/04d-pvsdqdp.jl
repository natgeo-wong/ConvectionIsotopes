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
		"$geoname-$(iso)dqdp$(dystr).nc"
	))
	p = ds["P"][:]
	dqdp = ds["$(iso)dqdp"][:,:]
	close(ds)

	return p,dqdp*1e7
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function plotdqdp(
	axes,ii,dqdpbin;
	ID, days=0, box = false, bnum = 1, cinfo = false
)

	binHDO = zeros(length(dqdpbin)-1,50)
	binO18 = zeros(length(dqdpbin)-1,50)
	binplt = (dqdpbin[1:(end-1)] .+ dqdpbin[2:end])/2

	for stn in ID
		stnstr = @sprintf("%02d",stn)
		if box
			for ibox = 1 : bnum
				boxstr = @sprintf("%02d",ibox)
				geoname = "OTREC_STN$(stnstr)_$(boxstr)"
				p,dqdp = extract(geoname,"HDO_",days)
				for jj = 1 : 50
					binHDO[:,jj] += fit(Histogram,dqdp[jj,:],dqdpbin).weights
				end
				p,dqdp = extract(geoname,"O18_",days)
				for jj = 1 : 50
					binO18[:,jj] += fit(Histogram,dqdp[jj,:],dqdpbin).weights
				end
			end
		else
			geoname = "OTREC_STN$stnstr"
			p,dqdp = extract(geoname,"HDO_",days)
			for jj = 1 : 50
				binHDO[:,jj] += fit(Histogram,dqdp[jj,:],dqdpbin).weights
			end
			p,dqdp = extract(geoname,"O18_",days)
			for jj = 1 : 50
				binO18[:,jj] += fit(Histogram,dqdp[jj,:],dqdpbin).weights
			end
		end
	end
	c1 = axes[2*ii].pcolormesh(binplt,1:50,(binHDO./sum(binHDO,dims=1))'*100,extend="both",levels=5:5:95)
	c2 = axes[2*ii-1].pcolormesh(binplt,1:50,(binO18./sum(binO18,dims=1))'*100,extend="both",levels=5:5:95)

	if cinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ 1343fbae-0ebd-4237-8273-0ebab8325424
function axesformat!(axes)

	for ax in axes
		ax.format(
			xlim=(-10,150)./100,ylim=(1,50),ylabel="Level",
			xlabel=L"$\partial_pq_h$ / 10$^{-7}$ kg kg$^{-1}$ Pa$^{-1}$",
			suptitle="7-Day Moving Average"
		)
	end

	for ii = 1 : 8
		axes[2*ii].format(lrtitle="HDO")
		axes[2*ii-1].format(lrtitle="O18")
	end

	axes[2].format(ultitle="(a) All Stations")
	axes[2].format(ultitle="(a) Colombia")
	axes[4].format(ultitle="(b) San Andres")
	axes[6].format(ultitle="(c) Buenaventura")
	axes[6].text(0.48,39.5,"Bahia Solano",fontsize=10)
	axes[8].format(ultitle="(d) Quibdo")
	axes[10].format(ultitle="(e) Costa Rica")
	axes[12].format(ultitle="(f) EEFMB")
	axes[12].text(0.41,39.5,"ADMQ",fontsize=10)
	axes[12].text(0.41,34.5,"CGFI",fontsize=10)
	axes[14].format(ultitle="(g) Cahuita")
	axes[14].text(0.48,39.5,"Bataan",fontsize=10)
	axes[14].text(0.48,34.5,"Limon",fontsize=10)
	axes[16].format(ultitle="(h) Liberia")
	axes[16].text(0.48,39.5,"OSA",fontsize=10)

	return

end

# ╔═╡ 8c211620-d632-4f23-85f5-a702faf82270
begin
	pplt.close(); fig,axs = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],aspect=0.5,axwidth=0.75,
		wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c1,_ = 
	plotdqdp(axs,1,(-10:5:150)/100,ID=1:4,box=true,bnum=4,days=7,cinfo=true)
	plotdqdp(axs,2,(-10:5:150)/100,ID=1,box=true,bnum=4,days=7)
	plotdqdp(axs,3,(-10:5:150)/100,ID=3:4,box=true,bnum=4,days=7)
	plotdqdp(axs,4,(-10:5:150)/100,ID=2,box=true,bnum=4,days=7)
	plotdqdp(axs,5,(-10:5:150)/100,ID=5:12,box=true,bnum=4,days=7)
	plotdqdp(axs,6,(-10:5:150)/100,ID=5:7,box=true,bnum=4,days=7)
	plotdqdp(axs,7,(-10:5:150)/100,ID=9:11,box=true,bnum=4,days=7)
	plotdqdp(axs,8,(-10:5:150)/100,ID=[8,12],box=true,bnum=4,days=7)
	
	axesformat!(axs)

	fig.colorbar(c1,length=0.75)
	fig.savefig(plotsdir("04d-pvdqdp.png"),transparent=false,dpi=400)
	load(plotsdir("04d-pvdqdp.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─1343fbae-0ebd-4237-8273-0ebab8325424
# ╟─8c211620-d632-4f23-85f5-a702faf82270
