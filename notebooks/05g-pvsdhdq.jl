### A Pluto.jl notebook ###
# v0.20.10

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
	pplt = pyimport("ultraplot")

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

# ╔═╡ 87415583-395d-4a39-b3d3-3f6a2dc95112
function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

# ╔═╡ b754ea7c-f907-47ba-87d3-c0df3a92dc3d
function smooth(data::AbstractVector,days)

	buffer,weights = calculatebufferweights(days)

	ndt = length(data)
	ndata = fill(NaN,ndt)
	smth  = zeros(1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		for ismth = 0 : (buffer*2)
			smth[ismth+1] = data[ii+ismth-buffer] * weights[ismth+1]
		end
		ndata[ii] = sum(smth)

	end

	return ndata

end

# ╔═╡ 7e114c52-4f6b-45e1-8081-abbc674d39c2
function extractbudget(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-∇decompose-20190801_20201231.nc"
	))
	advc = smooth(dropdims(mean(reshape(ds["ADV"][:],24,:),dims=1),dims=1),days) * 86400
	divg = smooth(dropdims(mean(reshape(ds["DIV"][:],24,:),dims=1),dims=1),days) * 86400
	close(ds)

	dsp = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-p_wwgt-$dystr.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	return pwgt,prcp,evap,advc,divg
	
end

# ╔═╡ 4319fd0e-fd9f-424e-9286-3b3b5a844b73
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf3","processed","$geoname-dhodq-$(dystr).nc"))
	p = ds["P"][:,:] ./ 100
	dhdq = ds["dHDOdH2O"][:,:]
	dodq = ds["dO18dH2O"][:,:]
	close(ds)

	return p,dhdq,dodq
	
end

# ╔═╡ 43e17868-ad3e-4e14-a925-2414b586558d
prfx = "CrossITCZ"

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function plotdqdp(
	axes,ii;
	ID, days=0, prfx = "", cinfo = false
)

	dhdqbin = 0 : 0.01 : 1.2
	dodqbin = 0.9 : 0.001 : 1.02
	pbin    = vcat(0 : 50 : 1000)
	binHDO = zeros(length(dhdqbin)-1,length(pbin)-1)
	binO18 = zeros(length(dodqbin)-1,length(pbin)-1)
	dhdqbinplt = (dhdqbin[1:(end-1)] .+ dhdqbin[2:end])/2
	dodqbinplt = (dodqbin[1:(end-1)] .+ dodqbin[2:end])/2
	pbinplt = (pbin[1:(end-1)] .+ pbin[2:end])/2

	for stn in ID
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_$(prfx)$stnstr"
		p,dhdq,dodq = extract(geoname,days)
		pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt))
		for jj = 1 : 49
			binHDO[:,:] += fit(
				Histogram,(dhdq[jj,it],p[jj,it]),
				(dhdqbin,pbin)
			).weights
			binO18[:,:] += fit(
				Histogram,(dodq[jj,it],p[jj,it]),
				(dodqbin,pbin)
			).weights
		end
	end

	binHDO ./= sum(binHDO,dims=1) / 100
	binO18 ./= sum(binO18,dims=1) / 100

	lvls = vcat(0.1,0.5,1,2:2:10,15,20)
	c1 = axes[2*ii].pcolormesh(dhdqbinplt,pbinplt,binHDO',extend="both",levels=lvls)
	c2 = axes[2*ii-1].pcolormesh(dodqbinplt,pbinplt,binO18',extend="both",levels=lvls)

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
			xlim=(-0.2,1.2),ylim=(1000,100),ylabel="Pressure / hPa",
			xlabel=L"$\partial_p(q_h/q)$ / $\perthousand$ 10$^{-1}$ hPa$^{-1}$",
			suptitle="7-Day Moving Average",
			rightlabels=["ITCZ","CrossITCZ","PAC2ATL"]
		)
	end

	for ii = 1 : 15
		axes[2*ii].format(xlim=(0.5,1.1),lltitle="HDO")
		axes[2*ii-1].format(xlim=(0.92,1.01),lltitle="O18")
	end

	for ii = 1 : 5
		axes[2*ii].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axes[2*ii+10].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axes[2*ii+20].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
	end

	return

end

# ╔═╡ 8c211620-d632-4f23-85f5-a702faf82270
begin
	pplt.close(); fig,axs = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],aspect=0.5,axwidth=0.5,
		wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c1,_ = 
	plotdqdp(axs,1,ID=1:5,  prfx="ITCZ",days=7,cinfo=true)
	plotdqdp(axs,2,ID=6:10, prfx="ITCZ",days=7)
	plotdqdp(axs,3,ID=11:15,prfx="ITCZ",days=7)
	plotdqdp(axs,4,ID=16:20,prfx="ITCZ",days=7)
	plotdqdp(axs,5,ID=20:25,prfx="ITCZ",days=7)

	plotdqdp(axs,6,ID=1:5,  prfx="CrossITCZ",days=7)
	plotdqdp(axs,7,ID=6:10, prfx="CrossITCZ",days=7)
	plotdqdp(axs,8,ID=11:15,prfx="CrossITCZ",days=7)
	plotdqdp(axs,9,ID=16:20,prfx="CrossITCZ",days=7)
	plotdqdp(axs,10,ID=20:25,prfx="CrossITCZ",days=7)

	plotdqdp(axs,11,ID=1:5,  prfx="PAC2ATL",days=7)
	plotdqdp(axs,12,ID=6:10, prfx="PAC2ATL",days=7)
	plotdqdp(axs,13,ID=11:15,prfx="PAC2ATL",days=7)
	plotdqdp(axs,14,ID=16:20,prfx="PAC2ATL",days=7)
	plotdqdp(axs,15,ID=20:25,prfx="PAC2ATL",days=7)
	
	axesformat!(axs)

	fig.colorbar(c1,loc="b",length=0.5,label="Number of Observations")
	fig.savefig(plotsdir("05g-pvdhqdp-smooth07.png"),transparent=false,dpi=400)
	load(plotsdir("05g-pvdhqdp-smooth07.png"))
end

# ╔═╡ c6660657-3b34-43a7-aaeb-67ceb66a89ec
begin
	pplt.close(); f2,a2 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],aspect=0.5,axwidth=0.5,
		wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c2,_ = 
	plotdqdp(a2,1,ID=1:5,  prfx="ITCZ",days=30,cinfo=true)
	plotdqdp(a2,2,ID=6:10, prfx="ITCZ",days=30)
	plotdqdp(a2,3,ID=11:15,prfx="ITCZ",days=30)
	plotdqdp(a2,4,ID=16:20,prfx="ITCZ",days=30)
	plotdqdp(a2,5,ID=20:25,prfx="ITCZ",days=30)

	plotdqdp(a2,6,ID=1:5,  prfx="CrossITCZ",days=30)
	plotdqdp(a2,7,ID=6:10, prfx="CrossITCZ",days=30)
	plotdqdp(a2,8,ID=11:15,prfx="CrossITCZ",days=30)
	plotdqdp(a2,9,ID=16:20,prfx="CrossITCZ",days=30)
	plotdqdp(a2,10,ID=20:25,prfx="CrossITCZ",days=30)

	plotdqdp(a2,11,ID=1:5,  prfx="PAC2ATL",days=30)
	plotdqdp(a2,12,ID=6:10, prfx="PAC2ATL",days=30)
	plotdqdp(a2,13,ID=11:15,prfx="PAC2ATL",days=30)
	plotdqdp(a2,14,ID=16:20,prfx="PAC2ATL",days=30)
	plotdqdp(a2,15,ID=20:25,prfx="PAC2ATL",days=30)
	
	axesformat!(a2)

	f2.colorbar(c2,loc="b",length=0.5,label="Number of Observations")
	f2.savefig(plotsdir("05g-pvdhqdp-smooth30.png"),transparent=false,dpi=400)
	load(plotsdir("05g-pvdhqdp-smooth30.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─87415583-395d-4a39-b3d3-3f6a2dc95112
# ╟─b754ea7c-f907-47ba-87d3-c0df3a92dc3d
# ╟─7e114c52-4f6b-45e1-8081-abbc674d39c2
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╠═43e17868-ad3e-4e14-a925-2414b586558d
# ╠═5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╠═1343fbae-0ebd-4237-8273-0ebab8325424
# ╟─8c211620-d632-4f23-85f5-a702faf82270
# ╟─c6660657-3b34-43a7-aaeb-67ceb66a89ec
