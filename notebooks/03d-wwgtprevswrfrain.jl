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
	using Statistics
	
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

# ╔═╡ 441f47a7-5757-4b24-8b52-a2877e0f0287
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

# ╔═╡ f1720645-69a8-4f45-a6c1-8c06279d3590
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

# ╔═╡ 4319fd0e-fd9f-424e-9286-3b3b5a844b73
function extract(geoname,iso,days)

	dystr  = @sprintf("%02d",days)

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-rain-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	prcp = ds["RAINNC"][:]
	hvyp = ds["$(iso)RAINNC"][:]
	close(ds)

	dsp = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-wcoeff-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	wc = dsp["wcoeff"][:,:]
	close(dsp)
	wcconv = sum(wc[:,[1,3,5,7,9]],dims=2)[:];
	wctpbt = sum(wc[:,[2,4,6,8,10]],dims=2)[:];
	wcabs  = sqrt.(wcconv.^2 .+ wctpbt.^2)[:];
	wctpbt = wctpbt ./ wcabs
	ii = wcconv .< 0

	dsp = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-p_wwgt-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	pw = dsp["p_wwgt"][:] ./ 100
	close(dsp)

	dsp = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-p_wwgt-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	pw = dsp["σ_wwgt"][:]
	close(dsp)

	return pw[ii],prcp[ii],hvyp[ii]
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "HDO", box = false, bnum = 1
)

	iso = "$(iso)_"
	stnstr = @sprintf("%02d",ID)

	if box
		for ibox = 1 : bnum
			boxstr = @sprintf("%d",ibox)
			geoname = "OTREC_wrf_stn$(stnstr)b_box$(boxstr)"
			pwgt,prcp,hvyp = extract(geoname,iso,days)
			nt = length(prcp)
		
			for it = 1 : nt
				if (prcp[it]>2.5) && !isnan(pwgt[it])
					rind = argmin(abs.(prcp[it].-rpnt))
					pind = argmin(abs.(pwgt[it].-ppnt))
					numstn[rind,pind] += 1
					binstn[rind,pind] += hvyp[it]
					prcstn[rind,pind] += prcp[it]
				end
			end
		end
	else
		geoname = "OTREC_wrf_stn$(stnstr)b"
		pwgt,prcp,hvyp = extract(geoname,iso,days)
		nt = length(prcp)
	
		for it = 1 : nt
			if (prcp[it]>2.5) && !isnan(pwgt[it])
				rind = argmin(abs.(prcp[it].-rpnt))
				pind = argmin(abs.(pwgt[it].-ppnt))
				numstn[rind,pind] += 1
				binstn[rind,pind] += hvyp[it]
				prcstn[rind,pind] += prcp[it]
			end
		end
	end

	return

end

# ╔═╡ 42c325e3-d1df-4e09-8d59-8e1505210c43
function threshold!(
	bin,num,prc;
	numthresh = 5
)

	bin[num.<numthresh] .= 0
	prc[num.<numthresh] .= 0
	# num[num.<numthresh] .= 0

	return
	
end

# ╔═╡ 7e66a747-056c-445e-a021-0d919ddc26bb
function plotbin!(
	axesnum,ii,
	rainbin,wgtpbin,
	iibin, iiprc, iinum,
	lvls;
	returncinfo = true
)

	c1 = axesnum[2*ii].pcolormesh(
		rainbin,wgtpbin,(iibin./iiprc .- 1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	c2 = axesnum[2*ii-1].pcolormesh(
		rainbin,wgtpbin,iinum' ./ sum(iinum)*100,
		cmap="fire",levels=0:10,extend="both"
	)

	ix = axesnum[2*ii].panel("l",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,(sum(iibin,dims=1)./sum(iiprc,dims=1).-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],(sum(iibin,dims=2)./sum(iiprc,dims=2).-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(ylocator=[],xlabel=L"$P$ / kg m$^{-2}$ day$^{-1}$")

	ix = axesnum[2*ii-1].panel("r",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,sum(iinum,dims=1)' ./ sum(iinum)*100,
		cmap="fire",levels=0:10,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii-1].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],sum(iinum,dims=2)' ./ sum(iinum)*100,
		cmap="fire",levels=0:10,extend="both"
	)
	ix.format(ylocator=[],xlabel=L"$P$ / kg m$^{-2}$ day$^{-1}$")

	if returncinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ 3bb9d01b-b214-44b1-975e-fcab56d8eb99
function axesformat!(axesnum)

	for ax in axesnum
		ax.format(
			ylim=(1,0),xlim=(0,60),xlocator=0:25:100,ylabel=L"$p_\omega$ / hPa",
			xlabel=L"$P$ / kg m$^{-2}$ day$^{-1}$"
		)
	end

	axesnum[2].format(ultitle="(a) All Stations")
	axesnum[2].format(ultitle="(a) Colombia")
	axesnum[4].format(ultitle="(b) San Andres")
	axesnum[6].format(ultitle="(c) Buenaventura")
	axesnum[6].text(11,0.220,"Bahia Solano",fontsize=10)
	axesnum[8].format(ultitle="(d) Quibdo")
	axesnum[10].format(ultitle="(e) Costa Rica")
	axesnum[12].format(ultitle="(f) EEFMB")
	axesnum[12].text(10,0.220,"ADMQ",fontsize=10)
	axesnum[12].text(10,0.320,"CGFI",fontsize=10)
	axesnum[14].format(ultitle="(g) Cahuita")
	axesnum[14].text(12,0.220,"Bataan",fontsize=10)
	axesnum[14].text(12,0.320,"Limon",fontsize=10)
	axesnum[16].format(ultitle="(h) Liberia")
	axesnum[16].text(10.5,0.220,"OSA",fontsize=10)

	return

end

# ╔═╡ 2fd946e2-bf3e-406f-9a19-5aa72b5d1640
begin
	rbin = 0 : 5 : 200; rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 0 : 0.05 : 1; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
	nr = length(rpnt); np = length(ppnt)
	abin = zeros(nr,np); anum = zeros(nr,np); aprc = zeros(nr,np)
	bbin = zeros(nr,np); bnum = zeros(nr,np); bprc = zeros(nr,np)
	cbin = zeros(nr,np); cnum = zeros(nr,np); cprc = zeros(nr,np)
	dbin = zeros(nr,np); dnum = zeros(nr,np); dprc = zeros(nr,np)
	ebin = zeros(nr,np); enum = zeros(nr,np); eprc = zeros(nr,np)
	fbin = zeros(nr,np); fnum = zeros(nr,np); fprc = zeros(nr,np)
	gbin = zeros(nr,np); gnum = zeros(nr,np); gprc = zeros(nr,np)
	hbin = zeros(nr,np); hnum = zeros(nr,np); hprc = zeros(nr,np)
	ibin = zeros(nr,np); inum = zeros(nr,np); iprc = zeros(nr,np)
	md"Preallocation of arrays ..."
end

# ╔═╡ d139bbdb-2e72-4259-a85d-ecb7addcdd45


# ╔═╡ b6500812-fd5e-4842-8855-655822d170f4
# ╠═╡ show_logs = false
begin
	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0
	gbin .= 0; gprc .= 0; gnum .= 0
	hbin .= 0; hprc .= 0; hnum .= 0
	ibin .= 0; iprc .= 0; inum .= 0
	
	for istn = 1
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end
	for istn = [3,4]
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end
	for istn = 2
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end
	for istn = 5 : 7
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end
	for istn = 9 : 11
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end
	for istn = [8,12]
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=7)
	end

	hbin .= bbin.+cbin.+dbin; ibin .= ebin.+fbin.+gbin
	hprc .= bprc.+cprc.+dprc; iprc .= eprc.+fprc.+gprc
	hnum .= bnum.+cnum.+dnum; inum .= enum.+fnum.+gnum

	threshold!(abin,anum,aprc)
	threshold!(bbin,bnum,bprc)
	threshold!(cbin,cnum,cprc)
	threshold!(dbin,dnum,dprc)
	threshold!(ebin,enum,eprc)
	threshold!(fbin,fnum,fprc)
	threshold!(gbin,gnum,gprc)
	threshold!(hbin,hnum,hprc)
	threshold!(ibin,inum,iprc)
	
	pplt.close(); f1,a1 = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],
		aspect=0.5,axwidth=0.75,wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c1_1,c1_2 = 
	plotbin!(a1,1,rbin,pbin,hbin,hprc,hnum,-75:2.5:-30,returncinfo=true)
	plotbin!(a1,2,rbin,pbin,bbin,bprc,bnum,-75:2.5:-30)
	plotbin!(a1,3,rbin,pbin,cbin,cprc,cnum,-75:2.5:-30)
	plotbin!(a1,4,rbin,pbin,dbin,dprc,dnum,-75:2.5:-30)
	plotbin!(a1,5,rbin,pbin,ibin,iprc,inum,-75:2.5:-30)
	plotbin!(a1,6,rbin,pbin,ebin,eprc,enum,-75:2.5:-30)
	plotbin!(a1,7,rbin,pbin,fbin,fprc,fnum,-75:2.5:-30)
	plotbin!(a1,8,rbin,pbin,gbin,gprc,gnum,-75:2.5:-30)

	axesformat!(a1)
	a1[1].format(suptitle="7-Day WRF Moving Average")

	f1.colorbar(c1_1,loc="r",rows=1,locator=-120:15:-30,label=L"$\delta^2$H / $\perthousand$")
	f1.colorbar(c1_2,loc="r",rows=2,locator=0:2:10,label="Probability Density / %")
	
	f1.savefig(plotsdir("03d-p_wwgtvswrfrain-07dyHDO.png"),transparent=false,dpi=400)
	load(plotsdir("03d-p_wwgtvswrfrain-07dyHDO.png"))
end

# ╔═╡ c57ae725-3056-481c-a004-a916192744be
# ╠═╡ show_logs = false
begin
	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0
	gbin .= 0; gprc .= 0; gnum .= 0
	hbin .= 0; hprc .= 0; hnum .= 0
	ibin .= 0; iprc .= 0; inum .= 0
	
	for istn = 1
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end
	for istn = [3,4]
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end
	for istn = 2
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end
	for istn = 5 : 7
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end
	for istn = 9 : 11
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end
	for istn = [8,12]
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=7)
	end

	hbin .= bbin.+cbin.+dbin; ibin .= ebin.+fbin.+gbin
	hprc .= bprc.+cprc.+dprc; iprc .= eprc.+fprc.+gprc
	hnum .= bnum.+cnum.+dnum; inum .= enum.+fnum.+gnum

	threshold!(abin,anum,aprc)
	threshold!(bbin,bnum,bprc)
	threshold!(cbin,cnum,cprc)
	threshold!(dbin,dnum,dprc)
	threshold!(ebin,enum,eprc)
	threshold!(fbin,fnum,fprc)
	threshold!(gbin,gnum,gprc)
	threshold!(hbin,hnum,hprc)
	threshold!(ibin,inum,iprc)

	pplt.close(); f2,a2 = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],
		aspect=0.5,axwidth=0.75,wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c2_1,c2_2 = 
	plotbin!(a2,1,rbin,pbin,hbin,hprc,hnum,-13:0.5:-6,returncinfo=true)
	plotbin!(a2,2,rbin,pbin,bbin,bprc,bnum,-13:0.5:-6)
	plotbin!(a2,3,rbin,pbin,cbin,cprc,cnum,-13:0.5:-6)
	plotbin!(a2,4,rbin,pbin,dbin,dprc,dnum,-13:0.5:-6)
	plotbin!(a2,5,rbin,pbin,ibin,iprc,inum,-13:0.5:-6)
	plotbin!(a2,6,rbin,pbin,ebin,eprc,enum,-13:0.5:-6)
	plotbin!(a2,7,rbin,pbin,fbin,fprc,fnum,-13:0.5:-6)
	plotbin!(a2,8,rbin,pbin,gbin,gprc,gnum,-13:0.5:-6)

	axesformat!(a2)
	a2[1].format(suptitle="7-Day WRF Moving Average")
	
	f2.colorbar(c2_1,loc="r",rows=1,locator=-15:-0,label=L"$\delta^{18}$O / $\perthousand$")
	f2.colorbar(c2_2,loc="r",rows=2,locator=0:2:10,label="Probability Density / %")
	
	f2.savefig(plotsdir("03d-p_wwgtvswrfrain-07dyO18.png"),transparent=false,dpi=400)
	load(plotsdir("03d-p_wwgtvswrfrain-07dyO18.png"))
end

# ╔═╡ 738c6dde-cc6b-489f-8251-849e4ca67d8c
# ╠═╡ show_logs = false
begin
	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0
	gbin .= 0; gprc .= 0; gnum .= 0
	hbin .= 0; hprc .= 0; hnum .= 0
	ibin .= 0; iprc .= 0; inum .= 0
	
	for istn = 1
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end
	for istn = [3,4]
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end
	for istn = 2
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end
	for istn = 5 : 7
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end
	for istn = 9 : 11
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end
	for istn = [8,12]
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,box=false,bnum=4,days=30)
	end

	hbin .= bbin.+cbin.+dbin; ibin .= ebin.+fbin.+gbin
	hprc .= bprc.+cprc.+dprc; iprc .= eprc.+fprc.+gprc
	hnum .= bnum.+cnum.+dnum; inum .= enum.+fnum.+gnum

	threshold!(abin,anum,aprc)
	threshold!(bbin,bnum,bprc)
	threshold!(cbin,cnum,cprc)
	threshold!(dbin,dnum,dprc)
	threshold!(ebin,enum,eprc)
	threshold!(fbin,fnum,fprc)
	threshold!(gbin,gnum,gprc)
	threshold!(hbin,hnum,hprc)
	threshold!(ibin,inum,iprc)
	
	pplt.close(); f3,a3 = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],
		aspect=0.5,axwidth=0.75,wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c3_1,c3_2 = 
	plotbin!(a3,1,rbin,pbin,hbin,hprc,hnum,-75:2.5:-30,returncinfo=true)
	plotbin!(a3,2,rbin,pbin,bbin,bprc,bnum,-75:2.5:-30)
	plotbin!(a3,3,rbin,pbin,cbin,cprc,cnum,-75:2.5:-30)
	plotbin!(a3,4,rbin,pbin,dbin,dprc,dnum,-75:2.5:-30)
	plotbin!(a3,5,rbin,pbin,ibin,iprc,inum,-75:2.5:-30)
	plotbin!(a3,6,rbin,pbin,ebin,eprc,enum,-75:2.5:-30)
	plotbin!(a3,7,rbin,pbin,fbin,fprc,fnum,-75:2.5:-30)
	plotbin!(a3,8,rbin,pbin,gbin,gprc,gnum,-75:2.5:-30)

	axesformat!(a3)
	a3[1].format(suptitle="30-Day WRF Moving Average")

	f3.colorbar(c3_1,loc="r",rows=1,locator=-150:15:-40,label=L"$\delta^2$H / $\perthousand$")
	f3.colorbar(c3_2,loc="r",rows=2,locator=0:2:10,label="Probability Density / %")
	
	f3.savefig(plotsdir("03d-p_wwgtvswrfrain-30dyHDO.png"),transparent=false,dpi=400)
	load(plotsdir("03d-p_wwgtvswrfrain-30dyHDO.png"))
end

# ╔═╡ 3d3a2954-fe7e-40da-a6d9-fe9affcc581b
# ╠═╡ show_logs = false
begin
	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0
	gbin .= 0; gprc .= 0; gnum .= 0
	hbin .= 0; hprc .= 0; hnum .= 0
	ibin .= 0; iprc .= 0; inum .= 0
	
	for istn = 1
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end
	for istn = [3,4]
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end
	for istn = 2
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end
	for istn = 5 : 7
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end
	for istn = 9 : 11
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end
	for istn = [8,12]
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,iso="O18",box=false,bnum=4,days=30)
	end

	hbin .= bbin.+cbin.+dbin; ibin .= ebin.+fbin.+gbin
	hprc .= bprc.+cprc.+dprc; iprc .= eprc.+fprc.+gprc
	hnum .= bnum.+cnum.+dnum; inum .= enum.+fnum.+gnum

	threshold!(abin,anum,aprc)
	threshold!(bbin,bnum,bprc)
	threshold!(cbin,cnum,cprc)
	threshold!(dbin,dnum,dprc)
	threshold!(ebin,enum,eprc)
	threshold!(fbin,fnum,fprc)
	threshold!(gbin,gnum,gprc)
	threshold!(hbin,hnum,hprc)
	threshold!(ibin,inum,iprc)

	pplt.close(); f4,a4 = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],
		aspect=0.5,axwidth=0.75,wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c4_1,c4_2 = 
	plotbin!(a4,1,rbin,pbin,hbin,hprc,hnum,-12:0.5:-6,returncinfo=true)
	plotbin!(a4,2,rbin,pbin,bbin,bprc,bnum,-12:0.5:-6)
	plotbin!(a4,3,rbin,pbin,cbin,cprc,cnum,-12:0.5:-6)
	plotbin!(a4,4,rbin,pbin,dbin,dprc,dnum,-12:0.5:-6)
	plotbin!(a4,5,rbin,pbin,ibin,iprc,inum,-12:0.5:-6)
	plotbin!(a4,6,rbin,pbin,ebin,eprc,enum,-12:0.5:-6)
	plotbin!(a4,7,rbin,pbin,fbin,fprc,fnum,-12:0.5:-6)
	plotbin!(a4,8,rbin,pbin,gbin,gprc,gnum,-12:0.5:-6)

	axesformat!(a4)
	a4[1].format(suptitle="30-Day WRF Moving Average")
	
	f4.colorbar(c4_1,loc="r",rows=1,locator=-13:-7,label=L"$\delta^{18}$O / $\perthousand$")
	f4.colorbar(c4_2,loc="r",rows=2,locator=0:2:10,label="Probability Density / %")
	
	f4.savefig(plotsdir("03d-p_wwgtvswrfrain-30dyO18.png"),transparent=false,dpi=400)
	load(plotsdir("03d-p_wwgtvswrfrain-30dyO18.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─441f47a7-5757-4b24-8b52-a2877e0f0287
# ╟─f1720645-69a8-4f45-a6c1-8c06279d3590
# ╠═4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╠═5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─42c325e3-d1df-4e09-8d59-8e1505210c43
# ╟─7e66a747-056c-445e-a021-0d919ddc26bb
# ╠═3bb9d01b-b214-44b1-975e-fcab56d8eb99
# ╠═2fd946e2-bf3e-406f-9a19-5aa72b5d1640
# ╠═d139bbdb-2e72-4259-a85d-ecb7addcdd45
# ╠═b6500812-fd5e-4842-8855-655822d170f4
# ╠═c57ae725-3056-481c-a004-a916192744be
# ╠═738c6dde-cc6b-489f-8251-849e4ca67d8c
# ╠═3d3a2954-fe7e-40da-a6d9-fe9affcc581b
