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
		"$geoname-$(iso)QBUDGET-20190801_20201231.nc"
	))
	hvyp = smooth(dropdims(sum(reshape(ds["$(iso)P"][:],24,:),dims=1),dims=1),days)
	hvye = smooth(dropdims(sum(reshape(ds["$(iso)E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf3","processed",
		"$geoname-$(iso)∇decompose-20190801_20201231.nc"
	))
	hvya = smooth(dropdims(mean(reshape(ds["$(iso)ADV"][:],24,:),dims=1),dims=1),days) * 86400
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
		"$geoname-p_wwgt2-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	
	# dsp = NCDataset(datadir(
	# 	"wrf3","processed",
	# 	"$geoname-wcoeff-daily-20190801_20201231-smooth_$(dystr)days.nc"
	# ))
	# wc = dsp["wcoeff"][:,:]
	# close(dsp)
	# wcconv = sum(wc[:,[1,3,5,7,9]],dims=2)[:];
	# wctpbt = sum(wc[:,[2,4,6,8,10]],dims=2)[:];
	# wcabs  = sqrt.(wcconv.^2 .+ wctpbt.^2)[:];
	# wctpbt = wctpbt ./ wcabs
	# ii = wcconv .< 0

	return pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "HDO", prfx = ""
)

	iso = "$(iso)_"
	stnstr = @sprintf("%02d",ID)
	
	geoname = "OTREC_wrf_$(prfx)$stnstr"
	pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg = extract(geoname,iso,days)
	nt = length(prcp)

	for it = 1 : nt
		if (prcp[it]+advc[it]-evap[it]>2.5) && !isnan(pwgt[it])
			rind = argmin(abs.(prcp[it]+advc[it]-evap[it].-rpnt))
			pind = argmin(abs.(pwgt[it].-ppnt))
			numstn[rind,pind] += 1
			binstn[rind,pind] += hvyp[it] + hvya[it] - hvye[it]
			prcstn[rind,pind] += prcp[it] + advc[it] - evap[it]
		end
	end

	return

end

# ╔═╡ 9d38e14e-7226-4d57-ba6f-3b3382dfce1c
function threshold!(
	bin,num,prc;
	numthresh = 15
)

	bin[num.<numthresh] .= 0
	prc[num.<numthresh] .= 0

	return
	
end

# ╔═╡ 6fc8d69a-81d1-47c4-8609-8ec7914bc935
function plotbin!(
	axesnum,ii,
	rainbin,wgtpbin,
	iibin, iiprc, iinum,
	lvls;
	returncinfo = true,
	doalpha = false
)

	tmpbin = sum(iibin,dims=1)
	tmpprc = sum(iiprc,dims=1)
	tmpnum = sum(iinum,dims=1)
	threshold!(tmpbin,tmpnum,tmpprc)
	ix = axesnum[2*ii].panel("l",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,(tmpbin./tmpprc.-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii-1].panel("r",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(xlocator=[])

	tmpbin = sum(iibin,dims=2)
	tmpprc = sum(iiprc,dims=2)
	tmpnum = sum(iinum,dims=2)
	threshold!(tmpbin,tmpnum,tmpprc)
	ix = axesnum[2*ii].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],(tmpbin./tmpprc.-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(ylocator=[])

	ix = axesnum[2*ii-1].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(ylocator=[])

	threshold!(iibin,iinum,iiprc)
	c1 = axesnum[2*ii].pcolormesh(
		rainbin,wgtpbin,(iibin./iiprc .- 1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	c2 = axesnum[2*ii-1].pcolormesh(
		rainbin,wgtpbin,iinum',
		cmap="fire",levels=0:5:50,extend="both"
	)

	if returncinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ 987ab0fb-a376-4470-bcad-0e5681c6ca84
function axesformat!(axesnum)

	for ax in axesnum
		ax.format(
			ylim=(1,0),xlim=(0,120),xlocator=0:100:200,ylocator=0:0.25:1,
			xlabel=L"$P - E + u\cdot\nabla q$ / kg m$^{-2}$ day$^{-1}$",
			leftlabels=["ITCZ","CrossITCZ","PAC2ATL"]
		)
	end

	for ii = 1 : 5
		axesnum[2*ii].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axesnum[2*ii+10].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axesnum[2*ii+20].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
	end
	
	naxs = length(axesnum)
	axesnum[29].format(ylabel=L"$\sigma_\omega$")

	return

end

# ╔═╡ 2fd946e2-bf3e-406f-9a19-5aa72b5d1640
begin
	rbin =  0 : 10 : 250; rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = -1 : 0.05 : 1; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
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
	jbin = zeros(nr,np); jnum = zeros(nr,np); jprc = zeros(nr,np)
	kbin = zeros(nr,np); knum = zeros(nr,np); kprc = zeros(nr,np)
	lbin = zeros(nr,np); lnum = zeros(nr,np); lprc = zeros(nr,np)
	mbin = zeros(nr,np); mnum = zeros(nr,np); mprc = zeros(nr,np)
	nbin = zeros(nr,np); nnum = zeros(nr,np); nprc = zeros(nr,np)
	obin = zeros(nr,np); onum = zeros(nr,np); oprc = zeros(nr,np)
	md"Preallocation of arrays ..."
end

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
	jbin .= 0; jprc .= 0; jnum .= 0
	kbin .= 0; kprc .= 0; knum .= 0
	lbin .= 0; lprc .= 0; lnum .= 0
	mbin .= 0; mprc .= 0; mnum .= 0
	nbin .= 0; nprc .= 0; nnum .= 0
	obin .= 0; oprc .= 0; onum .= 0
	
	for istn = 1 : 5
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=7)
	end
	for istn = 6 : 10
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=7)
	end
	for istn = 11 : 15
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=7)
	end
	for istn = 16 : 20
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=7)
	end
	for istn = 21 : 25
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=7)
	end
	for istn = 1 : 5
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=7)
	end
	for istn = 6 : 10
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=7)
	end
	for istn = 11 : 15
		binning!(hbin,hnum,hprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=7)
	end
	for istn = 16 : 20
		binning!(ibin,inum,iprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=7)
	end
	for istn = 21 : 25
		binning!(jbin,jnum,jprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=7)
	end
	for istn = 1 : 5
		binning!(kbin,knum,kprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=7)
	end
	for istn = 6 : 10
		binning!(lbin,lnum,lprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=7)
	end
	for istn = 11 : 14
		binning!(mbin,mnum,mprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=7)
	end
	for istn = 16 : 20
		binning!(nbin,nnum,nprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=7)
	end
	for istn = 21 : 25
		binning!(obin,onum,oprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=7)
	end
	
	pplt.close(); f1,a1 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],
		aspect=0.6,axwidth=0.6,wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c1_1,c1_2 = 
	plotbin!(a1,1,rbin,pbin,abin,aprc,anum,-150:5:-50,returncinfo=true)
	plotbin!(a1,2,rbin,pbin,bbin,bprc,bnum,-150:5:-50)
	plotbin!(a1,3,rbin,pbin,cbin,cprc,cnum,-150:5:-50)
	plotbin!(a1,4,rbin,pbin,dbin,dprc,dnum,-150:5:-50)
	plotbin!(a1,5,rbin,pbin,ebin,eprc,enum,-150:5:-50)
	plotbin!(a1,6,rbin,pbin,fbin,fprc,fnum,-150:5:-50)
	plotbin!(a1,7,rbin,pbin,gbin,gprc,gnum,-150:5:-50)
	plotbin!(a1,8,rbin,pbin,hbin,hprc,hnum,-150:5:-50)
	plotbin!(a1,9,rbin,pbin,ibin,iprc,inum,-150:5:-50)
	plotbin!(a1,10,rbin,pbin,jbin,jprc,jnum,-150:5:-50)
	plotbin!(a1,11,rbin,pbin,kbin,kprc,knum,-150:5:-50)
	plotbin!(a1,12,rbin,pbin,lbin,lprc,lnum,-150:5:-50)
	plotbin!(a1,13,rbin,pbin,mbin,mprc,mnum,-150:5:-50)
	plotbin!(a1,14,rbin,pbin,nbin,nprc,nnum,-150:5:-50)
	plotbin!(a1,15,rbin,pbin,obin,oprc,onum,-150:5:-50)

	axesformat!(a1)
	a1[1].format(suptitle="7-Day WRF Moving Average")

	f1.colorbar(c1_1,row=[1],locator=-150:50:-50,label=L"$\delta$D / $\perthousand$",minorlocator=-150:5:-45)
	f1.colorbar(c1_2,row=[2],locator=0:10:50,label="Number of Observations")
	
	f1.savefig(plotsdir("05h-p_wwgt2vswrfadvec-07dyHDO.png"),transparent=false,dpi=400)
	load(plotsdir("05h-p_wwgt2vswrfadvec-07dyHDO.png"))
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
	jbin .= 0; jprc .= 0; jnum .= 0
	kbin .= 0; kprc .= 0; knum .= 0
	lbin .= 0; lprc .= 0; lnum .= 0
	mbin .= 0; mprc .= 0; mnum .= 0
	nbin .= 0; nprc .= 0; nnum .= 0
	obin .= 0; oprc .= 0; onum .= 0
	
	for istn = 1 : 5
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 11 : 15
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 1 : 5
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 11 : 15
		binning!(hbin,hnum,hprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(ibin,inum,iprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(jbin,jnum,jprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 1 : 5
		binning!(kbin,knum,kprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(lbin,lnum,lprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 11 : 14
		binning!(mbin,mnum,mprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(nbin,nnum,nprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(obin,onum,oprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	
	pplt.close(); f2,a2 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],
		aspect=0.6,axwidth=0.6,wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c2_1,c2_2 = 
	plotbin!(a2,1,rbin,pbin,abin,aprc,anum,-20:-8,returncinfo=true)
	plotbin!(a2,2,rbin,pbin,bbin,bprc,bnum,-20:-8)
	plotbin!(a2,3,rbin,pbin,cbin,cprc,cnum,-20:-8)
	plotbin!(a2,4,rbin,pbin,dbin,dprc,dnum,-20:-8)
	plotbin!(a2,5,rbin,pbin,ebin,eprc,enum,-20:-8)
	plotbin!(a2,6,rbin,pbin,fbin,fprc,fnum,-20:-8)
	plotbin!(a2,7,rbin,pbin,gbin,gprc,gnum,-20:-8)
	plotbin!(a2,8,rbin,pbin,hbin,hprc,hnum,-20:-8)
	plotbin!(a2,9,rbin,pbin,ibin,iprc,inum,-20:-8)
	plotbin!(a2,10,rbin,pbin,jbin,jprc,jnum,-20:-8)
	plotbin!(a2,11,rbin,pbin,kbin,kprc,knum,-20:-8)
	plotbin!(a2,12,rbin,pbin,lbin,lprc,lnum,-20:-8)
	plotbin!(a2,13,rbin,pbin,mbin,mprc,mnum,-20:-8)
	plotbin!(a2,14,rbin,pbin,nbin,nprc,nnum,-20:-8)
	plotbin!(a2,15,rbin,pbin,obin,oprc,onum,-20:-8)

	axesformat!(a2)
	a2[1].format(suptitle="7-Day WRF Moving Average")

	f2.colorbar(c2_1,row=[1],locator=-20:4:-8,label=L"$\delta^{18}$O / $\perthousand$",minorlocator=-150:5:-45)
	f2.colorbar(c2_2,row=[2],locator=0:10:50,label="Number of Observations")
	
	f2.savefig(plotsdir("05h-p_wwgt2vswrfadvec-07dyO18.png"),transparent=false,dpi=400)
	load(plotsdir("05h-p_wwgt2vswrfadvec-07dyO18.png"))
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
	jbin .= 0; jprc .= 0; jnum .= 0
	kbin .= 0; kprc .= 0; knum .= 0
	lbin .= 0; lprc .= 0; lnum .= 0
	mbin .= 0; mprc .= 0; mnum .= 0
	nbin .= 0; nprc .= 0; nnum .= 0
	obin .= 0; oprc .= 0; onum .= 0
	
	for istn = 1 : 5
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=30)
	end
	for istn = 6 : 10
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=30)
	end
	for istn = 11 : 15
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=30)
	end
	for istn = 16 : 20
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=30)
	end
	for istn = 21 : 25
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="HDO",days=30)
	end
	for istn = 1 : 5
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=30)
	end
	for istn = 6 : 10
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=30)
	end
	for istn = 11 : 15
		binning!(hbin,hnum,hprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=30)
	end
	for istn = 16 : 20
		binning!(ibin,inum,iprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=30)
	end
	for istn = 21 : 25
		binning!(jbin,jnum,jprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="HDO",days=30)
	end
	for istn = 1 : 5
		binning!(kbin,knum,kprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=30)
	end
	for istn = 6 : 10
		binning!(lbin,lnum,lprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=30)
	end
	for istn = 11 : 14
		binning!(mbin,mnum,mprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=30)
	end
	for istn = 16 : 20
		binning!(nbin,nnum,nprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=30)
	end
	for istn = 21 : 25
		binning!(obin,onum,oprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="HDO",days=30)
	end
	
	pplt.close(); f3,a3 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],
		aspect=0.6,axwidth=0.6,wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c3_1,c3_2 = 
	plotbin!(a3,1,rbin,pbin,abin,aprc,anum,-150:5:-50,returncinfo=true)
	plotbin!(a3,2,rbin,pbin,bbin,bprc,bnum,-150:5:-50)
	plotbin!(a3,3,rbin,pbin,cbin,cprc,cnum,-150:5:-50)
	plotbin!(a3,4,rbin,pbin,dbin,dprc,dnum,-150:5:-50)
	plotbin!(a3,5,rbin,pbin,ebin,eprc,enum,-150:5:-50)
	plotbin!(a3,6,rbin,pbin,fbin,fprc,fnum,-150:5:-50)
	plotbin!(a3,7,rbin,pbin,gbin,gprc,gnum,-150:5:-50)
	plotbin!(a3,8,rbin,pbin,hbin,hprc,hnum,-150:5:-50)
	plotbin!(a3,9,rbin,pbin,ibin,iprc,inum,-150:5:-50)
	plotbin!(a3,10,rbin,pbin,jbin,jprc,jnum,-150:5:-50)
	plotbin!(a3,11,rbin,pbin,kbin,kprc,knum,-150:5:-50)
	plotbin!(a3,12,rbin,pbin,lbin,lprc,lnum,-150:5:-50)
	plotbin!(a3,13,rbin,pbin,mbin,mprc,mnum,-150:5:-50)
	plotbin!(a3,14,rbin,pbin,nbin,nprc,nnum,-150:5:-50)
	plotbin!(a3,15,rbin,pbin,obin,oprc,onum,-150:5:-50)

	axesformat!(a3)
	a3[1].format(suptitle="30-Day WRF Moving Average")

	f3.colorbar(c3_1,row=[1],locator=-150:50:-50,label=L"$\delta$D / $\perthousand$",minorlocator=-150:5:-45)
	f3.colorbar(c3_2,row=[2],locator=0:10:50,label="Number of Observations")
	
	f3.savefig(plotsdir("05h-p_wwgt2vswrfadvec-30dyHDO.png"),transparent=false,dpi=400)
	load(plotsdir("05h-p_wwgt2vswrfadvec-30dyHDO.png"))
end

# ╔═╡ e2d4c2ec-8855-401e-98ba-84124a549c9a
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
	jbin .= 0; jprc .= 0; jnum .= 0
	kbin .= 0; kprc .= 0; knum .= 0
	lbin .= 0; lprc .= 0; lnum .= 0
	mbin .= 0; mprc .= 0; mnum .= 0
	nbin .= 0; nprc .= 0; nnum .= 0
	obin .= 0; oprc .= 0; onum .= 0
	
	for istn = 1 : 5
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=30)
	end
	for istn = 6 : 10
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=30)
	end
	for istn = 11 : 15
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=30)
	end
	for istn = 16 : 20
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=30)
	end
	for istn = 21 : 25
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=30)
	end
	for istn = 1 : 5
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=30)
	end
	for istn = 6 : 10
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=30)
	end
	for istn = 11 : 15
		binning!(hbin,hnum,hprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=30)
	end
	for istn = 16 : 20
		binning!(ibin,inum,iprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=30)
	end
	for istn = 21 : 25
		binning!(jbin,jnum,jprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=30)
	end
	for istn = 1 : 5
		binning!(kbin,knum,kprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=30)
	end
	for istn = 6 : 10
		binning!(lbin,lnum,lprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=30)
	end
	for istn = 11 : 14
		binning!(mbin,mnum,mprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=30)
	end
	for istn = 16 : 20
		binning!(nbin,nnum,nprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=30)
	end
	for istn = 21 : 25
		binning!(obin,onum,oprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=30)
	end
	
	pplt.close(); f4,a4 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],
		aspect=0.6,axwidth=0.6,wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c4_1,c4_2 = 
	plotbin!(a4,1,rbin,pbin,abin,aprc,anum,-20:-8,returncinfo=true)
	plotbin!(a4,2,rbin,pbin,bbin,bprc,bnum,-20:-8)
	plotbin!(a4,3,rbin,pbin,cbin,cprc,cnum,-20:-8)
	plotbin!(a4,4,rbin,pbin,dbin,dprc,dnum,-20:-8)
	plotbin!(a4,5,rbin,pbin,ebin,eprc,enum,-20:-8)
	plotbin!(a4,6,rbin,pbin,fbin,fprc,fnum,-20:-8)
	plotbin!(a4,7,rbin,pbin,gbin,gprc,gnum,-20:-8)
	plotbin!(a4,8,rbin,pbin,hbin,hprc,hnum,-20:-8)
	plotbin!(a4,9,rbin,pbin,ibin,iprc,inum,-20:-8)
	plotbin!(a4,10,rbin,pbin,jbin,jprc,jnum,-20:-8)
	plotbin!(a4,11,rbin,pbin,kbin,kprc,knum,-20:-8)
	plotbin!(a4,12,rbin,pbin,lbin,lprc,lnum,-20:-8)
	plotbin!(a4,13,rbin,pbin,mbin,mprc,mnum,-20:-8)
	plotbin!(a4,14,rbin,pbin,nbin,nprc,nnum,-20:-8)
	plotbin!(a4,15,rbin,pbin,obin,oprc,onum,-20:-8)

	axesformat!(a4)
	a4[1].format(suptitle="30-Day WRF Moving Average")

	f4.colorbar(c4_1,row=[1],locator=-20:4:-8,label=L"$\delta^{18}$O / $\perthousand$")
	f4.colorbar(c4_2,row=[2],locator=0:10:50,label="Number of Observations")
	
	f4.savefig(plotsdir("05h-p_wwgt2vswrfadvec-30dyO18.png"),transparent=false,dpi=400)
	load(plotsdir("05h-p_wwgt2vswrfadvec-30dyO18.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╠═bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╠═441f47a7-5757-4b24-8b52-a2877e0f0287
# ╠═f1720645-69a8-4f45-a6c1-8c06279d3590
# ╠═4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─9d38e14e-7226-4d57-ba6f-3b3382dfce1c
# ╟─6fc8d69a-81d1-47c4-8609-8ec7914bc935
# ╟─987ab0fb-a376-4470-bcad-0e5681c6ca84
# ╟─2fd946e2-bf3e-406f-9a19-5aa72b5d1640
# ╟─b6500812-fd5e-4842-8855-655822d170f4
# ╟─c57ae725-3056-481c-a004-a916192744be
# ╟─738c6dde-cc6b-489f-8251-849e4ca67d8c
# ╟─e2d4c2ec-8855-401e-98ba-84124a549c9a
