### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 3c4b037f-c854-4148-ae6e-fff50bbb51cf
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
begin
	@quickactivate "ColombiaIsotope"
	using Dates
	using DelimitedFiles
	using Dierckx
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ColombiaIsotope project..."
end

# ╔═╡ c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
md"
# 04g. Isotope vs W-Weighted Pressure LandSea Mask
"

# ╔═╡ 1ae21e79-5b6d-4f31-8c62-9ed7d3ac6e79
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	clon = coast[:,1]
	clat = coast[:,2]
	md"Loading global coastlines ..."
end

# ╔═╡ 22dfb33e-276a-4285-81c8-5f3300c50c77
TableOfContents()

# ╔═╡ c853d30e-5a51-47f4-a6be-88ab857253c1
md"
### A. Loading the data from the Pacific Coast
"

# ╔═╡ 18e518f8-df88-4e57-9123-7852546af717
geo = [GeoRegion("OTREC_PCS"),GeoRegion("OTREC_PAC"),GeoRegion("OTREC_ATL")]

# ╔═╡ 2489ca4e-0765-40b0-b131-49d8a1563727
begin
	rainb = 5  : 5 : 95; nr = length(rainb)
	rainu = 10 : 5 : 100
	rainm = 5  : 5 : 100
	pwgtb = 100 : 25 : 925; np = length(pwgtb)
	pwgtu = 125 : 25 : 950
	pwgtm = 100 : 25 : 950'
	md"Defining ranges for binning of isotopes and rainfall ..."
end

# ╔═╡ 6b55be2c-2f4e-4c41-8eb6-24e6d769df00
begin
	lsm   = Vector{Array}(undef,1)
	ginfo = Vector{RegionGrid}(undef,3)
	for ii = 1 : 3
		ds  = NCDataset(datadir("flsm","flsm_wrf.nc"))
		wln = ds["longitude"][:]
		wlt = ds["latitude"][:]
		lsm[1] = ds["flsm"][:]
		close(ds)
		ginfo[ii] = RegionGrid(geo[ii],wln,wlt)
	end
end

# ╔═╡ c26cdc06-70f1-444a-b190-e2beb4b00ef7
lsmvec = [0,1e-3,0.1,0.95,1]

# ╔═╡ 6a1623cc-1531-46bd-91fc-7077a73e4682
md"
### B. Binning by Land-Sea Mask Values
"

# ╔═╡ f43548df-2e32-456a-93d1-767eedcfb27e
begin
	binmat = zeros(nr,np,4,4)
	binstd = zeros(nr,np,4,4)
	binnum = zeros(nr,np,4,4)

	for igeo = 1 : 4
		for ilsm = 1 : (length(lsmvec)-1)
			rain = []
			rHDO = []
			pwgt = []
			if isone(igeo)
				ind = (lsm[1].>=lsmvec[ilsm]) .& (lsm[1].<=lsmvec[ilsm+1])
			else
				ind = isone.(ginfo[igeo-1].mask) .& (lsm[1].>=lsmvec[ilsm]) .& (lsm[1].<=lsmvec[ilsm+1])
			end
			for imo = 8 : 12
				mostr = uppercase(monthabbr(imo))
				ids  = NCDataset(datadir("wrf","2D","$(mostr)-RAINNC"))
				lon  = ids["longitude"][:]
				lat  = ids["latitude"][:]
				irin = ids["RAINNC"][ind]
				close(ids)
				ids  = NCDataset(datadir("wrf","2D","$(mostr)-HDO_RAINNC"))
				iHDO = ids["HDO_RAINNC"][ind]
				iHDO = (iHDO./irin .-1) * 1000
				close(ids)
				ids  = NCDataset(datadir("wrf","2D","$(mostr)-p_wwgt"))
				ipwg = ids["p_wwgt"][ind]
				close(ids)
				rain = vcat(rain,irin[:])
				rHDO = vcat(rHDO,iHDO[:])
				pwgt = vcat(pwgt,ipwg[:])
			end
			rainl = rain[:]
			pwgtl = pwgt[:] ./ 100
			rHDOl = rHDO[:]
	
			for ip = 1 : np, ir = 1 : nr
				try
					binmat[ir,ip,ilsm,igeo] = mean(rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir])])
					binnum[ir,ip,ilsm,igeo] = sum((pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir]))
				catch
					binmat[ir,ip,ilsm,igeo] = NaN
					binnum[ir,ip,ilsm,igeo] = NaN
				end
			end
			for ip = 1 : np, ir = 1 : nr
				try
					binstd[ir,ip,ilsm,igeo] = sum((rHDOl[(pwgtl.>pwgtb[ip]) .& (pwgtl.<pwgtu[ip]) .& (rainl.>rainb[ir]) .& (rainl.<rainu[ir])] .- binmat[ir,ip,ilsm,igeo]).^2)
				catch
					binstd[ir,ip,ilsm,igeo] = NaN
				end
			end
		end
	end
end

# ╔═╡ 4b428535-33c6-4b1d-919a-3a56de5b96b3
begin
	pplt.close()
	f4,a4 = pplt.subplots(
		ncols=8,nrows=4,axwidth=0.75,aspect=2/3,
		wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c4 = a4[1].pcolormesh(
		rainm,pwgtm,binmat[:,:,1,1]',
		cmap="viridis",levels=-100:5:-40,extend="both"
	)
	c4_2 = a4[1].pcolormesh(
		rainm,pwgtm,sqrt.(binstd[:,:,1,1])'.*NaN,
		cmap="fire",levels=0:2:20,extend="both"
	)

	for igeo = 1 : 4, ilsd = 1 : 4
		binni = binnum[:,:,ilsd,igeo]
		binii = binmat[:,:,ilsd,igeo]
		binii[binni.<=5] .= NaN
		a4[(ilsd+(igeo-1)*4)*2-1].pcolormesh(
			rainm,pwgtm,binii',
			cmap="viridis",levels=-100:5:-40,extend="both"
		)
		# a4[ii].contour(
		# 	rainm,pwgtm,binni'./sum(binni[.!isnan.(binni)])*np*nr,linestyle="--",
		# 	cmap="fire",levels=[1,2,5,10,20],extend="both",
		# )

		binii = sqrt.(binstd[:,:,ilsd,igeo]./binnum[:,:,ilsd,igeo])
		binni = binnum[:,:,ilsd,igeo]
		binii[binni.<=5] .= NaN
		a4[(ilsd+(igeo-1)*4)*2].pcolormesh(
			rainm,pwgtm,binni'./sum(binni[.!isnan.(binni)])*(nr-1)*(np-1),
			levels=(0:2:20),extend="both"
		)
		if (ilsd+(igeo-1)*4)*2 <= 8
			a4[(ilsd+(igeo-1)*4)*2].format(rtitle="$(@sprintf("%.1e",lsmvec[ilsd])) < lsm < $(@sprintf("%.1e",lsmvec[ilsd+1]))")
		end
	end

	for ax = a4
		ax.format(
			ylim=(maximum(pwgtm),minimum(pwgtm)),ylabel=L"p$_w$ / hPa",
			xlabel=L"Rain Rate / mm day$^{-1}$",
			xlocator=0:20:100,xlim=(5,75)
		)
	end

	a4[8].format(urtitle="(a) All")
	a4[16].format(urtitle="(b) PCS")
	a4[24].format(urtitle="(c) PAC")
	a4[32].format(urtitle="(d) ATL")

	lbl = L"$\delta^{2}$H / $\perthousand$"
	f4.colorbar(c4,label= L"$\mu\delta^{2}$H / $\perthousand$ " * lbl,row=2)
	f4.colorbar(c4_2,label="Probability Density",row=3)
	f4.savefig(plotsdir("04g-wprofilelandsea.png"),transparent=false,dpi=400)
	load(plotsdir("04g-wprofilelandsea.png"))
end

# ╔═╡ 825965a9-4bf6-4269-ad65-859ca220475c
(nr-1)*(np-1)

# ╔═╡ f359ef94-6519-4b56-a072-102541fc086c
md"
### C. Binning Isotopic Profiles by Land-Sea Mask Values
"

# ╔═╡ cea5964b-2065-41c2-9905-2ac7c9de8486
lsmvec2 = vcat(0,collect(10. .^(-4:0.125:-1)),0.15,0.2,0.25,0.5,0.9,0.95,0.99,1)

# ╔═╡ 5a9309ff-663a-4120-a480-508615317309
begin
	qvap_profile = zeros(50,length(lsmvec2)-1)
	δ18O_profile = zeros(50,length(lsmvec2)-1)
	δHDO_profile = zeros(50,length(lsmvec2)-1)
	pres_profile = zeros(50,length(lsmvec2)-1)
	qvap_mean    = zeros(50)
	δ18O_mean    = zeros(50)
	δHDO_mean    = zeros(50)
	pres_mean    = zeros(50)
	for imo = 8 : 12
	
	    f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-QVAPOR_ISO")
	    d3D = NCDataset(f3D)
	
		qvap = d3D["QVAPOR"][:]
		δHDO = d3D["HDO_QVAPOR"][:]
		δ18O = d3D["O18_QVAPOR"][:]
	
	    close(d3D)
	
	    f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-W")
	    d3D = NCDataset(f3D)
	
		pres = d3D["P"][:]
	
	    close(d3D)
	
		for ilsm = 1 : (length(lsmvec2)-1)
	
			ii = isone.(ginfo.mask) .& (lsm .< lsmvec2[ilsm+1]) .& (lsm .> lsmvec2[ilsm])
	
			for ilvl = 1 : 50
	
				pres_profile[ilvl,ilsm] += mean((pres[:,:,ilvl])[ii])
				qvap_profile[ilvl,ilsm] += mean((qvap[:,:,ilvl])[ii])
				δHDO_profile[ilvl,ilsm] += mean((δHDO[:,:,ilvl])[ii])
				δ18O_profile[ilvl,ilsm] += mean((δ18O[:,:,ilvl])[ii])
	
			end
	
		end

		jj = isone.(ginfo.mask) .& (lsm .< 0.5)
		for ilvl = 1 : 50
	
			pres_mean[ilvl] += mean((pres[:,:,ilvl])[jj])
			qvap_mean[ilvl] += mean((qvap[:,:,ilvl])[jj])
			δHDO_mean[ilvl] += mean((δHDO[:,:,ilvl])[jj])
			δ18O_mean[ilvl] += mean((δ18O[:,:,ilvl])[jj])

		end
	
	end
	
	pres_profile /= 500
	qvap_profile /= 5
	δ18O_profile /= 5
	δHDO_profile /= 5
	pres_mean /= 500
	qvap_mean /= 5
	δ18O_mean /= 5
	δHDO_mean /= 5
	md"Binning the vertical profiles for isotopes ..."
end

# ╔═╡ 2a1c1e7b-003f-4315-a4a7-ff50731ba17d
begin
	pres_tmp = zeros(size(pres_profile,1)+1,size(pres_profile,2))
	pres_tmp[2:(end-1),:] .= (pres_profile[1:(end-1),:] .+ pres_profile[2:(end),:])/2
	pres_tmp[1,:] .= pres_profile[1,:] .+ (pres_profile[1,:].-pres_tmp[2,:])
	pres_tmp[end,:] .= pres_profile[end,:] .+ (pres_profile[end,:].-pres_tmp[end-1,:])
	pres_plot = zeros(size(pres_profile,1)+1,size(pres_profile,2)+1)
	pres_plot[:,2:(end-1)] .= (pres_tmp[:,1:(end-1)] .+ pres_tmp[:,2:(end)])/2
	pres_plot[:,1] .= pres_tmp[:,1] .+ (pres_tmp[:,1].-pres_plot[:,2])
	pres_plot[:,end] .= pres_tmp[:,end] .+ (pres_tmp[:,end].-pres_plot[:,end-1])
	md"preparing pressure axis for plotting"
end

# ╔═╡ 5438cba9-af9c-4e04-bcb7-e83a1d1cd69d
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=2,aspect=3,axwidth=4)

	lsmplot = ones(51,length(lsmvec2)) .* vcat(0,0.5:0.25:6.5,6.5+1/3,6.5+2/3,7.5,8,8.5,9,9.5,10)'
	
	c1 = a2[1].contourf(
		lsmplot,pres_plot,(δ18O_profile ./ qvap_profile .-1) * 1000,
		levels=-100:5:-20,
		cmap="viridis",extend="both"
	)
	
	c2 = a2[2].contourf(
		lsmplot,pres_plot,(δHDO_profile ./ qvap_profile .-1) * 1000,
		levels=-625:25:-125,
		cmap="viridis",extend="both"
	)

	for ax in a2
		ax.plot([8,8],[1000,50],c="w",linestyle="--")
		ax.text(8.2,80,"Coastline",c="w")
		ax.format(
			ylim=(1000,60),yscale="log",ylabel="Pressure / hPa",
			xlim=(0,10),xlocator=xlocator=vcat(0,0.5:9.5,10),
			xlabel="Smoothed Land-Sea Mask",
			xticklabels=["0","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","0.1","0.25","0.9","0.99","1"]
		)
	end
	
	a2[1].colorbar(c1,label=L"$\delta^{18}$O / $\perthousand$")
	a2[2].colorbar(c2,label=L"$\delta^{2}$H / $\perthousand$")
	
	f2.savefig(plotsdir("04g-depletion.png"),transparent=false,dpi=400)
	load(plotsdir("04g-depletion.png"))
end

# ╔═╡ 6f975f7c-5aff-4497-b154-97b1330f912c
δ18O_spl = Spline1D(reverse(pres_mean), reverse(δ18O_mean ./ qvap_mean))

# ╔═╡ c66d845f-0900-4454-a958-f5410221fb92
δHDO_spl = Spline1D(reverse(pres_mean), reverse(δHDO_mean ./ qvap_mean))

# ╔═╡ 0b54387f-74f8-43ee-9d94-ece0baf6290c
begin
	δ18O_diff = zeros(size(δ18O_profile))
	δHDO_diff = zeros(size(δHDO_profile))
	for ii in eachindex(δ18O_diff)
		δ18O_diff[ii] = δ18O_profile[ii] / qvap_profile[ii] - δ18O_spl(pres_profile[ii])
		δHDO_diff[ii] = δHDO_profile[ii] / qvap_profile[ii] - δHDO_spl(pres_profile[ii])
	end
end

# ╔═╡ fc4bb9df-16b3-4406-873b-442bd9421b0d
begin
	pplt.close(); f3,a3 = pplt.subplots(nrows=2,aspect=4,axwidth=5)
	
	c3_1 = a3[1].pcolormesh(
		lsmplot,pres_plot,δ18O_diff * 1000,cmap="RdBu",extend="both",
		levels=vcat(-10,-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5,10),
	)
	
	c3_2 = a3[2].pcolormesh(
		lsmplot,pres_plot,δHDO_diff * 1000,cmap="RdBu",extend="both",
		levels=vcat(-10,-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5,10)*10,
	)

	for ax in a3
		ax.plot([8,8],[1000,50],c="k")
		ax.text(8.15,75,"Coastline",c="k")
		ax.format(
			ylim=(1000,50),yscale="log",ylabel="Pressure / hPa",
			xlim=(0,10),xlocator=vcat(0,0.5:9.5,10),xlabel="Smoothed Land-Sea Mask",
			xticklabels=["0","1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","0.1","0.25","0.9","0.99","1"]
		)
	end

	a3[1].format(ultitle=L"(a) $\delta^{18}$O / $\perthousand$")
	a3[2].format(ultitle=L"(b) $\delta^{2}$H / %")
	
	f3.colorbar(c3_1)
	f3.savefig(plotsdir("04g-depletion_againstmean.png"),transparent=false,dpi=400)
	load(plotsdir("04g-depletion_againstmean.png"))
end

# ╔═╡ Cell order:
# ╟─c06f21d0-e50a-11ec-2cf3-8d8ec5d3bcbb
# ╟─3c4b037f-c854-4148-ae6e-fff50bbb51cf
# ╟─4a2f8ac0-9ca0-4d33-b5f8-983c7efb0f3d
# ╟─1ae21e79-5b6d-4f31-8c62-9ed7d3ac6e79
# ╟─22dfb33e-276a-4285-81c8-5f3300c50c77
# ╟─c853d30e-5a51-47f4-a6be-88ab857253c1
# ╠═18e518f8-df88-4e57-9123-7852546af717
# ╟─2489ca4e-0765-40b0-b131-49d8a1563727
# ╟─6b55be2c-2f4e-4c41-8eb6-24e6d769df00
# ╠═c26cdc06-70f1-444a-b190-e2beb4b00ef7
# ╟─6a1623cc-1531-46bd-91fc-7077a73e4682
# ╟─f43548df-2e32-456a-93d1-767eedcfb27e
# ╟─4b428535-33c6-4b1d-919a-3a56de5b96b3
# ╠═825965a9-4bf6-4269-ad65-859ca220475c
# ╟─f359ef94-6519-4b56-a072-102541fc086c
# ╟─cea5964b-2065-41c2-9905-2ac7c9de8486
# ╟─5a9309ff-663a-4120-a480-508615317309
# ╟─2a1c1e7b-003f-4315-a4a7-ff50731ba17d
# ╟─5438cba9-af9c-4e04-bcb7-e83a1d1cd69d
# ╠═6f975f7c-5aff-4497-b154-97b1330f912c
# ╠═c66d845f-0900-4454-a958-f5410221fb92
# ╟─0b54387f-74f8-43ee-9d94-ece0baf6290c
# ╟─fc4bb9df-16b3-4406-873b-442bd9421b0d
