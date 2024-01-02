### A Pluto.jl notebook ###
# v0.19.36

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
	using ERA5Reanalysis
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# Figure 1. Theoretical $p_\omega$ Schematic
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 75925166-d1e8-450e-bae6-29affd25d635
md"
### A. Defining Datasets, Variables and Regions
"

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Monthly(start=Date(1979,1,1),stop=Date(2021,12,31),path=datadir())

# ╔═╡ 9473b30b-787a-495b-bdfe-b386c9995761
evar = SingleVariable("p_wwgt")

# ╔═╡ c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
evar_pre = SingleVariable("sp")

# ╔═╡ c7d0ae59-8d78-4946-8348-9fa570197b0b
egeo = ERA5Region("TRP")

# ╔═╡ cb7c6118-e25b-462a-84a5-751ea0682b52
elsd = getLandSea(e5ds,egeo,smooth=true,σlon=2,σlat=2)

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
begin
	tmp = zeros(length(elsd.lon),length(elsd.lat),12)
	wp  = zeros(length(elsd.lon),length(elsd.lat),12)
	sp  = zeros(length(elsd.lon),length(elsd.lat),12)
	cnt = zeros(length(elsd.lon),length(elsd.lat),12)
	for dt in e5ds.start : Year(1) : e5ds.stop
		wpds = read(e5ds,evar,egeo,dt,quiet=true)
		spds = read(e5ds,evar_pre,egeo,dt,quiet=true)
		tmp[:,:,:] .= nomissing(wpds[evar.ID][:,:,:],0) / 100
		wp[:,:,:]  += tmp
		cnt[:,:,:] += Int.(.!iszero.(tmp))
		sp[:,:,:]  += nomissing(spds[evar_pre.ID][:,:,:]) / 100
		close(wpds)
		close(spds)
	end
	wp = wp ./ cnt; wp_yr = dropdims(mean(wp,dims=3),dims=3)
	sp = sp ./ 43;  sp_yr = dropdims(mean(sp,dims=3),dims=3)
	wp_yr[elsd.z.>500] .= NaN
	md"Loading and filtering w-weighted pressure and sigma data and counts ..."
end

# ╔═╡ 459cca84-c2bb-4ebe-8f2d-52e576d62906
md"
### C. Calculating Sample Depletion Profiles from WRF
"

# ╔═╡ 7c954c2b-bfdc-4efd-a3b4-8038e66e1ee5
begin
	wds = NCDataset(datadir("fig2wrfdata.nc"))
	q   = wds["QVAPOR"][:,:]
	qh  = wds["HDO_QVAPOR"][:,:]
	p   = wds["P"][:,:]
	ps  = wds["PSFC"][:]
	close(wds)
	md"Loading the WRF 2D and 3D datasets"
end

# ╔═╡ 6231902f-897b-4fd1-b6f6-856e11f83168
begin
	μq  = dropdims(mean( q[ps.>100000,:],dims=1),dims=1)
	μqh = dropdims(mean(qh[ps.>100000,:],dims=1),dims=1)
	μp  = dropdims(mean( p[ps.>100000,:],dims=1),dims=1)
	μδ  = (μqh ./ μq .- 1) * 1000
	md"Calculating the average depletion over ocean/lowland regions ..."
end

# ╔═╡ 00f069bb-0834-4859-94fa-7b0c146e3a11
md"
### D. Plotting the climatology
"

# ╔═╡ 4342521a-6a0a-4d30-a5f3-975514d9006f
begin
	pplt.close(); fig,axs = pplt.subplots([
		[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
		[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
		[0,2,2,2,2,2,2,0,0,3,3,3,3,3,3,0],
		[0,2,2,2,2,2,2,0,0,3,3,3,3,3,3,0],
		[0,2,2,2,2,2,2,0,0,3,3,3,3,3,3,0],
		[0,2,2,2,2,2,2,0,0,3,3,3,3,3,3,0],
		[6,4,4,4,4,4,4,0,0,5,5,5,5,5,5,7],
		[6,4,4,4,4,4,4,0,0,5,5,5,5,5,5,7],
		[6,4,4,4,4,4,4,0,0,5,5,5,5,5,5,7],
	],aspect=7,axwidth=5,wspace=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,1],sharex=0)
	
	c1 = axs[1].pcolormesh(elsd.lon,elsd.lat,wp_yr',levels=(9:0.5:15).*50,extend="both",cmap="drywet_r")

	for ii in 1 : 3
		axs[ii].pcolormesh(elsd.lon,elsd.lat,wp_yr',levels=(9:0.5:15).*50,extend="both",cmap="drywet_r")
		axs[ii].plot(x,y,c="k",lw=0.5)
		axs[ii].format(ylabel=L"Latitude / $\degree$",ylocator=-10:10:20,yminorlocator=-10:5:20,xminorlocator=90:5:300)
	end

	axs[1].plot([100,130,130,100,100],[-5,-5,15,15,-5],c="k",linestyle="--")
	axs[1].plot([260,290,290,260,260],[-5,-5,15,15,-5],c="k",linestyle="--")
	axs[1].text(133,10,"(b)")
	axs[1].text(251,-3,"(c)")

	axs[2].plot([0,360],[0,0],c="r",linestyle="--")
	axs[3].plot([0,360],[5,5],c="r",linestyle="--")
	
	axs[1].format(xlim=(90,300),xlocator=60:30:330,ylim=(-10,20),ltitle="(a) Tropical Pacific Basin")
	axs[2].format(xlim=(100,130),xlocator=100:5:130,ylim=(-5,15),ylocator=-15:5:30,yminorlocator=-15:20,xminorlocator=100:130,ltitle="(b) West Pacific",ultitle="(i)")
	axs[3].format(xlim=(260,290),xlocator=260:5:290,ylim=(-5,15),ylocator=-15:5:20,yminorlocator=-15:20,xminorlocator=260:290,ltitle="(c) East Pacific",ultitle="(i)")

	c2 = axs[4].contourf(elsd.lon,μp/100,repeat(μδ,outer=(1,1440)),levels=-600:50:-150,cmap="viridis",extend="both",cmap_kw=Dict("left"=>0.2))
	axs[4].fill_between(elsd.lon,1000,sp_yr[:,121],c="brown")
	axs[4].plot([0,360],ones(2)*10. ^2.65,c="r",linestyle=":")
	axs[4].format(xlim=(100,130),xlabel=L"Longitude / $\degree$",xlocator=100:5:130,xminorlocator=100:130,ultitle="(ii)")

	axs[5].contourf(elsd.lon,μp/100,repeat(μδ,outer=(1,1440)),levels=-600:50:-150,cmap="viridis",extend="both",cmap_kw=Dict("left"=>0.2))
	axs[5].fill_between(elsd.lon,1000,sp_yr[:,101],c="brown")
	axs[5].plot([0,360],ones(2)*10. ^2.8,c="r",linestyle=":")
	axs[5].format(xlim=(260,290),xlabel=L"Longitude / $\degree$",xlocator=260:5:290,xminorlocator=260:290,ultitle="(ii)")

	axs[6].plot(0.8sin.((0:0.01:1)*pi),10. .^(LinRange(2.3,3,101)),c="k")
	axs[6].plot([0,0],10. .^[2,2.3],c="k")
	axs[6].plot([-5,5],ones(2)*10. ^2.65,c="r",linestyle=":")
	axs[6].format(ultitle="(iii)")
	
	axs[7].plot(0.2sin.((0:0.01:1)*pi).-0.5sin.((0:0.01:1)*2pi),10. .^(LinRange(2.3,3,101)),c="k")
	axs[7].plot([0,0],10. .^[2,2.3],c="k")
	axs[7].plot([-5,5],ones(2)*10. ^2.8,c="r",linestyle=":")
	axs[7].format(ultitle="(iii)")

	for ii in 4 : 5
		axs[ii].format(ylim=(1000,100),yscale="log",xlocator=0:5:360,ylabel="Pressure / hPa",xlabel=L"Longitude / $\degree$")
	end

	for ii in 6 : 7
		axs[ii].format(xlim=(-2,2),xticks=[0],yloc="none",xloc="b",xlabel=L"w")
	end

	fig.colorbar(c1,locator=400:50:800,minorlocator=(5:0.5:15).*50,label=L"$p_\omega$",length=0.8,rows=(1,6))
	fig.colorbar(c2,locator=[-600,-150],label=L"$\delta_H$",ticklabels=[L"-",L"+"],rows=(7,9))
	fig.savefig(plotsdir("fig1-pomegaclimatology.png"),transparent=false,dpi=400)
	load(plotsdir("fig1-pomegaclimatology.png"))
end

# ╔═╡ cc7368e3-9591-4bf6-a907-67c585106a8c
md"
### E. Loading some 2D Cloud images
"

# ╔═╡ ee25250d-ddad-4242-9536-d3d0f72e0f85
begin
	cds = NCDataset(datadir("cloud.nc"))
	xx  = cds["x"][:]
	zz  = cds["z"][:]
	qt  = cds["QN"][:,:,end] .+ cds["QP"][:,:,end]
	qt[qt.<0.1] .= NaN
	close(cds)
end

# ╔═╡ c07e34e0-01f5-4e0e-afa1-bd3d5f157a87
begin
	pplt.close(); ftmp,atmp = pplt.subplots(ncols=2,aspect=2,axwidth=4,sharex=0)
	
	atmp[1].contourf(xx[8800:9600]/1000,zz/1000,qt[8800:9600,:]',cmap="greys")
	atmp[2].contourf(xx[15950:16250]/1000,zz/1000,qt[15950:16250,:]',cmap="greys")

	for ax in atmp
		ax.format(ylim=(0,12),grid=false,xticks="null",yticks="null")
	end

	ftmp.savefig(plotsdir("fig1-clouds.png"),transparent=true,dpi=400)
	load(plotsdir("fig1-clouds.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─9473b30b-787a-495b-bdfe-b386c9995761
# ╟─c9d1e411-a65b-41f7-949b-cf2e32b2dfb4
# ╟─c7d0ae59-8d78-4946-8348-9fa570197b0b
# ╟─cb7c6118-e25b-462a-84a5-751ea0682b52
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╟─acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
# ╟─459cca84-c2bb-4ebe-8f2d-52e576d62906
# ╟─7c954c2b-bfdc-4efd-a3b4-8038e66e1ee5
# ╟─6231902f-897b-4fd1-b6f6-856e11f83168
# ╟─00f069bb-0834-4859-94fa-7b0c146e3a11
# ╟─4342521a-6a0a-4d30-a5f3-975514d9006f
# ╟─cc7368e3-9591-4bf6-a907-67c585106a8c
# ╟─ee25250d-ddad-4242-9536-d3d0f72e0f85
# ╟─c07e34e0-01f5-4e0e-afa1-bd3d5f157a87
