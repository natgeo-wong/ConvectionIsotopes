### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ColombiaIsotope"
	using DelimitedFiles
	using ERA5Reanalysis
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ColumbiaIsotope project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 04a. Exploring the Weighted Column-Mean Pressure over OTREC

We have calculated gridded values for the vertical-velocity-weighted column-mean pressure ($p_{wwgt}$) over the OTREC region we defined for the years 2013 to 2021.  Now in this notebook, we do a basic exploration of $p_{wwgt}$ on a monthly basis.  We also explore $\sigma_{wwgt} = p_{wwgt} / p_s$, where $p_s$ is the surface pressure, which allows us to adjust for topography.
"

# ╔═╡ 5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 30d2be4a-9bdd-4939-8785-413cd7a59f78
md"
### A. Defining the Datasets, Variables and Regions
"

# ╔═╡ 5981e960-e7b4-40a9-9bfe-2ef90220e225
e5ds = ERA5Monthly(
	dtbeg=Date(2013,1,1),dtend=Date(2013,12,31),
	eroot="/n/kuangdss01/lab"
)

# ╔═╡ a9b15e68-ae2d-4698-a650-2f5f24c8888a
evar_wp = SingleVariable("p_wwgt")

# ╔═╡ 0d6e5046-273d-4e8b-aa4e-8ffa8a46fc9e
evar_sp = SingleVariable("sp")

# ╔═╡ 43cd84a6-53b8-4e5c-b327-3d0c30f90e5c
egeo = GeoRegion("OTREC")

# ╔═╡ 228728ed-bfab-405e-83f8-8236d6952c90
ereg = ERA5Region(egeo)

# ╔═╡ f0751d18-0d6b-4d83-a79b-21d463972edf
md"
### B. Loading a Sample Dataset
"

# ╔═╡ 41157e32-a476-44aa-92a5-b52242c7d3de
begin
	ds_wp,lon,lat = read(e5ds,evar_wp,ereg,Date(2018,1,1),lonlat=true)
	ds_sp,_,_ = read(e5ds,evar_sp,ereg,Date(2018,1,1),lonlat=true)
	mo = 4
	p_wwgt = nomissing(ds_wp[evar_wp.varID][:,:,mo],NaN)
	p_sfc  = nomissing(ds_sp[evar_sp.varID][:,:,mo],NaN)
	close(ds_wp)
	close(ds_sp)
	md"We load the surface pressure and weighted column-mean pressure"
end

# ╔═╡ 92d2d9bf-9f10-4349-9a87-2d2e072a20db
elon = 282

# ╔═╡ b401447c-acf6-4eb3-b50b-24e5afaa850d
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=2)

	c = axs[1].contourf(
		lon,lat,p_wwgt' ./ 100,
		levels=(3:0.5:9) .*100,cmap="drywet_r",extend="both"
	)
	axs[1].plot(x,y,c="k")
	axs[1].plot([elon,elon],[-10,10],c="r")
	axs[1].scatter(ones(12)*elon,vcat(-10:2:10,7),s=10,zorder=4)
	axs[1].colorbar(c,loc="r",label=L"$p_w$ / hPa")
	axs[1].format(ltitle="(a) Weighted Pressure")
	
	c = axs[2].contourf(
		lon,lat,(p_wwgt./p_sfc)',
		levels=(3:0.5:9) ./ 10,cmap="drywet_r",extend="both"
	)
	axs[2].plot(x,y,c="k")
	axs[2].plot([elon,elon],[-10,10],c="r")
	axs[2].scatter(ones(12)*elon,vcat(-10:2:10,7),s=10,zorder=4)
	axs[2].colorbar(c,loc="r",label=L"$\sigma_{p_w}$")
	axs[2].format(ltitle=L"(b) Weighted $\sigma$")

	for ax in axs
		ax.format(
			xlim=(260,305),ylim=(-15,30),
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$"
		)
	end
	
	fig.savefig(plotsdir("04a-psig_wwtg.png"),dpi=300,transparent=false)
	load(plotsdir("04a-psig_wwtg.png"))
end

# ╔═╡ e84c3b61-e52b-4c52-b7ca-92f281436237
md"Because many of our station datasets occur inland where the topography is non-negligible, as opposed to Torri et al. (2015) where they use the raw values of weighted column-mean pressure, I propose to normalize these values by the surface pressure.  This would allow us to do generalize our comparisons over land and ocean."

# ╔═╡ 4cea704d-822f-4aac-bfdf-6c63f98b2dc9
md"
### C. Sample profiles of vertical velocity
"

# ╔═╡ 1a899eaa-e28a-4260-b960-bf84e05beb43
begin
	nlon = length(lon)
	nlat = length(lat)
	wair = zeros(nlon,nlat,32)
	md"Preallocating wair arrays ..."
end

# ╔═╡ f66d4346-b937-4d2c-8132-c159c8bf7738
begin
	plvl = sort(era5Pressures())
	plvl = plvl[plvl.>=10]
	for i = 1 : 32
		evar_w = PressureVariable("w",hPa=plvl[i])
		ds_w = read(e5ds,evar_w,ereg,Date(2013,1,1))
		wair[:,:,i] = nomissing(ds_w[evar_w.varID][:,:,mo],NaN)
	end
	md"Loading vertical velocity data ..."
end

# ╔═╡ 8b7eeea3-2059-4857-87a5-29e1363d759c
begin
	wprf = zeros(34,12)
	pprf = zeros(34,12)
	pprf[2:(end-1),:] .= plvl 
	for elat = vcat(-10:2:10)
		ilon = argmin(abs.(lon.-elon))
		ilat = argmin(abs.(lat.-elat))
		ind  = Int((elat+10)/2)
		wprf[2:(end-1),ind+1] = wair[ilon,ilat,:]
		pprf[end,ind+1] = p_sfc[ilon,ilat] / 100
		pprf[:,ind+1] = pprf[:,ind+1] / pprf[end,ind+1]
	end

	for elat = 7
		ilon = argmin(abs.(lon.-elon))
		ilat = argmin(abs.(lat.-elat))
		wprf[2:(end-1),end] = wair[ilon,ilat,:]
		pprf[end,end] = p_sfc[ilon,ilat] / 100
		pprf[:,end] = pprf[:,end] / pprf[end,end]
	end
	md"Retrieving vertical velocity profile for (lon,lat) coordinates (270,-10 to 10)"
end

# ╔═╡ d58d64a6-99fb-43a6-973e-c67e8a13746a
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=2,ncols=6,aspect=1/2,axwidth=0.75)

	vlat = vcat(-10:2:10,7)
	for ii = 1 : 12
		a2[ii].plot(wprf[:,ii],pprf[:,ii])
		a2[ii].plot([0.3,-0.3],ones(2) * (p_wwgt./p_sfc)[
			argmin(abs.(lon.-elon)),
			argmin(abs.(lat.-vlat[ii]))
		])
		a2[ii].format(ultitle="lat = $(vlat[ii])" * L"$\degree$")
	end

	for ax in a2
		ax.format(
			xlim=(0.3,-0.3),xlabel=L"$\omega$ / hPa s$^{-1}$",
			ylim=(1,0),ylabel="Pressure / hPa",
			suptitle=L"Vertical Velocity Profile for lon = 270$\degree$"
		)
	end
	
	f2.savefig(plotsdir("04a-sig_wwtgprofile.png"),transparent=false,dpi=300)
	load(plotsdir("04a-sig_wwtgprofile.png"))
end

# ╔═╡ 39bff452-2d1b-44b6-b98a-632ee558cad4
md"
These profiles show that even in regions where there is heavy rainfall and strong ascent, that the atmosphere is bottom-heavy in nature.  This corroborates with results by Back and Bretherton (2006) in the eastern Pacific.
"

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─30d2be4a-9bdd-4939-8785-413cd7a59f78
# ╟─5981e960-e7b4-40a9-9bfe-2ef90220e225
# ╟─a9b15e68-ae2d-4698-a650-2f5f24c8888a
# ╟─0d6e5046-273d-4e8b-aa4e-8ffa8a46fc9e
# ╠═43cd84a6-53b8-4e5c-b327-3d0c30f90e5c
# ╟─228728ed-bfab-405e-83f8-8236d6952c90
# ╟─f0751d18-0d6b-4d83-a79b-21d463972edf
# ╠═41157e32-a476-44aa-92a5-b52242c7d3de
# ╠═92d2d9bf-9f10-4349-9a87-2d2e072a20db
# ╟─b401447c-acf6-4eb3-b50b-24e5afaa850d
# ╟─e84c3b61-e52b-4c52-b7ca-92f281436237
# ╟─4cea704d-822f-4aac-bfdf-6c63f98b2dc9
# ╟─1a899eaa-e28a-4260-b960-bf84e05beb43
# ╟─f66d4346-b937-4d2c-8132-c159c8bf7738
# ╟─8b7eeea3-2059-4857-87a5-29e1363d759c
# ╟─d58d64a6-99fb-43a6-973e-c67e8a13746a
# ╟─39bff452-2d1b-44b6-b98a-632ee558cad4
