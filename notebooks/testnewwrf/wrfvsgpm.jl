### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ eb55ab76-4877-11ef-2539-2549bd9d61aa
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ 77f31220-ef16-4b4c-aa82-a89d603d23da
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	using GeoRegions
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow

	pplt = pyimport("proplot")
end

# ╔═╡ 255b4912-b4f3-420c-a9b0-5854da29038a
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ ddf80924-937f-4f62-afbe-b33764a36db4
begin
	ds  = NCDataset(datadir("wrf","grid.nc"))
	wln = ds["longitude"][:,:]; nx,ny = size(wln)
	wlt = ds["latitude"][:,:]
	close(ds)
end

# ╔═╡ 3942c663-ed45-4ec8-92e2-568ca89ea013
begin
	fnc1 = "wrfout_d02_2019-08-01_00:00:00"
	fnc2 = "wrfout_d02_2019-08-11_00:00:00"
	pds1  = NCDataset(datadir("wrf","wrfout","AUGp01",fnc1))
	pds2  = NCDataset(datadir("wrf","wrfout","AUGp01",fnc2))
	rain  = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O18   = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
	close(pds1)
	close(pds2)
	pds1  = NCDataset(datadir("wrftest3","zTest04",fnc1))
	pds2  = NCDataset(datadir("wrftest3","zTest04",fnc2))
	rain2 = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O182  = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
	close(pds1)
	close(pds2)
	pds1  = NCDataset(datadir("wrftest3","zTest05",fnc1))
	pds2  = NCDataset(datadir("wrftest3","zTest05",fnc2))
	rain3 = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O183  = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
	close(pds1)
	close(pds2)
	pds1  = NCDataset(datadir("wrftest3","zTest06",fnc1))
	pds2  = NCDataset(datadir("wrftest3","zTest06",fnc2))
	rain4 = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O184  = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
	close(pds1)
	close(pds2)
	pds1  = NCDataset(datadir("wrftest3","zTest07",fnc1))
	pds2  = NCDataset(datadir("wrftest3","zTest07",fnc2))
	rain5 = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O185  = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
	close(pds1)
	close(pds2)
end

# ╔═╡ 70258162-c565-4fd6-b9fc-80140ded9336
geo = GeoRegion("OTREC",path=srcdir())

# ╔═╡ a8393b07-cef0-482d-8f5a-6db938e3c87b
npd = IMERGFinalHH(start=Date(2019,8),stop=Date(2019,8,10),path=datadir(),v6=true)

# ╔═╡ a968a93b-323b-4f82-b530-019ad35e8b30
lsd = getLandSea(npd,geo); nlon = length(lsd.lon); nlat = length(lsd.lat)

# ╔═╡ 7ee9555c-7bc3-483c-8f46-f56dbda31f91
begin
	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wln[ilon,ilat].-lsd.lon.+360))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlt[ilon,ilat].-lsd.lat))
	end
	md"Finding closest IMERG points to each of the WRF points ..."
end

# ╔═╡ 5fe1f788-f439-413d-9342-a63cfe1bbd40
begin
	rain_grd1 = zeros(nlon,nlat)
	rain_grd2 = zeros(nlon,nlat)
	rain_grd3 = zeros(nlon,nlat)
	rain_grd4 = zeros(nlon,nlat)
	rain_grd5 = zeros(nlon,nlat)
	for ilat = 1 : nlat, ilon = 1 : nlon
		ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
		iprcp = @view rain[ind];  rain_grd1[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
		iprcp = @view rain2[ind]; rain_grd2[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
		iprcp = @view rain3[ind]; rain_grd3[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
		iprcp = @view rain4[ind]; rain_grd4[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
		iprcp = @view rain5[ind]; rain_grd5[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
	end
	md"Regridding/binning precipitation into IMERG Grid"
end

# ╔═╡ fa732405-8811-4ce1-92c4-f089fb9c7ae3
begin
	prcp_imerg = zeros(nlon,nlat)
	dtvec = Date(2019,8) : Day(1) : Date(2019,8,10)
	for idt in dtvec
		ids = read(npd,geo,idt)
		prcp_imerg[:,:] += sum(ids["precipitation"][:,:,:],dims=3) * 1800
		close(ids)
	end
end

# ╔═╡ d76c7c13-b3e6-49bf-88fa-fa48989689ed
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=2,ncols=3,axwidth=1.5)
	
	c = 
	axs[1].pcolormesh(wln,wlt,rain./10, levels=0:5:50,cmap="Blues",extend="both")
	axs[2].pcolormesh(wln,wlt,rain2./10,levels=0:5:50,cmap="Blues",extend="both")
	axs[3].pcolormesh(wln,wlt,rain3./10,levels=0:5:50,cmap="Blues",extend="both")
	axs[4].pcolormesh(wln,wlt,rain4./10,levels=0:5:50,cmap="Blues",extend="both")
	axs[5].pcolormesh(wln,wlt,rain5./10,levels=0:5:50,cmap="Blues",extend="both")
	axs[6].pcolormesh(lsd.lon.-360,lsd.lat,prcp_imerg'./10,levels=0:5:50,cmap="Blues",extend="both")
	
	
	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-91,-74),ylim=(-1,16))
	end

	axs[1].format(title="Original")
	axs[2].format(title="Test4")
	axs[3].format(title="Test5")
	axs[4].format(title="Test6")
	axs[5].format(title="Test7")
	axs[6].format(title="GPM IMERG")

	fig.colorbar(c,label=L"Rain Rate / mm day$^{-1}$")
	fig.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 311585c5-0763-43b1-a78c-8c7b15b2939d
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=3,ncols=2,axwidth=1.5)

	lvls = 5:5:50
	
	a2[1].pcolormesh(lsd.lon,lsd.lat,rain_grd1' ./10,levels=lvls,cmap="Blues",extend="both")
	a2[3].pcolormesh(lsd.lon,lsd.lat,rain_grd4' ./10,levels=lvls,cmap="Blues",extend="both")
	a2[5].pcolormesh(lsd.lon,lsd.lat,rain_grd5' ./10,levels=lvls,cmap="Blues",extend="both")
	
	c2 =
	a2[2].pcolormesh(lsd.lon,lsd.lat,(rain_grd1.-prcp_imerg)'./10,levels=-50:5:50,cmap="drywet",extend="both")
	a2[4].pcolormesh(lsd.lon,lsd.lat,(rain_grd4.-prcp_imerg)'./10,levels=-50:5:50,cmap="drywet",extend="both")
	a2[6].pcolormesh(lsd.lon,lsd.lat,(rain_grd5.-prcp_imerg)'./10,levels=-50:5:50,cmap="drywet",extend="both")
	
	
	for ax in a2
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-91,-74).+360,ylim=(-1,16))
	end

	sig = [
		mean(abs.(rain_grd1.-prcp_imerg)[.!isnan.(rain_grd1)]) ./ 10
		mean(abs.(rain_grd2.-prcp_imerg)[.!isnan.(rain_grd2)]) ./ 10
		mean(abs.(rain_grd3.-prcp_imerg)[.!isnan.(rain_grd3)]) ./ 10
		mean(abs.(rain_grd4.-prcp_imerg)[.!isnan.(rain_grd4)]) ./ 10
		mean(abs.(rain_grd5.-prcp_imerg)[.!isnan.(rain_grd5)]) ./ 10
	]

	a2[1].format(ultitle="Original")
	a2[2].format(urtitle=L"$\sigma = $" * "$(@sprintf("%.2f",sig[1]))")
	a2[3].format(ultitle="Test6")
	a2[4].format(urtitle=L"$\sigma = $" * "$(@sprintf("%.2f",sig[4]))")
	a2[5].format(ultitle="Test7")
	a2[6].format(urtitle=L"$\sigma = $" * "$(@sprintf("%.2f",sig[5]))")

	f2.colorbar(c, label=L"WRF Rainfall / mm day$^{-1}$",rows=[1])
	f2.colorbar(c2,label=L"WRF - IMERG / mm day$^{-1}$",rows=[2,3])
	f2.savefig("test2.png",transparent=false,dpi=200)
	load("test2.png")
end

# ╔═╡ df1c949c-d8ad-4087-aadd-e0fe0be1b098
begin
	ds2 = NCDataset(datadir("wrf2","grid.nc"))
	wln2 = ds2["longitude"][:,:]; nx2,ny2 = size(wln2)
	wlt2 = ds2["latitude"][:,:]
	close(ds2)
end

# ╔═╡ 06fbba85-c01a-4a90-bb21-22ce7c25139c
geo2 = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 3f502514-d299-4cc2-ba46-3902bfbea811
lsd2 = getLandSea(npd,geo2); nlon2 = length(lsd2.lon); nlat2 = length(lsd2.lat)

# ╔═╡ fe4aee33-f52c-47e6-8020-dd6575ff0bae
begin
	ipnt2_lon = zeros(Int,nx2,ny2)
	ipnt2_lat = zeros(Int,nx2,ny2)
	for ilat = 1 : ny2, ilon = 1 : nx2
		ipnt2_lon[ilon,ilat] = argmin(abs.(wln2[ilon,ilat].-lsd2.lon))
		ipnt2_lat[ilon,ilat] = argmin(abs.(wlt2[ilon,ilat].-lsd2.lat))
	end
	md"Finding closest IMERG points to each of the WRF points ..."
end

# ╔═╡ 75951346-2b62-4515-9a9a-dc1e86531e38
begin
	fnc3 = "wrfout_d02_2020-01-02_00:00:00"
	fnc4 = "wrfout_d02_2020-01-11_00:00:00"
	pds3  = NCDataset(datadir("wrftest3","zTest08",fnc3))
	pds4  = NCDataset(datadir("wrftest3","zTest08",fnc4))
	rain6 = pds4["RAINNC"][:,:,end] .- pds3["RAINNC"][:,:,1]
	O186  = pds4["O18_RAINNC"][:,:,end] .- pds3["O18_RAINNC"][:,:,1]
	close(pds3)
	close(pds4)
end

# ╔═╡ 3814bdb3-281f-48f2-919a-d2da96e7f8ce
begin
	rain_grd6 = zeros(nlon2,nlat2)
	for ilat = 1 : nlat2, ilon = 1 : nlon2
		ind = (ipnt2_lon.==ilon).&(ipnt2_lat.==ilat)
		iprcp = @view rain6[ind]; rain_grd6[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
	end
	md"Regridding/binning precipitation into IMERG Grid"
end

# ╔═╡ aac270c7-8641-4c2e-bf7d-c264a446ecc0
begin
	prcp2_imerg = zeros(nlon2,nlat2)
	dtvec2 = Date(2020,1,2) : Day(1) : Date(2020,1,11)
	for idt in dtvec2
		ids = read(npd,geo2,idt)
		prcp2_imerg[:,:] += sum(ids["precipitation"][:,:,:],dims=3) * 1800
		close(ids)
	end
end

# ╔═╡ b59848b0-b3ef-4b88-9832-a13e50f6ef29
begin
	pplt.close(); f3,a3 = pplt.subplots(ncols=3,nrows=2,axwidth=1.5)
	
	a3[1].pcolormesh(lsd2.lon,lsd2.lat,rain_grd6'./10,levels=lvls./2,cmap="Blues",extend="both")
	a3[2].pcolormesh(lsd2.lon,lsd2.lat,prcp2_imerg'./10,levels=lvls./2,cmap="Blues",extend="both")
	c3 = a3[3].pcolormesh(lsd2.lon,lsd2.lat,(rain_grd6-prcp2_imerg)'./10,levels=-20:2:20,cmap="drywet",extend="both")

	δ1 = (O18./rain .-1).*1000; δ1[rain .< 50] .= NaN
	δ2 = (O186./rain6 .-1).*1000; δ2[rain6 .< 50] .= NaN
	
	c4 = a3[4].pcolormesh(wln,wlt,δ1,levels=-15:0,cmap="viridis",extend="both")
	a3[5].pcolormesh(wln2,wlt2,δ2,levels=-15:0,cmap="viridis",extend="both")
	
	for ax in a3
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-95,-65),ylim=(-10,20),suptitle="1-10 Jan 2020")
	end

	f3.colorbar(c3,rows=[1],label=L"WRF - IMERG / mm day$^{-1}$")
	f3.colorbar(c4,rows=[2],label=L"$\delta / \perthousand$")
	f3.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─eb55ab76-4877-11ef-2539-2549bd9d61aa
# ╠═77f31220-ef16-4b4c-aa82-a89d603d23da
# ╠═255b4912-b4f3-420c-a9b0-5854da29038a
# ╠═ddf80924-937f-4f62-afbe-b33764a36db4
# ╠═3942c663-ed45-4ec8-92e2-568ca89ea013
# ╟─70258162-c565-4fd6-b9fc-80140ded9336
# ╟─a8393b07-cef0-482d-8f5a-6db938e3c87b
# ╠═a968a93b-323b-4f82-b530-019ad35e8b30
# ╠═7ee9555c-7bc3-483c-8f46-f56dbda31f91
# ╟─5fe1f788-f439-413d-9342-a63cfe1bbd40
# ╟─d76c7c13-b3e6-49bf-88fa-fa48989689ed
# ╠═fa732405-8811-4ce1-92c4-f089fb9c7ae3
# ╠═311585c5-0763-43b1-a78c-8c7b15b2939d
# ╠═df1c949c-d8ad-4087-aadd-e0fe0be1b098
# ╠═06fbba85-c01a-4a90-bb21-22ce7c25139c
# ╠═3f502514-d299-4cc2-ba46-3902bfbea811
# ╠═fe4aee33-f52c-47e6-8020-dd6575ff0bae
# ╠═75951346-2b62-4515-9a9a-dc1e86531e38
# ╠═3814bdb3-281f-48f2-919a-d2da96e7f8ce
# ╠═aac270c7-8641-4c2e-bf7d-c264a446ecc0
# ╟─b59848b0-b3ef-4b88-9832-a13e50f6ef29
