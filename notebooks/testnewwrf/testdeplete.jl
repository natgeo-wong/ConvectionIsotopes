### A Pluto.jl notebook ###
# v0.20.3

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
	wln1 = ds["longitude"][:,:]
	wlt1 = ds["latitude"][:,:]
	close(ds)
	begin
	ds  = NCDataset(datadir("wrf2","grid.nc"))
	wln2 = ds["longitude"][:,:]
	wlt2 = ds["latitude"][:,:]
	close(ds)
end
end

# ╔═╡ 3942c663-ed45-4ec8-92e2-568ca89ea013
begin
	fnc1 = "wrfout_d02_2019-08-01_00:00:00"
	fnc2 = "wrfout_d02_2019-08-11_00:00:00"
	pds1  = NCDataset(datadir("wrf","wrfout","AUGp01",fnc1))
	pds2  = NCDataset(datadir("wrf","wrfout","AUGp01",fnc2))
	rain1 = pds2["RAINNC"][:,:,end] .- pds1["RAINNC"][:,:,1]
	O181  = pds2["O18_RAINNC"][:,:,end] .- pds1["O18_RAINNC"][:,:,1]
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

# ╔═╡ d76c7c13-b3e6-49bf-88fa-fa48989689ed
begin
	pplt.close(); fig,axs = pplt.subplots([[1,2,3],[0,4,5]],axwidth=1.5)

	δ1 = (O181./rain1 .-1).*1000; δ1[rain1 .< 50] .= NaN
	δ2 = (O182./rain2 .-1).*1000; δ2[rain2 .< 50] .= NaN
	δ3 = (O183./rain3 .-1).*1000; δ3[rain3 .< 50] .= NaN
	δ4 = (O184./rain4 .-1).*1000; δ4[rain4 .< 50] .= NaN
	δ5 = (O185./rain5 .-1).*1000; δ5[rain5 .< 50] .= NaN
	
	c = 
	axs[1].pcolormesh(wln1,wlt1,δ1,levels=-15:15,cmap="drywet",extend="both")
	axs[2].pcolormesh(wln1,wlt1,δ2,levels=-15:15,cmap="drywet",extend="both")
	axs[3].pcolormesh(wln1,wlt1,δ3,levels=-15:15,cmap="drywet",extend="both")
	axs[4].pcolormesh(wln1,wlt1,δ4,levels=-15:15,cmap="drywet",extend="both")
	axs[5].pcolormesh(wln1,wlt1,δ5,levels=-15:15,cmap="drywet",extend="both")
	
	
	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-91,-74),ylim=(-1,16),xlocator=-90:5:-75,suptitle="1-10 Aug 2019")
	end

	axs[1].format(title="Original")
	axs[2].format(title="Test4")
	axs[3].format(title="Test5")
	axs[4].format(title="Test6")
	axs[5].format(title="Test7")

	fig.colorbar(c,label=L"$\delta^{18}$O / $\perthousand$")
	fig.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─eb55ab76-4877-11ef-2539-2549bd9d61aa
# ╠═77f31220-ef16-4b4c-aa82-a89d603d23da
# ╠═255b4912-b4f3-420c-a9b0-5854da29038a
# ╠═ddf80924-937f-4f62-afbe-b33764a36db4
# ╠═3942c663-ed45-4ec8-92e2-568ca89ea013
# ╠═d76c7c13-b3e6-49bf-88fa-fa48989689ed
