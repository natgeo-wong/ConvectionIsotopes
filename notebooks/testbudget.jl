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
	using Dates
	using GeoRegions
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
# 04a. Moisture Budget from WRF
"

# ╔═╡ 5245d696-7392-402c-9364-8bdb45305b80
function smooth(data::AbstractArray;days=0)

	if !iszero(days)
		hrs = days * 24
		nt = length(data) - hrs
		ndata = zeros(nt)
	
		for it = 1 : nt
			for ii = 0 : (hrs-1)
				ndata[it] += data[it+ii]  / hrs
			end
		end
	else
		ndata = deepcopy(data)
	end

	return ndata

end

# ╔═╡ d447fdfd-b5c5-44f7-a2a2-161d5bb157b5
geo = GeoRegion("OTREC_wrf_stn12",path=srcdir())

# ╔═╡ 43bb0b2c-995a-4829-a30e-b953651305de
begin
	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-QBUDGET-20190801_20201231.nc"))
	t  = ds["time"][:]
	p  = ds["P"][:]; nt = length(p)
	e  = ds["E"][:]
	∇  = ds["∇"][:] * 3600
	Δ  = ds["ΔWVP"][:]
	p  = smooth(ds["P"][:],days=7); nt = length(p)
	e  = smooth(ds["E"][:],days=7)
	∇  = smooth(ds["∇"][:],days=7) * 3600
	Δ  = smooth(ds["ΔWVP"][:],days=7)
	t  = t[1:nt]
	close(ds)
	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-HDO_QBUDGET-20190801_20201231.nc"))
	HDOp = ds["HDO_P"][:]
	HDOe = ds["HDO_E"][:]
	HDOp = smooth(ds["HDO_P"][:],days=7)
	HDOe = smooth(ds["HDO_E"][:],days=7)
	close(ds)
	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-O18_QBUDGET-20190801_20201231.nc"))
	O18p = ds["O18_P"][:]
	O18e = ds["O18_E"][:]
	O18p = smooth(ds["O18_P"][:],days=7)
	O18e = smooth(ds["O18_E"][:],days=7)
	close(ds)
	# ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-∇decompose-20190801_20201231.nc"))
	# ADV  = ds["ADV"][:] * 3600
	# DIV  = ds["DIV"][:] * 3600
	# ADV  = smooth(ds["ADV"][:],days=30) * 3600
	# DIV  = smooth(ds["DIV"][:],days=30) * 3600
	close(ds)
	ds = NCDataset(datadir("wrf3","processed","$(geo.ID)-p_wwgt-daily-20190801_20201231-smooth_30days.nc"))
	dt = ds["time"][:]
	σ  = ds["σ_wwgt"][:]
	close(ds)
end

# ╔═╡ 509f6ffc-2d24-4a85-9e64-d0dd49108e9e
HDOe

# ╔═╡ 839bf2a5-0f04-45d1-b683-99f275619ff5
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=3,axwidth=6)
	
	axs[1].plot(HDOe./e)
	# axs[1].format(xlim=(Date(2019,8,1),Date(2020,12,31)))/
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 5a230aea-8970-4d96-8aaa-85806f9a1ab8
begin
	maxx = zeros(length(∇))
	for it = 1 : length(∇)
		maxx[it] = max(abs(p[it]),abs(e[it]),abs(∇[it]),abs(Δ[it]))
	end
	maxx = abs.(p .- e .+ ∇ .+ Δ) ./ maxx
end

# ╔═╡ a38d959a-38ae-437f-9af0-86b805449174
begin
	pplt.close(); f3,a3 = pplt.subplots(axwidth=1.5)
	
	# a3[1].scatter(∇[1:length(∇calc)],∇calc,s=5)
	# a3[1].scatter(∇,ADV.+DIV,s=5)
	# a3[1].scatter(p,∇,s=5)
	# a3[1].scatter(p,DIV,s=5)
	# a3[1].scatter(DIV,ADV,s=1)
	# a3[1].scatter(∇[1:744],∇calc,s=1)
	a3[1].scatter(∇,(ADV.+DIV),s=1)
	a3[1].plot([-1000,1000],[-1000,1000],lw=0.5,c="k")
	a3[1].plot([-1000,1000],[1000,-1000],lw=0.5,c="k")
	a3[1].format(xlim=(-5,5),ylim=(-5,5),xlocator=(-2:2),ylocator=(-2:2))
	
	f3.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 2c7dc5d6-00b5-4136-ba46-ce3b524e2df3
DIV

# ╔═╡ ffb80048-2557-405a-aab8-5cc2fef9bfad
length(ADV.+DIV)

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╠═bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─5245d696-7392-402c-9364-8bdb45305b80
# ╠═d447fdfd-b5c5-44f7-a2a2-161d5bb157b5
# ╠═509f6ffc-2d24-4a85-9e64-d0dd49108e9e
# ╠═43bb0b2c-995a-4829-a30e-b953651305de
# ╠═839bf2a5-0f04-45d1-b683-99f275619ff5
# ╠═5a230aea-8970-4d96-8aaa-85806f9a1ab8
# ╠═a38d959a-38ae-437f-9af0-86b805449174
# ╠═2c7dc5d6-00b5-4136-ba46-ce3b524e2df3
# ╠═ffb80048-2557-405a-aab8-5cc2fef9bfad
