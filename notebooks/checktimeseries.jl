### A Pluto.jl notebook ###
# v0.19.45

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
	using NCDatasets
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
	ds  = NCDataset(datadir("wrf2","grid.nc"))
	lon = ds["longitude"][:,:]
	lat = ds["latitude"][:,:]
	close(ds)
end

# ╔═╡ 1d428bbb-8c71-4c45-a2b9-2c44a6e18afb
md"
### Calculation of Vertical Modes
"

# ╔═╡ 9fe2ac23-2e92-4a3b-9436-7f2f2f6876b0
function calcvertmodes(p,w,n=2,ptop=100e2)

	nt = size(w,2)
	c  = zeros(n,nt)

	for it = 1 : nt

		ip = @views p[:,it]; ii = ip .> ptop; nz = sum(ii)
		iip = ip[ii]
		iiw = w[ii,it]

		pd = iip[end] - iip[1]

		for iz = 2 : nz, i = 1 : n

			c[i,it] +=  (iiw[iz-1] + iiw[iz]) / 2 * (iip[iz] - iip[iz-1]) * sin(pi*(iip[iz] + iip[iz-1] - iip[1])*i/2/pd)

		end

		c[:,it] ./= pd

	end

	return c

end

# ╔═╡ c2c8bc39-af9c-4e09-a1fc-ec4f274182bb
md"
### Rainfall extraction
"

# ╔═╡ a0f24335-cd66-42d1-92d4-1e63115f3af8
begin
	rds = NCDataset(datadir("wrf2","2D","RAINNC-daily-smooth_07days.nc"))
	prcp = rds["RAINNC"][:,:,4:(end-3)]
	close(rds)
end

# ╔═╡ ce3de4d4-f9b0-4e06-be86-68e1e1d1d1c4
begin
	ods = NCDataset(datadir("wrf2","2D","O18_RAINNC-daily-smooth_07days.nc"))
	O18 = ods["O18_RAINNC"][:,:,4:(end-3)]
	close(ods)
end

# ╔═╡ e240d50c-a2c9-4510-8847-83de036b7f28
geo = GeoRegion("OTREC_wrf_stn07",path=srcdir())

# ╔═╡ 3942c663-ed45-4ec8-92e2-568ca89ea013
begin
	pds   = NCDataset(datadir("wrf2","processed","$(geo.ID)-p_wwgt-daily-smooth_07days.nc"))
	w = pds["W"][:,4:(end-3)]
	p = pds["P"][:,4:(end-3)]
	pwwgt = pds["p_wwgt"][4:(end-3)]
	t = pds["time"].var[4:(end-3)]' .* ones(size(p)); 
	close(pds)
end

# ╔═╡ 745ce483-74f7-4fc8-9944-c89d3abe098a
coef = calcvertmodes(p,w)

# ╔═╡ abfe3dff-e445-4426-b963-3776ba310101
ggrd = RegionGrid(geo,lon,lat)

# ╔═╡ 2f48892c-31b2-443d-b788-1fbc6e3371ac
begin
	prcp_geo = zeros(size(prcp,3))
	O18_geo = zeros(size(prcp,3))
	mask = ggrd.mask; mask[isnan.(mask)] .= 0
	for it = 1 : size(prcp,3)
		prcp_geo[it] = sum(prcp[:,:,it].*mask) / sum(mask)
		O18_geo[it] = sum(O18[:,:,it].*mask) / sum(mask)
	end
end

# ╔═╡ d76c7c13-b3e6-49bf-88fa-fa48989689ed
begin
	pplt.close(); fig,axs = pplt.subplots([[1],[1],[1],[2]],aspect=2)

	c = axs[1].pcolormesh(t,p/100,w)
	axs[1].format(ylim=(1000,0))
	
	# axs[2].plot(t[1,:],pwwgt/100/500 .-1)

	pwwgt_filter = deepcopy(pwwgt)
	pwwgt_filter[prcp_geo.<(5)] .= NaN
	# axs[1].plot(t[1,:],pwwgt_filter/100)

	axs[2].plot(t[1,:],(O18_geo./prcp_geo .- 1)*500)
	axs[2].plot(t[1,:],prcp_geo)
	# axs[2].plot(t[1,:],O18_geo .- prcp_geo)
	# ang = atan.(coef[2,:]./coef[1,:])./pi*2
	# a = coef[2,:]
	# ang[a.<0] .= NaN
	# axs[2].plot(t[1,:],coef[2,:]./(abs.(coef[1,:]).+abs.(coef[2,:])))
	axs[1].plot(t[1,:],coef[2,:]./(abs.(coef[1,:]).+abs.(coef[2,:]))*450 .+ 550)
	# axs[2].plot(t[1,:],ang)
	# axs[2].plot(t[1,:],coef[2,:])
	# axs[2].format(ylim=(-1.5,1.5))

	fig.colorbar(c)
	fig.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 4ffbd55d-93e8-4956-9489-a34fe3249fa5


# ╔═╡ 692ac40e-bd85-4c77-858f-a2415731afac
abs.(coef[1,:].+coef[2,:])

# ╔═╡ 743ca220-78de-4878-98b9-785611bd56f7
mean(pwwgt_filter[.!isnan.(pwwgt_filter)])

# ╔═╡ Cell order:
# ╠═eb55ab76-4877-11ef-2539-2549bd9d61aa
# ╠═77f31220-ef16-4b4c-aa82-a89d603d23da
# ╠═255b4912-b4f3-420c-a9b0-5854da29038a
# ╠═ddf80924-937f-4f62-afbe-b33764a36db4
# ╠═3942c663-ed45-4ec8-92e2-568ca89ea013
# ╟─1d428bbb-8c71-4c45-a2b9-2c44a6e18afb
# ╠═9fe2ac23-2e92-4a3b-9436-7f2f2f6876b0
# ╠═745ce483-74f7-4fc8-9944-c89d3abe098a
# ╟─c2c8bc39-af9c-4e09-a1fc-ec4f274182bb
# ╟─a0f24335-cd66-42d1-92d4-1e63115f3af8
# ╟─ce3de4d4-f9b0-4e06-be86-68e1e1d1d1c4
# ╠═e240d50c-a2c9-4510-8847-83de036b7f28
# ╟─abfe3dff-e445-4426-b963-3776ba310101
# ╟─2f48892c-31b2-443d-b788-1fbc6e3371ac
# ╠═d76c7c13-b3e6-49bf-88fa-fa48989689ed
# ╠═4ffbd55d-93e8-4956-9489-a34fe3249fa5
# ╠═692ac40e-bd85-4c77-858f-a2415731afac
# ╠═743ca220-78de-4878-98b9-785611bd56f7
