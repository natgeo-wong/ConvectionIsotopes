using DelimitedFiles
using GeoRegions
using Statistics

function stninfoload()

    return readdlm(srcdir("stninfo.csv"),',',header=true)

end

function stninfoall()

    return stninfoload()[1][:,1:3]

end

function stninfody()

    stninfo = stninfoload()[1]
    tf = stninfo[:,4]
    return stninfo[tf.==true,1:3]

end

function stninfomo()

    stninfo = stninfoload()[1]
    tf = stninfo[:,5]
    return stninfo[tf.==true,1:3]

end

function stninfocostarica()

    stninfo = stninfoload()[1]
    tf = stninfo[:,8]
    return stninfo[tf.==true,1:3]

end

function extractregclimate(
    data_yr :: Array{<:Real,2},
    data_mo :: Array{<:Real,3},
    geo :: GeoRegion,
    lon :: Vector,
    lat :: Vector
)

	ggrd = RegionGrid(geo,lon,lat)
    nlon = length(ggrd.ilon)
    nlat = length(ggrd.ilat)
    
	data_yr_tmp = zeros(nlon,nlat)
	data_mo_tmp = zeros(nlon,nlat,12)
	data_mo_new = zeros(12)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for glat in 1 : nlat, glon in 1 : nlon
		data_yr_tmp[glon,glat] = data_yr[
			ggrd.ilon[glon],ggrd.ilat[glat]
		] * mask[glon,glat]
	end

	for imo = 1 : 12
		for glat in 1 : nlat, glon in 1 : nlon
			data_mo_tmp[glon,glat,imo] = data_mo[
				ggrd.ilon[glon],ggrd.ilat[glat],imo
			] * mask[glon,glat]
		end
	end

	data_yr_new = mean(data_yr_tmp[.!isnan.(data_yr_tmp)])
	for imo = 1 : 12
		data_mo_ii   = @view data_mo_tmp[:,:,imo]
		data_mo_new[imo] = mean(data_mo_ii[.!isnan.(data_mo_ii)])
	end

	return data_yr_new,data_mo_new

end

SMOW18O(x :: Real = 1.) = x * 2.0052e-3
SMOWHDO(x :: Real = 1.) = x * 1.5576e-4

function SMOW(x::Real, iso::AbstractString)
	if iso == "HDO_"
		return SMOWHDO(x)
	elseif iso == "18O_"
		return SMOW18O(x)
	end
end