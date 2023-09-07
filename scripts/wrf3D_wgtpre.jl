using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using Logging
using NCDatasets
using Trapz

@info "$(now()) - ConvectionIsotopes - Initializing scripts to perform monthly mean on 3D WRF data ..."

v3D_p = zeros(Float32,543,543,50)
v3D_z = zeros(Float32,543,543,51)
v3D_w = zeros(Float32,543,543,51)

tmp_z = zeros(Float32,50)
tmp_w = zeros(Float32,50)
tmpwp = zeros(Float32,50)
w_pre = zeros(Float32,543,543)

for imo = 8 : 12

    @info "$(now()) - ConvectionIsotopes - Initializing temporary arrays for 3D data ..."

    f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-W")
    d3D = NCDataset(f3D)

    lon = d3D["longitude"][:]; nlon = size(lon,1)
    lat = d3D["latitude"][:];  nlat = size(lat,2)
    NCDatasets.load!(d3D["W"].var,v3D_w,:,:,:)
    NCDatasets.load!(d3D["Z"].var,v3D_z,:,:,:)
    NCDatasets.load!(d3D["P"].var,v3D_p,:,:,:)

    close(d3D)

    for ilat = 1 : nlat, ilon = 1 : nlon

        for iz = 1 : 50
            tmp_z[iz] = (v3D_z[ilon,ilat,iz] + v3D_z[ilon,ilat,iz+1]) / 2
            tmp_w[iz] = (v3D_w[ilon,ilat,iz] + v3D_w[ilon,ilat,iz+1]) / 2
        end

        for iz = 2 : 49
            tmp_w[iz] = tmp_w[iz] *
                    (v3D_p[ilon,ilat,iz+1] - v3D_p[ilon,ilat,iz-1]) /
                    ((tmp_z[iz+1] - tmp_z[iz-1]) / 9.81)
            tmpwp[iz] = tmp_w[iz] * v3D_p[ilon,ilat,iz]
        end

        ptmp = @view v3D_p[ilon,ilat,:]
        wp = trapz(ptmp,tmpwp) / trapz(ptmp,tmp_w)
        if (wp > ptmp[1]) || (wp < ptmp[end])
              w_pre[ilon,ilat] = NaN
        else; w_pre[ilon,ilat] = wp
        end
        
    end

    @info "$(now()) - ConvectionIsotopes - Saving WRF W-weighted pressure"

    f2D = datadir("wrf","2D","$(uppercase(monthabbr(imo)))-p_wwgt")
    if isfile(f2D); rm(f2D,force=true) end
    d2D = NCDataset(f2D,"c")

    d2D.dim["longitude"] = size(lon,1)
    d2D.dim["latitude"]  = size(lat,2)

    nclon = defVar(d2D,"longitude",Float32,("longitude","latitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(d2D,"latitude",Float32,("longitude","latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncvar = defVar(d2D,"p_wwgt",Float32,("longitude","latitude",),attrib = Dict(
        "long_name"     => "column_mean_lagrangian_tendency_of_air_pressure",
        "full_name"     => "Vertical Wind Weighted Column Pressure",
        "units"         => "Pa",
    ))

    nclon[:] = lon
    nclat[:] = lat
    ncvar[:] = w_pre

    close(d2D)

end