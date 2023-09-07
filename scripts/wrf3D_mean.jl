using DrWatson
@quickactivate "ColombiaIsotope"

using Dates
using Glob
using Logging
using NCDatasets

@info "$(now()) - ColombiaIsotope - Initializing scripts to perform monthly mean on 3D WRF data ..."

gds  = NCDataset(datadir("wrf","AUGp02","wrfout_d02_2019-08-11_00:00:00"))
lon = gds["XLONG"][:,:,1]
lat = gds["XLAT"][:,:,1]
close(gds)

for imo = 8 : 12

    @info "$(now()) - ColombiaIsotope - Initializing temporary arrays for 3D data ..."

    v3D_1 = zeros(Float32,543,543,50,8)
    v3D_2 = zeros(Float32,543,543,51,8)
    
    v3D_p = zeros(Float32,543,543,50)
    v3D_z = zeros(Float32,543,543,51)
    v3D_w = zeros(Float32,543,543,51)

    for ifol = 1 : 3
        dsvec = glob(
            "wrfout3D_*",
            datadir("wrf","wrfout","$(uppercase(monthabbr(imo)))p0$ifol")
        )
        for inc = 2 : (length(dsvec) - 1)
            fnc = dsvec[inc]
            @info "$(now()) - ColombiaIsotope - Extracting 3D data from $(fnc) ..."
            d3D = NCDataset(fnc)
            
            NCDatasets.load!(d3D["P"].var,v3D_1,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 50, iy = 1 : 543, ix = 1 : 543
                v3D_p[ix,iy,ilvl] += v3D_1[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["PB"].var,v3D_1,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 50, iy = 1 : 543, ix = 1 : 543
                v3D_p[ix,iy,ilvl] += v3D_1[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["PH"].var,v3D_2,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 51, iy = 1 : 543, ix = 1 : 543
                v3D_z[ix,iy,ilvl] += v3D_2[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["PHB"].var,v3D_2,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 51, iy = 1 : 543, ix = 1 : 543
                v3D_z[ix,iy,ilvl] += v3D_2[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["W"].var,v3D_2,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 51, iy = 1 : 543, ix = 1 : 543
                v3D_w[ix,iy,ilvl] += v3D_2[ix,iy,ilvl,it]
            end
            
            close(d3D)
        end
    end

    v3D_p = v3D_p / daysinmonth(Date(2019,imo)) / 8
    v3D_z = v3D_z / daysinmonth(Date(2019,imo)) / 8
    v3D_w = v3D_w / daysinmonth(Date(2019,imo)) / 8

    f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-W")
    if isfile(f3D); rm(f3D,force=true) end
    d3D = NCDataset(f3D,"c")

    @info "$(now()) - ColombiaIsotope - Saving monthly mean 3D data into $(f3D) ..."

    defDim(d3D,"longitude",543)
    defDim(d3D,"latitude",543)
    defDim(d3D,"level",50)
    defDim(d3D,"level_staggered",51)

    nclon_3D = defVar(d3D,"longitude",Float32,("longitude","latitude",))
    nclat_3D = defVar(d3D,"latitude",Float32,("longitude","latitude",))
    ncp_3D = defVar(d3D,"P",Float32,("longitude","latitude","level",))
    ncz_3D = defVar(d3D,"W",Float32,("longitude","latitude","level_staggered",))
    ncw_3D = defVar(d3D,"Z",Float32,("longitude","latitude","level_staggered",))

    nclon_3D[:] = lon
    nclat_3D[:] = lat
    ncp_3D[:] = v3D_p
    ncz_3D[:] = v3D_z
    ncw_3D[:] = v3D_w

    close(d3D)

end