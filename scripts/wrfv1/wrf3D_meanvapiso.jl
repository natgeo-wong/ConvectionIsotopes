using DrWatson
@quickactivate "ConvectionIsotopes"

using Dates
using Glob
using Logging
using NCDatasets

@info "$(now()) - ConvectionIsotopes - Initializing scripts to perform monthly mean on 3D WRF data ..."

gds  = NCDataset(datadir("wrf","wrfout","AUGp02","wrfout_d02_2019-08-11_00:00:00"))
lon = gds["XLONG"][:,:,1]
lat = gds["XLAT"][:,:,1]
close(gds)

for imo = 8 : 12

    @info "$(now()) - ConvectionIsotopes - Initializing temporary arrays for 3D data ..."

    v3D = zeros(Float32,543,543,50,8)
    
    v3D_qvap = zeros(Float32,543,543,50)
    v3D_δ18O = zeros(Float32,543,543,50)
    v3D_δHDO = zeros(Float32,543,543,50)

    for ifol = 1 : 3
        dsvec = glob(
            "wrfout3D_*",
            datadir("wrf","wrfout","$(uppercase(monthabbr(imo)))p0$ifol")
        )
        for inc = 2 : (length(dsvec) - 1)
            fnc = dsvec[inc]
            @info "$(now()) - ConvectionIsotopes - Extracting 3D data from $(fnc) ..."
            d3D = NCDataset(fnc)
            
            NCDatasets.load!(d3D["QVAPOR"].var,v3D,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 50, iy = 1 : 543, ix = 1 : 543
                v3D_qvap[ix,iy,ilvl] += v3D[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["HDO_QVAPOR"].var,v3D,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 50, iy = 1 : 543, ix = 1 : 543
                v3D_δHDO[ix,iy,ilvl] += v3D[ix,iy,ilvl,it]
            end
            
            NCDatasets.load!(d3D["O18_QVAPOR"].var,v3D,:,:,:,:)
            for it = 1 : 8, ilvl = 1 : 50, iy = 1 : 543, ix = 1 : 543
                v3D_δ18O[ix,iy,ilvl] += v3D[ix,iy,ilvl,it]
            end
            
            close(d3D)
        end
    end

    v3D_qvap = v3D_qvap / daysinmonth(Date(2019,imo)) / 8
    v3D_δHDO = v3D_δHDO / daysinmonth(Date(2019,imo)) / 8
    v3D_δ18O = v3D_δ18O / daysinmonth(Date(2019,imo)) / 8

    f3D = datadir("wrf","3D","$(uppercase(monthabbr(imo)))-QVAPOR_ISO")
    if isfile(f3D); rm(f3D,force=true) end
    d3D = NCDataset(f3D,"c")

    @info "$(now()) - ConvectionIsotopes - Saving monthly mean 3D data into $(f3D) ..."

    defDim(d3D,"longitude",543)
    defDim(d3D,"latitude",543)
    defDim(d3D,"level",50)

    nclon_3D = defVar(d3D,"longitude",Float32,("longitude","latitude",))
    nclat_3D = defVar(d3D,"latitude",Float32,("longitude","latitude",))
    ncqvap = defVar(d3D,"QVAPOR",Float32,("longitude","latitude","level",))
    ncδHDO = defVar(d3D,"HDO_QVAPOR",Float32,("longitude","latitude","level",))
    ncδ18O = defVar(d3D,"O18_QVAPOR",Float32,("longitude","latitude","level",))

    nclon_3D[:] = lon
    nclat_3D[:] = lat
    ncqvap[:] = v3D_qvap
    ncδHDO[:] = v3D_δHDO
    ncδ18O[:] = v3D_δ18O

    close(d3D)

end