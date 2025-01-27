using Dates
using GeoRegions
using RegionGrids
using NCDatasets
using Printf
using Statistics
using Trapz

include(srcdir("backend.jl"))

function wrfwdecompose(
    geo  :: GeoRegion;
    start :: Date,
    stop  :: Date,
	days  :: Int = 0,
    nc    :: Int = 2
)

    dtbegstr = Dates.format(start,dateformat"yyyymmdd")
    dtbegend = Dates.format(stop,dateformat"yyyymmdd")
    timestr = "$(dtbegstr)_$(dtbegend)"
    smthstr = "smooth_$(@sprintf("%02d",days))days"

    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-p_wwgt-daily-$timestr-$smthstr.nc")
    end

    ds = NCDataset(fnc)
    ω  = ds["W"][:,:]
    p  = ds["P"][:,:] ./ 100
    close(ds)

    nt = size(p,2)
    c  = zeros(nt,2)

    for it = 1 : nt

        iω = @views ω[:,it]
        ip = @views p[:,it]; ps = ip[end]
        iiω = @views iω[ip.>=100]
        iip = @views ip[ip.>=100]; np = length(iip)

        if np > 1
            for ic = 1 : nc

                c[it,1] += iiω[1] * sin(ic*pi*(iip[1]-100)/(ps-100)) * (iip[1]-100)

                for jp = 2 : np

                    c[it,ic] += iiω[jp] * sin(ic*pi*(iip[jp]-100)/(ps-100)) *
                                          (iip[jp]-iip[jp-1])

                end

                c[it,ic] = c[it,ic] / (ps-100)

            end
        end

    end

    mkpath(datadir("wrf3","processed"))
    if iszero(days)
        fnc = datadir("wrf3","processed","$(geo.ID)-wcoeff-daily-$timestr.nc")
    else
        fnc = datadir("wrf3","processed","$(geo.ID)-wcoeff-daily-$timestr-$smthstr.nc")
    end

    if isfile(fnc); rm(fnc,force=true) end

    ds = NCDataset(fnc,"c")
    ds.dim["date"]  = nt
    ds.dim["coeff"] = nc

    nct = defVar(ds,"time",Int32,("date",),attrib=Dict(
        "units"     => "days since $(start) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian"
    ))

    ncc = defVar(ds,"wcoeff",Float32,("date","coeff"),attrib=Dict(
        "long_name" => "vertical_mode_coefficients",
        "full_name" => "Vertical Mode Coefficients of Pressure Velocity",
        "units"     => "m s**-1",
    ))

    nct.var[:] = collect(0 : (nt-1))
    ncc[:] = c

    close(ds)

end