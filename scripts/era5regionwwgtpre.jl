using DrWatson
@quickactivate "ConvectionIsotopes"

using ERA5Reanalysis

include(srcdir("common.jl"))

e5ds = ERA5Daily(start=Date(2019,2),stop=Date(2021,6),path=datadir())
egeo = ERA5Region("OTREC")
evar = SingleVariable("p_wwgt")
elsd = getLandSea(e5ds,egeo)

dtvec = e5ds.start : Day(1) : e5ds.stop
ndt   = length(dtvec)
pω07  = zeros(ndt,12)
pω30  = zeros(ndt,12)
stndy = stninfody()
lon   = stndy[:,2]; ilon = zeros(Int,12)
lat   = stndy[:,3]; ilat = zeros(Int,12)

for istn = 1 : 12
    ilon[istn] = argmin(abs.(elsd.lon.-lon[istn]))
    ilat[istn] = argmin(abs.(elsd.lat.-lat[istn]))
end

for idt in e5ds.start : Month(1) : e5ds.stop

    it = (dtvec .>= Date(year(idt),month(idt))) .& (dtvec .<= Date(year(idt),month(idt),daysinmonth(idt)))

    ds = read(e5ds,evar,egeo,idt,smooth=true,smoothtime=7)
    for istn = 1 : 12
        pω07[it,istn] .= nomissing(ds[evar.ID][ilon[istn],ilat[istn],:],NaN)
    end
    close(ds)

    ds = read(e5ds,evar,egeo,idt,smooth=true,smoothtime=30)
    for istn = 1 : 12
        pω30[it,istn] .= nomissing(ds[evar.ID][ilon[istn],ilat[istn],:],NaN)
    end
    close(ds)

end

fnc = datadir("processed-p_wwgt.nc")
if isfile(fnc); rm(fnc,force=true) end
pds = NCDataset(fnc,"c")

pds.dim["station"] = 12
pds.dim["time"] = ndt

nclon = defVar(pds,"longitude",Float32,("station",),attrib = Dict(
    "units"     => "degrees_east",
    "long_name" => "longitude",
))

nclat = defVar(pds,"latitude",Float32,("station",),attrib = Dict(
    "units"     => "degrees_north",
    "long_name" => "latitude",
))

nctime = defVar(pds,"time",Int32,("time",),attrib = Dict(
    "units"     => "Days since $(dtvec[1]) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncpω07 = defVar(pds,"pω07",Float32,("time","station",),attrib = Dict(
    "units"     => "Pa",
    "long_name" => "7-Day Moving Average of ERA5 p_wwgt",
))

ncpω30 = defVar(pds,"pω30",Float32,("time","station",),attrib = Dict(
    "units"     => "Pa",
    "long_name" => "30-Day Moving Average of ERA5 p_wwgt",
))

nclon[:]   = lon
nclat[:]   = lat
nctime[:]  = 0 : (length(dtvec) - 1)
ncpω07[:,:]  = pω07
ncpω30[:,:]  = pω30

close(pds)