using ERA5Reanalysis
using Logging
using NCDatasets
using Trapz

include(srcdir("backend.jl"))

function create_wp(
    e5ds :: ERA5Monthly,
    ereg :: ERA5Region
)
    
    @info "$(now()) - ColombiaIsotope - Preliminary preparation for calculation of vertical velocity weighted column-mean pressure  ..."

    dtbeg = e5ds.dtbeg
    dtend = e5ds.dtend
    dtvec = dtbeg : Year(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)
    mask = lsd.mask

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)

    @info "$(now()) - ColombiaIsotope - Preallocating arrays ..."

    wp = Array{Float32,3}(undef,nlon,nlat,12)
    sp = Array{Float32,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    wa = Array{Float32,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wpii = Vector{Float32}(undef,np+2)
    waii = Vector{Float32}(undef,np+2)
    
    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ColombiaIsotope - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_wp = SingleVariable("p_wwgt");
    evar_sp = SingleVariable("sp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    plvl = Float32.(vcat(0,plvl,0))

    for dtii in dtvec

        @info "$(now()) - ColombiaIsotope - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        disable_logging(Logging.Warn)
        sp_ds = read(e5ds,evar_sp,ereg,dtii)
        sp_sc = sp_ds[evar_sp.varID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.varID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.varID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.varID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["_FillValue"]
        end
        disable_logging(Logging.Debug)

        @info "$(now()) - ColombiaIsotope - Calculating the $(evar_wp.vname) data over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        for it = 1 : 12

            NCDatasets.load!(sp_ds[evar_sp.varID].var,sptmp,:,:,it)
            int2real!(sp,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            for ip = 1 : np
                waip = @view wa[:,:,ip]
                NCDatasets.load!(wa_ds[ip][evar_wa[ip].varID].var,watmp,:,:,it)
                int2real!(
                    waip,watmp,
                    scale=wa_sc[ip],offset=wa_of[ip],
                    mvalue=wa_mv[ip],fvalue=wa_fv[ip]
                )
            end

            for ilat = 1 : nlat, ilon = 1 : nlon

                if isone(mask[ilon,ilat])
                    spii = sp[ilon,ilat] / 100
                    plvl[end] = spii
                    for ip = 1 : np
                        waii[ip+1] = abs(wa[ilon,ilat,ip])
                        wpii[ip+1] = waii[ip+1] * plvl[ip+1]
                    end

                    pre = @view plvl[plvl.<=spii]
                    wan = @view waii[plvl.<=spii]
                    wpn = @view wpii[plvl.<=spii]

                    if pre[end-1] != pre[end]
                        inttop = integrate(pre,wpn)
                        intbot = integrate(pre,wan)
                    else
                        pre2 = @view pre[1:(end-1)]
                        wan2 = @view wan[1:(end-1)]
                        wpn2 = @view wpn[1:(end-1)]
                        inttop = trapz(pre2,wpn2)
                        intbot = trapz(pre2,wan2)
                    end

                    wp[ilon,ilat,it] = inttop / intbot * 100
                else
                    wp[ilon,ilat,it] = NaN
                end

            end

            @info "$(now()) - ColombiaIsotope - Extraction, Calculation and Integration completed over Month $it of $(year(dtii)) ..."

        end

        ERA5Reanalysis.save(wp,dtii,e5ds,evar_wp,ereg,lsd)

    end

end