using ERA5Reanalysis
using Logging
using NCDatasets
using Trapz

include(srcdir("backend.jl"))

function create_wp(
    e5ds :: ERA5Daily,
    ereg :: ERA5Region;
    smooth :: Bool = false,
    smoothtime :: Int = 30,
)
    
    @info "$(now()) - ConvectionIsotopes - Preliminary preparation for calculation of vertical velocity weighted column-mean pressure  ..."

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Month(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    wp = Array{Float64,3}(undef,nlon,nlat,31)
    sp = Array{Float64,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    tp = Array{Float64,2}(undef,nlon,nlat);    tptmp = Array{Int16,2}(undef,nlon,nlat)
    wa = Array{Float64,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wpii = Vector{Float64}(undef,np+2)
    waii = Vector{Float64}(undef,np+2)
    
    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_wp = SingleVariable("p_wwgt",path=srcdir());
    evar_sp = SingleVariable("sp");
    evar_tp = SingleVariable("tp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    plvl = Float32.(vcat(0,plvl,0))

    for dtii in dtvec

        ndy = daysinmonth(dtii)

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        disable_logging(Logging.Warn)
        sp_ds = read(e5ds,evar_sp,ereg,dtii,smooth=smooth,smoothtime=smoothtime)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]

        tp_ds = read(e5ds,evar_tp,ereg,dtii,smooth=smooth,smoothtime=smoothtime)
        tp_sc = tp_ds[evar_tp.ID].attrib["scale_factor"]
        tp_of = tp_ds[evar_tp.ID].attrib["add_offset"]
        tp_mv = tp_ds[evar_tp.ID].attrib["missing_value"]
        tp_fv = tp_ds[evar_tp.ID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii,smooth=smooth,smoothtime=smoothtime)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["_FillValue"]
        end
        disable_logging(Logging.Debug)

        @info "$(now()) - ConvectionIsotopes - Calculating the $(evar_wp.name) data over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        for it = 1 : ndy

            NCDatasets.load!(sp_ds[evar_sp.ID].var,sptmp,:,:,it)
            int2real!(sp,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            NCDatasets.load!(tp_ds[evar_tp.ID].var,tptmp,:,:,it)
            int2real!(tp,tptmp,scale=tp_sc,offset=tp_of,mvalue=tp_mv,fvalue=tp_fv)

            for ip = 1 : np
                waip = @view wa[:,:,ip]
                NCDatasets.load!(wa_ds[ip][evar_wa[ip].ID].var,watmp,:,:,it)
                int2real!(
                    waip,watmp,
                    scale=wa_sc[ip],offset=wa_of[ip],
                    mvalue=wa_mv[ip],fvalue=wa_fv[ip]
                )
            end

            for ilat = 1 : nlat, ilon = 1 : nlon

                if !isnan(lsd.z[ilon,ilat]) && (tp[ilon,ilat] > 0.005/24)
                    spii = sp[ilon,ilat] / 100
                    plvl[end] = spii
                    for ip = 1 : np
                        waii[ip+1] = wa[ilon,ilat,ip]
                        wpii[ip+1] = waii[ip+1] * plvl[ip+1]
                    end

                    pre = @view plvl[plvl.<=spii]
                    wan = @view waii[plvl.<=spii]
                    wpn = @view wpii[plvl.<=spii]

                    if pre[end-1] != pre[end]
                        inttop = trapz(pre,wpn)
                        intbot = trapz(pre,wan)
                    else
                        pre2 = @view pre[1:(end-1)]
                        wan2 = @view wan[1:(end-1)]
                        wpn2 = @view wpn[1:(end-1)]
                        inttop = trapz(pre2,wpn2)
                        intbot = trapz(pre2,wan2)
                    end

                    wp[ilon,ilat,it] = inttop / intbot * 100
                    if (wp[ilon,ilat,it] > sp[ilon,ilat]) || (wp[ilon,ilat,it] < 0)
                        wp[ilon,ilat,it] = NaN
                    end
                else
                    wp[ilon,ilat,it] = NaN
                end

            end

        end

        @info "$(now()) - ConvectionIsotopes - Extraction, Calculation and Integration completed over Month $(dtii) ..."

        ERA5Reanalysis.save(
            view(wp,:,:,1:ndy),dtii,e5ds,evar_wp,ereg,lsd,
            smooth=smooth,smoothtime=smoothtime
        )

    end

end

function create_wp(
    e5ds :: ERA5Monthly,
    ereg :: ERA5Region
)
    
    @info "$(now()) - ConvectionIsotopes - Preliminary preparation for calculation of vertical velocity weighted column-mean pressure  ..."

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Year(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    wp = Array{Float32,3}(undef,nlon,nlat,12)
    sp = Array{Float64,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    wa = Array{Float64,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wpii = Vector{Float64}(undef,np+2)
    waii = Vector{Float64}(undef,np+2)
    
    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_wp = SingleVariable("p_wwgt",path=srcdir());
    evar_sp = SingleVariable("sp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    plvl = Float32.(vcat(0,plvl,0))

    for dtii in dtvec

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        disable_logging(Logging.Warn)
        sp_ds = read(e5ds,evar_sp,ereg,dtii)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["_FillValue"]
        end
        disable_logging(Logging.Debug)

        @info "$(now()) - ConvectionIsotopes - Calculating the $(evar_wp.name) data over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        for it = 1 : 12

            NCDatasets.load!(sp_ds[evar_sp.ID].var,sptmp,:,:,it)
            int2real!(sp,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            for ip = 1 : np
                waip = @view wa[:,:,ip]
                NCDatasets.load!(wa_ds[ip][evar_wa[ip].ID].var,watmp,:,:,it)
                int2real!(
                    waip,watmp,
                    scale=wa_sc[ip],offset=wa_of[ip],
                    mvalue=wa_mv[ip],fvalue=wa_fv[ip]
                )
            end

            for ilat = 1 : nlat, ilon = 1 : nlon

                if !isnan(lsd.z[ilon,ilat])
                    spii = sp[ilon,ilat] / 100
                    plvl[end] = spii
                    for ip = 1 : np
                        waii[ip+1] = wa[ilon,ilat,ip]
                        wpii[ip+1] = waii[ip+1] * plvl[ip+1]
                    end

                    pre = @view plvl[plvl.<=spii]
                    wan = @view waii[plvl.<=spii]
                    wpn = @view wpii[plvl.<=spii]

                    if pre[end-1] != pre[end]
                        inttop = trapz(pre,wpn)
                        intbot = trapz(pre,wan)
                    else
                        pre2 = @view pre[1:(end-1)]
                        wan2 = @view wan[1:(end-1)]
                        wpn2 = @view wpn[1:(end-1)]
                        inttop = trapz(pre2,wpn2)
                        intbot = trapz(pre2,wan2)
                    end

                    wp[ilon,ilat,it] = inttop / intbot * 100
                    if (wp[ilon,ilat,it] > sp[ilon,ilat]) || (wp[ilon,ilat,it] < 0)
                        wp[ilon,ilat,it] = NaN
                    end
                else
                    wp[ilon,ilat,it] = NaN
                end

            end

            @info "$(now()) - ConvectionIsotopes - Extraction, Calculation and Integration completed over Month $it of $(year(dtii)) ..."

        end

        ERA5Reanalysis.save(wp,dtii,e5ds,evar_wp,ereg,lsd)

    end

end

function compiled_wp(
    e5ds :: ERA5Daily,
    ereg :: ERA5Region
)
    
    @info "$(now()) - ConvectionIsotopes - Preliminary preparation for calculation of vertical velocity weighted column-mean pressure  ..."

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Month(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    sp_it = Array{Float64,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    tp_it = Array{Float64,2}(undef,nlon,nlat);    tptmp = Array{Int16,2}(undef,nlon,nlat)
    wa_it = Array{Float64,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wp = Array{Float64,2}(undef,nlon,nlat)
    sp = Array{Float64,2}(undef,nlon,nlat)
    tp = Array{Float64,2}(undef,nlon,nlat)
    wa = Array{Float64,3}(undef,nlon,nlat,np)

    wpii = Vector{Float64}(undef,np+2)
    waii = Vector{Float64}(undef,np+2)
    
    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_sp = SingleVariable("sp");
    evar_tp = SingleVariable("tp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    plvl = Float32.(vcat(0,plvl,0))

    for dtii in dtvec

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) $(monthname(dtii)) ..."

        ndy = daysinmonth(dtii)

        sp_ds = read(e5ds,evar_sp,ereg,dtii,quiet=true)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]

        tp_ds = read(e5ds,evar_tp,ereg,dtii,quiet=true)
        tp_sc = tp_ds[evar_tp.ID].attrib["scale_factor"]
        tp_of = tp_ds[evar_tp.ID].attrib["add_offset"]
        tp_mv = tp_ds[evar_tp.ID].attrib["missing_value"]
        tp_fv = tp_ds[evar_tp.ID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii,quiet=true)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["_FillValue"]
        end

        for it = 1 : ndy

            NCDatasets.load!(sp_ds[evar_sp.ID].var,sptmp,:,:,it)
            int2real!(sp_it,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            NCDatasets.load!(tp_ds[evar_tp.ID].var,tptmp,:,:,it)
            int2real!(tp_it,tptmp,scale=tp_sc,offset=tp_of,mvalue=tp_mv,fvalue=tp_fv)

            for ip = 1 : np
                waip = @view wa_it[:,:,ip]
                NCDatasets.load!(wa_ds[ip][evar_wa[ip].ID].var,watmp,:,:,it)
                int2real!(
                    waip,watmp,
                    scale=wa_sc[ip],offset=wa_of[ip],
                    mvalue=wa_mv[ip],fvalue=wa_fv[ip]
                )
            end

            sp .+= sp_it
            tp .+= tp_it
            wa .+= wa_it

        end

    end

    sp ./= (Dates.value(dtend - dtbeg) + 1)
    tp ./= (Dates.value(dtend - dtbeg) + 1)
    wa ./= (Dates.value(dtend - dtbeg) + 1)

    for ilat = 1 : nlat, ilon = 1 : nlon

        if !isnan(lsd.z[ilon,ilat]) && (tp[ilon,ilat] > 0.005)
            spii = sp[ilon,ilat] / 100
            plvl[end] = spii
            for ip = 1 : np
                waii[ip+1] = wa[ilon,ilat,ip]
                wpii[ip+1] = waii[ip+1] * plvl[ip+1]
            end

            pre = @view plvl[plvl.<=spii]
            wan = @view waii[plvl.<=spii]
            wpn = @view wpii[plvl.<=spii]

            if pre[end-1] != pre[end]
                inttop = trapz(pre,wpn)
                intbot = trapz(pre,wan)
            else
                pre2 = @view pre[1:(end-1)]
                wan2 = @view wan[1:(end-1)]
                wpn2 = @view wpn[1:(end-1)]
                inttop = trapz(pre2,wpn2)
                intbot = trapz(pre2,wan2)
            end

            wp[ilon,ilat] = inttop / intbot * 100
            if (wp[ilon,ilat] > sp[ilon,ilat]) || (wp[ilon,ilat] < 0)
                wp[ilon,ilat] = NaN
            end
        else
            wp[ilon,ilat] = NaN
        end

    end

    wσ = wp ./sp

    fnc = datadir("$(ereg.string)-p_wwgt-compiled-$(Dates.format(e5ds.start,dateformat"yyyymmdd"))_$(Dates.format(e5ds.stop,dateformat"yyyymmdd")).nc")
    if isfile(fnc); rm(fnc,force=true) end
    ds  = NCDataset(fnc,"c")

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
	    "units"     => "degrees_east",
	    "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
	    "units"     => "degrees_north",
	    "long_name" => "latitude",
	))

	nctp = defVar(ds,evar_tp.ID,Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => evar_tp.units,
	    "long_name" => evar_tp.long,
		"full_name" => evar_tp.name,
	))

	ncsp = defVar(ds,evar_sp.ID,Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => evar_sp.units,
	    "long_name" => evar_sp.long,
		"full_name" => evar_sp.name,
	))

	ncpω = defVar(ds,"p_wwgt",Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => "Pa",
	    "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
		"full_name" => "Vertical Wind Weighted Column Pressure",
	))

	ncσω = defVar(ds,"σ_wwgt",Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => "0-1",
	    "long_name" => "column_mean_lagrangian_tendency_of_sigma",
		"full_name" => "Vertical Wind Weighted Sigma",
	))

    nclon[:]  = lsd.lon
    nclat[:]  = lsd.lat
    nctp[:,:] = tp
    ncsp[:,:] = sp
    ncpω[:,:] = wp
    ncσω[:,:] = wσ

    close(ds)

end

function climatology_wp(
    e5ds :: ERA5Monthly,
    ereg :: ERA5Region
)
    
    @info "$(now()) - ConvectionIsotopes - Preliminary preparation for calculation of vertical velocity weighted column-mean pressure  ..."

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Year(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    sp_it = Array{Float64,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    tp_it = Array{Float64,2}(undef,nlon,nlat);    tptmp = Array{Int16,2}(undef,nlon,nlat)
    wa_it = Array{Float64,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wp = Array{Float64,2}(undef,nlon,nlat)
    sp = Array{Float64,2}(undef,nlon,nlat)
    tp = Array{Float64,2}(undef,nlon,nlat)
    wa = Array{Float64,3}(undef,nlon,nlat,np)

    wpii = Vector{Float64}(undef,np+2)
    waii = Vector{Float64}(undef,np+2)
    
    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_sp = SingleVariable("sp");
    evar_tp = SingleVariable("tp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    plvl = Float32.(vcat(0,plvl,0))

    for dtii in dtvec

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        sp_ds = read(e5ds,evar_sp,ereg,dtii,quiet=true)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]

        tp_ds = read(e5ds,evar_tp,ereg,dtii,quiet=true)
        tp_sc = tp_ds[evar_tp.ID].attrib["scale_factor"]
        tp_of = tp_ds[evar_tp.ID].attrib["add_offset"]
        tp_mv = tp_ds[evar_tp.ID].attrib["missing_value"]
        tp_fv = tp_ds[evar_tp.ID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii,quiet=true)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].ID].attrib["_FillValue"]
        end

        for it = 1 : 12

            NCDatasets.load!(sp_ds[evar_sp.ID].var,sptmp,:,:,it)
            int2real!(sp_it,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            NCDatasets.load!(tp_ds[evar_tp.ID].var,tptmp,:,:,it)
            int2real!(tp_it,tptmp,scale=tp_sc,offset=tp_of,mvalue=tp_mv,fvalue=tp_fv)

            for ip = 1 : np
                waip = @view wa_it[:,:,ip]
                NCDatasets.load!(wa_ds[ip][evar_wa[ip].ID].var,watmp,:,:,it)
                int2real!(
                    waip,watmp,
                    scale=wa_sc[ip],offset=wa_of[ip],
                    mvalue=wa_mv[ip],fvalue=wa_fv[ip]
                )
            end

            sp .+= sp_it
            tp .+= tp_it
            wa .+= wa_it

        end

    end

    sp ./= (12 * length(dtvec))
    tp ./= (12 * length(dtvec))
    wa ./= (12 * length(dtvec))

    for ilat = 1 : nlat, ilon = 1 : nlon

        if !isnan(lsd.z[ilon,ilat]) && (tp[ilon,ilat] > 0.005)
            spii = sp[ilon,ilat] / 100
            plvl[end] = spii
            for ip = 1 : np
                waii[ip+1] = wa[ilon,ilat,ip]
                wpii[ip+1] = waii[ip+1] * plvl[ip+1]
            end

            pre = @view plvl[plvl.<=spii]
            wan = @view waii[plvl.<=spii]
            wpn = @view wpii[plvl.<=spii]

            if pre[end-1] != pre[end]
                inttop = trapz(pre,wpn)
                intbot = trapz(pre,wan)
            else
                pre2 = @view pre[1:(end-1)]
                wan2 = @view wan[1:(end-1)]
                wpn2 = @view wpn[1:(end-1)]
                inttop = trapz(pre2,wpn2)
                intbot = trapz(pre2,wan2)
            end

            wp[ilon,ilat] = inttop / intbot * 100
            if (wp[ilon,ilat] > sp[ilon,ilat]) || (wp[ilon,ilat] < 0)
                wp[ilon,ilat] = NaN
            end
        else
            wp[ilon,ilat] = NaN
        end

    end

    wσ = wp ./sp

    fnc = datadir("$(ereg.string)-p_wwgt-climatology-$(year(e5ds.start))_$(year(e5ds.stop)).nc")
    if isfile(fnc); rm(fnc,force=true) end
    ds  = NCDataset(fnc,"c")

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
	    "units"     => "degrees_east",
	    "long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
	    "units"     => "degrees_north",
	    "long_name" => "latitude",
	))

	nctp = defVar(ds,evar_tp.ID,Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => evar_tp.units,
	    "long_name" => evar_tp.long,
		"full_name" => evar_tp.name,
	))

	ncsp = defVar(ds,evar_sp.ID,Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => evar_sp.units,
	    "long_name" => evar_sp.long,
		"full_name" => evar_sp.name,
	))

	ncpω = defVar(ds,"p_wwgt",Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => "Pa",
	    "long_name" => "column_mean_lagrangian_tendency_of_air_pressure",
		"full_name" => "Vertical Wind Weighted Column Pressure",
	))

	ncσω = defVar(ds,"σ_wwgt",Float64,("longitude","latitude"),attrib = Dict(
	    "units"     => "0-1",
	    "long_name" => "column_mean_lagrangian_tendency_of_sigma",
		"full_name" => "Vertical Wind Weighted Sigma",
	))

    nclon[:]  = lsd.lon
    nclat[:]  = lsd.lat
    nctp[:,:] = tp
    ncsp[:,:] = sp
    ncpω[:,:] = wp
    ncσω[:,:] = wσ

    close(ds)

end