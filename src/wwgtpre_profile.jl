using Dierckx
using ERA5Reanalysis
using Logging
using NCDatasets
using Statistics

include(srcdir("backend.jl"))

function wprofile(
    e5ds :: ERA5Monthly,
    ereg :: ERA5Region;
    σbin :: AbstractRange = 0 : 0.01 : 1,
)

    dtbeg = e5ds.dtbeg
    dtend = e5ds.dtend
    dtvec = dtbeg : Year(1) : dtend

    lsd  = getLandSea(e5ds,ereg);
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)
    mask = lsd.mask

    plvl = sort(era5Pressures())
    plvl = plvl[plvl.>=10]; np = length(plvl)
    ptmp = Vector{Float64}(undef,np+2)
    waii = Vector{Float64}(undef,np+2)
    ind  = ones(Bool,np+2)

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."
    
    wp = Array{Float32,2}(undef,nlon,nlat);    wptmp = Array{Int16,2}(undef,nlon,nlat)
    sp = Array{Float32,2}(undef,nlon,nlat);    sptmp = Array{Int16,2}(undef,nlon,nlat)
    wa = Array{Float32,3}(undef,nlon,nlat,np); watmp = Array{Int16,2}(undef,nlon,nlat)

    wa_ds = Vector{Any}(undef,np)
    wa_sc = Vector{Float64}(undef,np)
    wa_of = Vector{Float64}(undef,np)
    wa_mv = Vector{Int16}(undef,np)
    wa_fv = Vector{Int16}(undef,np)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_wp = SingleVariable("p_wwgt");
    evar_sp = SingleVariable("sp");
    evar_wa = Vector{PressureVariable}(undef,np)
    for ip = 1 : np
        evar_wa[ip] = PressureVariable("w",hPa=plvl[ip]);
    end
    disable_logging(Logging.Debug)

    @info "$(now()) - ConvectionIsotopes - Binning the $(evar_wp.vname) data over the $(ereg.geo.name) Region ..."
    σlvls = vcat(0:0.01:1); nσ = length(σlvls)
    waσ   = zeros(nσ)
    pfreq = zeros(Int,length(σbin)-1)
    plvlb = zeros( nσ,length(σbin)-1)

    ds  = NCDataset(datadir("flsm","flsm-$(ereg.geoID).nc"))
    lsm = ds["flsm"][:]
    close(ds)

    for dtii in dtvec

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure, Weighted Column Pressure, and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        disable_logging(Logging.Warn)
        sp_ds = read(e5ds,evar_sp,ereg,dtii)
        sp_sc = sp_ds[evar_sp.varID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.varID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.varID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.varID].attrib["_FillValue"]

        wp_ds = read(e5ds,evar_wp,ereg,dtii)
        wp_sc = wp_ds[evar_wp.varID].attrib["scale_factor"]
        wp_of = wp_ds[evar_wp.varID].attrib["add_offset"]
        wp_mv = wp_ds[evar_wp.varID].attrib["missing_value"]
        wp_fv = wp_ds[evar_wp.varID].attrib["_FillValue"]

        for ip = 1 : np
            wa_ds[ip] = read(e5ds,evar_wa[ip],ereg,dtii)
            wa_sc[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["scale_factor"]
            wa_of[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["add_offset"]
            wa_mv[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["missing_value"]
            wa_fv[ip] = wa_ds[ip][evar_wa[ip].varID].attrib["_FillValue"]
        end
        disable_logging(Logging.Debug)

        for it = 1 : 12

            NCDatasets.load!(sp_ds[evar_sp.varID].var,sptmp,:,:,it)
            int2real!(sp,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

            NCDatasets.load!(wp_ds[evar_wp.varID].var,wptmp,:,:,it)
            int2real!(wp,wptmp,scale=wp_sc,offset=wp_of,mvalue=wp_mv,fvalue=wp_fv)

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

                if isone(mask[ilon,ilat]) && !isnan(wp[ilon,ilat]) && (lsm[ilon,ilat]<0.9)

                    spii = sp[ilon,ilat]
                    wpii = wp[ilon,ilat] / spii
                    for ip = 1 : np
                        waii[ip+1] = wa[ilon,ilat,ip]
                        ptmp[ip+1] = plvl[ip] / spii * 100
                    end
                    ptmp[end] = 1.

                    for ip = 1 : (np+1)
                        ind[ip] = ptmp[ip] < 1
                    end

                    pσt = @view ptmp[ind]
                    wat = @view waii[ind]
                    spl = Spline1D(pσt,wat)

                    for iσ = 2 : (nσ-1)
                        waσ[iσ] = evaluate(spl,σlvls[iσ])
                    end

                    if wpii > 1; wpii = 1 end
                    if wpii < 0; wpii = 0 end

                    binwprofile(plvlb,pfreq,wpii,waσ,σbin)

                end
                
            end

        end

    end

    savebin(plvlb,pfreq,σbin,σlvls,ereg)

end

function binwprofile(
    plvlb :: AbstractArray,
    pfreq :: AbstractVector,
    wpii :: Real,
    waσ  :: AbstractVector,
    σbin :: AbstractRange
)

    nσ = length(waσ); nb = length(σbin)
    jj = (wpii .<= view(σbin,2:nb)) .& (wpii .> view(σbin,1:(nb-1)))
    jj = findfirst(isone,jj)
    pfreq[jj] += 1
    for iσ = 2 : (nσ-1)
        plvlb[iσ,jj] += waσ[iσ]
    end
    
    return

end

function savebin(
    plvlb :: AbstractArray,
    pfreq :: AbstractVector,
    σbin :: AbstractVector,
    σlvl :: AbstractVector,
    ereg :: ERA5Region
)

    fnc = datadir("wprofile","wwgtpre-profile-$(ereg.gstr).nc")
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end
    if isfile(fnc)
        @info "$(now()) - ConvectionIsotopes - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) with ConvectionIsotopes scripts",
    ))

    ds.dim["σ_wgt"] = length(σbin)
    ds.dim["σ_bin"] = length(σbin) - 1
    ds.dim["level"] = length(σlvl)

    ncσ = defVar(ds,"σ",Float32,("level",),attrib = Dict(
        "long_name"     => "sigma_level",
        "full_name"     => "σ",
        "units"         => "0-1",
    ))

    ncσw = defVar(ds,"σ_wgt",Float32,("σ_wgt",),attrib = Dict(
        "long_name"     => "sigma_bin",
        "full_name"     => "σ Bin",
        "units"         => "0-1",
    ))

    ncpf = defVar(ds,"p_freq",Float32,("σ_bin",),attrib = Dict(
        "long_name"     => "sigma_bin_frequency",
        "full_name"     => "σ Bin Frequency"
    ))

    ncpl = defVar(ds,"w_profile",Float32,("level","σ_bin",),attrib = Dict(
        "long_name"     => "binned_vertical_profiles_of_vertical_wind",
        "full_name"     => "Binned Vertical Profiles of Vertical Wind",
        "units"         => "Pa s**-1"
    ))

    ncσ[:]  = σlvl
    ncσw[:] = σbin
    ncpf[:] = pfreq
    ncpl[:] = plvlb

    close(ds)

    return

end
