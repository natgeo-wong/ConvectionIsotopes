using ERA5Reanalysis
using Logging
using NCDatasets
using Trapz

include(srcdir("backend.jl"))

function create_wσ(
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

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    wp = Array{Float32,3}(undef,nlon,nlat,12); wptmp = Array{Int16,3}(undef,nlon,nlat,12)
    sp = Array{Float32,3}(undef,nlon,nlat,12); sptmp = Array{Int16,3}(undef,nlon,nlat,12)
    wσ = Array{Float32,3}(undef,nlon,nlat,12)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_sp = SingleVariable("sp")
    evar_wp = SingleVariable("p_wwgt")
    evar_wσ = SingleVariable("σ_wwgt")
    disable_logging(Logging.Debug)

    for dtii in dtvec

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        sp_ds = read(e5ds,evar_sp,ereg,dtii,quiet=true)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]
        wp_ds = read(e5ds,evar_wp,ereg,dtii,quiet=true)
        wp_sc = wp_ds[evar_wp.ID].attrib["scale_factor"]
        wp_of = wp_ds[evar_wp.ID].attrib["add_offset"]
        wp_mv = wp_ds[evar_wp.ID].attrib["missing_value"]
        wp_fv = wp_ds[evar_wp.ID].attrib["_FillValue"]

        @info "$(now()) - ConvectionIsotopes - Calculating the $(evar_wp.name) data over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        NCDatasets.load!(sp_ds[evar_sp.ID].var,sptmp,:,:,:)
        int2real!(sp,sptmp,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

        NCDatasets.load!(wp_ds[evar_wp.ID].var,wptmp,:,:,:)
        int2real!(wp,wptmp,scale=wp_sc,offset=wp_of,mvalue=wp_mv,fvalue=wp_fv)

        for it = 1 : 12, ilat = 1 : nlat, ilon = 1 : nlon
            wσ[ilon,ilat,it] = wp[ilon,ilat,it] / sp[ilon,ilat,it]
        end

        ERA5Reanalysis.save(wσ,dtii,e5ds,evar_wσ,ereg,lsd)

    end

end

function create_wσ(
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

    @info "$(now()) - ConvectionIsotopes - Preallocating arrays ..."

    wp = Array{Float32,3}(undef,nlon,nlat,31); wptmp = Array{Int16,3}(undef,nlon,nlat,31)
    sp = Array{Float32,3}(undef,nlon,nlat,31); sptmp = Array{Int16,3}(undef,nlon,nlat,31)
    wσ = Array{Float32,3}(undef,nlon,nlat,31)

    @info "$(now()) - ConvectionIsotopes - Extracting the ERA5Variable Information for Surface Pressure, Vertical Winds and Vertical Wind Weighted Column Pressure ..."

    disable_logging(Logging.Warn)
    evar_sp = SingleVariable("sp")
    evar_wp = SingleVariable("p_wwgt")
    evar_wσ = SingleVariable("σ_wwgt")
    disable_logging(Logging.Debug)

    for dtii in dtvec

        ndy  = daysinmonth(dtii)
        spii = view(sp,:,:,1:ndy); spti = view(sptmp,:,:,1:ndy)
        wpii = view(wp,:,:,1:ndy); wpti = view(wptmp,:,:,1:ndy)

        @info "$(now()) - ConvectionIsotopes - Extracting the Surface Pressure and Vertical Wind dataset over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        sp_ds = read(e5ds,evar_sp,ereg,dtii,smooth=smooth,smoothtime=smoothtime,quiet=true)
        sp_sc = sp_ds[evar_sp.ID].attrib["scale_factor"]
        sp_of = sp_ds[evar_sp.ID].attrib["add_offset"]
        sp_mv = sp_ds[evar_sp.ID].attrib["missing_value"]
        sp_fv = sp_ds[evar_sp.ID].attrib["_FillValue"]
        
        wp_ds = read(e5ds,evar_wp,ereg,dtii,smooth=smooth,smoothtime=smoothtime,quiet=true)
        wp_sc = wp_ds[evar_wp.ID].attrib["scale_factor"]
        wp_of = wp_ds[evar_wp.ID].attrib["add_offset"]
        wp_mv = wp_ds[evar_wp.ID].attrib["missing_value"]
        wp_fv = wp_ds[evar_wp.ID].attrib["_FillValue"]

        @info "$(now()) - ConvectionIsotopes - Calculating the $(evar_wp.name) data over the $(ereg.geo.name) Region for $(year(dtii)) ..."

        NCDatasets.load!(sp_ds[evar_sp.ID].var,spti,:,:,:)
        int2real!(spii,spti,scale=sp_sc,offset=sp_of,mvalue=sp_mv,fvalue=sp_fv)

        NCDatasets.load!(wp_ds[evar_wp.ID].var,wpti,:,:,:)
        int2real!(wpii,wpti,scale=wp_sc,offset=wp_of,mvalue=wp_mv,fvalue=wp_fv)

        for it = 1 : ndy, ilat = 1 : nlat, ilon = 1 : nlon
            wσ[ilon,ilat,it] = wpii[ilon,ilat,it] / spii[ilon,ilat,it]
        end

        ERA5Reanalysis.save(
            view(wσ,:,:,1:ndy),dtii,e5ds,evar_wσ,ereg,lsd,
            smooth=smooth,smoothtime=smoothtime
        )

    end

end