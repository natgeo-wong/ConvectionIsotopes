using DelimitedFiles
using DrWatson
using Glob

function int2real!(
    oarray :: AbstractArray{FT},
    iarray :: AbstractArray{Int16};
    scale  :: Real,
    offset :: Real,
    fvalue :: Int16,
    mvalue :: Int16
) where FT <: Real

    for ii = 1 : length(iarray)

        if (iarray[ii] == fvalue) || (iarray[ii] == mvalue)
              oarray[ii] = FT(NaN)
        else; oarray[ii] = iarray[ii] * scale + offset
        end

    end

    return

end

function inplaceadd!(
    oarray :: AbstractArray{<:Real},
    iarray :: AbstractArray{<:Real},
    scale  :: Real = 1
)

    if length(orray) == length(iarray)
        for ii in eachindex(iarray)
            oarray[ii] += iarray[ii] * scale
        end
    end

    return

end

function loadflightpaths()
    data = Vector{Array{Float64,2}}(undef,22)
    fvec = glob("OTRECrf*",datadir("flightpaths"))
    for ii in 1 : length(fvec)
        data[ii] = readdlm(fvec[ii],',')
    end
    return data
end