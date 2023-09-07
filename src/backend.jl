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