"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}})::T where T <: DssOptions
    options = _apply_property_pairs(T(), property_pairs)

    return options
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssLinecode
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    linecode = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :b1 ∈ raw_fields
        linecode.c1 = linecode.b1 / (2 * pi * linecode.basefreq)
    else
        linecode.b1 = linecode.c1 * (2 * pi * linecode.basefreq)
    end

    if linecode.nphases == 1
        linecode.r0 = linecode.r1
        linecode.x0 = linecode.x1
        linecode.c0 = linecode.c1
        linecode.b0 = linecode.b1
    else
        if :b0 ∈ raw_fields
            linecode.c0 = linecode.b0 / (2 * pi * linecode.basefreq)
        else
            linecode.b0 = linecode.c0 * (2 * pi * linecode.basefreq)
        end
    end

    Zs = (complex(linecode.r1, linecode.x1) * 2.0 + complex(linecode.r0, linecode.x0)) / 3.0
    Zm = (complex(linecode.r0, linecode.x0) - complex(linecode.r1, linecode.x1)) / 3.0

    Ys = (complex(0.0, 2 * pi * linecode.basefreq * linecode.c1) * 2.0 + complex(0.0, 2 * pi * linecode.basefreq * linecode.c0)) / 3.0
    Ym = (complex(0.0, 2 * pi * linecode.basefreq * linecode.c0) - complex(0.0, 2 * pi * linecode.basefreq * linecode.c1)) / 3.0

    Z = zeros(Complex{Float64}, linecode.nphases, linecode.nphases)
    Yc = zeros(Complex{Float64}, linecode.nphases, linecode.nphases)
    for i in 1:linecode.nphases
        Z[i,i] = Zs
        Yc[i,i] = Ys
        for j in 1:i-1
            Z[i,j] = Z[j,i] = Zm
            Yc[i,j] = Yc[j,i] = Ym
        end
    end

    linecode.rmatrix = :rmatrix ∈ raw_fields ? linecode.rmatrix : real(Z)
    linecode.xmatrix = :xmatrix ∈ raw_fields ? linecode.xmatrix : imag(Z)
    linecode.cmatrix = :cmatrix ∈ raw_fields ? linecode.cmatrix : imag(Yc) / (2 * pi * linecode.basefreq)

    linecode.r1 = linecode.r1 / _convert_to_meters[linecode.units]
    linecode.x1 = linecode.x1 / _convert_to_meters[linecode.units]
    linecode.r0 = linecode.r0 / _convert_to_meters[linecode.units]
    linecode.x0 = linecode.x0 / _convert_to_meters[linecode.units]
    linecode.c1 = linecode.c1 / _convert_to_meters[linecode.units]
    linecode.c0 = linecode.c0 / _convert_to_meters[linecode.units]
    linecode.rmatrix = linecode.rmatrix / _convert_to_meters[linecode.units]
    linecode.xmatrix = linecode.xmatrix / _convert_to_meters[linecode.units]
    linecode.cmatrix = linecode.cmatrix / _convert_to_meters[linecode.units]
    linecode.b1 = linecode.b1 / _convert_to_meters[linecode.units]
    linecode.b0 = linecode.b0 / _convert_to_meters[linecode.units]
    linecode.units = "m"

    return linecode
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssLoadshape
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    loadshape = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :minterval ∈ raw_fields
        interval = loadshape.minterval / 60
    elseif :sinterval ∈ raw_fields
        interval = loadshape.sinterval / 60 / 60
    else
        interval = loadshape.interval
    end

    # pmult = loadshape.pmult
    loadshape.qmult = :qmult ∉ raw_fields ? loadshape.pmult : loadshape.qmult

    # TODO
    npts = :npts ∈ raw_fields ? loadshape.npts : length(loadshape.pmult)
    # npts = get(kwargs, :npts, length(pmult) == 0 && length(qmult) == 0 ? 0 : minimum(Int[length(a) for a in [pmult, qmult] if length(a) > 0]))

    # loadshape.pmult = loadshape.pmult[1:npts]
    # loadshape.qmult = loadshape.qmult[1:npts]

    # loadshape.hour = loadshape.hour[1:npts]

    return loadshape
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: OpenDssDataObject
    dataobject = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: OpenDssControlObject
    controlobject = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end
