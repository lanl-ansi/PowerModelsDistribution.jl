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
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssDataObject
    dataobject = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssControlObject
    controlobject = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssGrowthshape
    growthshape = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssXycurve
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    xycurve = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :points ∈ raw_fields
        xarray = Float64[]
        yarray = Float64[]

        i = 1
        for point in xycurve.points
            if i % 2 == 1
                push!(xarray, point)
            else
                push!(yarray, point)
            end
            i += 1
        end
    else
        xarray = getproperty(xycurve, :xarray, Float64[])
        yarray = getproperty(xycurve, :yarray, Float64[])
    end

    npts = min(length(xarray), length(yarray))

    points = Float64[]
    for (x, y) in zip(xarray, yarray)
        push!(points, x)
        push!(points, y)
    end

    xycurve.npts = npts
    xycurve.points = points
    xycurve.xarray = xarray
    xycurve.yarray = yarray

    return xycurve
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssSpectrum
    spectrum = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssCapcontrol
    capcontrol = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssRegcontrol
    regcontrol = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssLinegeometry
    linegeometry = _apply_property_pairs(T(), property_pairs, dss, dss_raw)
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssWiredata
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    wiredata = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :diam ∈ raw_fields || :radius ∈ raw_fields
        radunits = wiredata.radunits
        if :diam ∈ raw_fields
            diam = wiredata.diam
            radius = wiredata.diam / 2
        elseif :radius ∈ raw_fields
            radius = wiredata.radius
            diam = 2 * wiredata.radius
        end

        if :gmrac ∉ raw_fields
            gmrac = 0.7788 * wiredata.radius
            gmrunits = wiredata.radunits
        else
            gmrac = wiredata.gmrac
            gmrunits = wiredata.gmrunits
        end
    elseif :gmrac ∈ raw_fields && :diag ∉ raw_fields && :radius ∉ raw_fields
        radius = wiredata.gmrac / 0.7788
        diam = 2 * wiredata.radius
        radunits = wiredata.gmrunits
        gmrac = wiredata.gmrac
        gmrunits = radunits
    else
        radius = 1
        diam = 2
        gmrac = 0.7788
        radunits = "none"
        gmrunits = "none"
    end

    rdc = :rdc ∈ raw_fields ? wiredata.rdc : :rac ∈ raw_fields ? wiredata.rac : 0.0
    rac = :rac ∈ raw_fields ? wiredata.rac : :rdc ∈ raw_fields ? wiredata.rdc : 0.0

    capradius = :capradius ∈ raw_fields ? wiredata.capradius : radius

    normamps = :normamps ∈ raw_fields ? wiredata.normamps : :emergamps ∈ raw_fields ? wiredata.emergamps / 1.5 : 400.0
    emergamps = :emergamps ∈ raw_fields ? wiredata.emergamps : :normamps ∈ raw_fields ? wiredata.normamps * 1.5 : 600.0

    wiredata.rdc = rdc / _convert_to_meters[wiredata.runits]
    wiredata.rac = rac / _convert_to_meters[wiredata.runits]
    wiredata.runits = "m"

    wiredata.gmrac = gmrac * _convert_to_meters[gmrunits]
    wiredata.gmrunits = "m"

    wiredata.radius = radius * _convert_to_meters[radunits]
    wiredata.capradius = capradius * _convert_to_meters[radunits]
    wiredata.diam = diam * _convert_to_meters[radunits]
    wiredata.radunits = "m"

    wiredata.normamps = normamps
    wiredata.emergamps = emergamps

    return wiredata
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssLinespacing
    linespacing = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    linespacing.x = linespacing.x .* _convert_to_meters[linespacing.units]
    linespacing.h = linespacing.h .* _convert_to_meters[linespacing.units]
    linespacing.units = "m"

    linespacing.fx = linespacing.x
    linespacing.fh = linespacing.h
    linespacing.funits = fill("m", linespacing.nconds)

    return linespacing
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssCndata
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    cndata = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :diam ∈ raw_fields || :radius ∈ raw_fields
        radunits = cndata.radunits
        if :diam ∈ raw_fields
            diam = cndata.diam
            radius = cndata.diam / 2
        elseif :radius ∈ raw_fields
            radius = cndata.radius
            diam = 2 * cndata.radius
        end

        if :gmrac ∉ raw_fields
            gmrac = 0.7788 * cndata.radius
            gmrunits = cndata.radunits
        else
            gmrac = cndata.gmrac
            gmrunits = cndata.gmrunits
        end
    elseif :gmrac ∈ raw_fields && :diam ∉ raw_fields && :radius ∉ raw_fields
        radius = cndata.gmrac / 0.7788
        diam = 2 * cndata.radius
        radunits = cndata.gmrunits
        gmrac = cndata.gmrac
        gmrunits = radunits
    else
        radius = 1.0
        diam = 2.0
        gmrac = 0.7788
        radunits = "none"
        gmrunits = "none"
    end

    if :diastrand ∈ raw_fields && :gmrstrand ∉ raw_fields
        diastrand = cndata.diastrand
        gmrstrand = diastrand / 2 * 0.7788
    elseif :gmrstrand ∈ raw_fields && :diastrand ∉ raw_fields
        gmrstrand = cndata.gmrstrand
        diastrand = gmrstrand / 0.7788 * 2
    elseif :diastrand ∈ raw_fields && :gmrstrand ∈ raw_fields
        diastrand = cndata.diastrand
        gmrstrand = cndata.gmrstrand
    else
        diastrand = 2.0
        gmrstrand = 0.7788
    end

    rac = :rac ∈ raw_fields ? cndata.rac : :rdc ∈ raw_fields ? cndata.rdc * 1.02 : 1.02
    rdc = :rdc ∈ raw_fields ? cndata.rdc : :rac ∈ raw_fields ? cndata.rac / 1.02 : 1.0

    normamps = :normamps ∈ raw_fields ? cndata.normamps : :emergamps ∈ raw_fields ? cndata.emergamps / 1.5 : 400.0
    emergamps = :emergamps ∈ raw_fields ? cndata.emergamps : :normamps ∈ raw_fields ? cndata.normamps * 1.5 : 600.0

    cndata.diacable = cndata.diacable * _convert_to_meters[radunits]
    cndata.diains = cndata.diains * _convert_to_meters[radunits]
    cndata.diam = diam * _convert_to_meters[radunits]
    cndata.diastrand = diastrand * _convert_to_meters[radunits]
    cndata.inslayer = cndata.inslayer * _convert_to_meters[radunits]
    cndata.radius = radius * _convert_to_meters[radunits]
    cndata.radunits = "m"

    cndata.gmrac = gmrac * _convert_to_meters[gmrunits]
    cndata.gmrstrand = gmrstrand * _convert_to_meters[gmrunits]
    cndata.gmrunits = "m"

    cndata.rac = rac / _convert_to_meters[cndata.runits]
    cndata.rdc = rdc / _convert_to_meters[cndata.runits]
    cndata.rstrand = cndata.rstrand / _convert_to_meters[cndata.runits]
    cndata.runits = "m"

    cndata.normamps = normamps
    cndata.emergamps = emergamps

    return cndata
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel) where T <: DssTsdata
    raw_fields = collect(Symbol(x.first) for x in property_pairs)

    tsdata = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :diam ∈ raw_fields || :radius ∈ raw_fields
        radunits = tsdata.radunits
        if :diam ∈ raw_fields
            diam = tsdata.diam
            radius = tsdata.diam / 2
        elseif :radius ∈ raw_fields
            radius = tsdata.radius
            diam = 2 * tsdata.radius
        end

        if :gmrac ∉ raw_fields
            gmrac = 0.7788 * tsdata.radius
            gmrunits = tsdata.radunits
        else
            gmrac = tsdata.gmrac
            gmrunits = tsdata.gmrunits
        end
    elseif :gmrac ∈ raw_fields && :diag ∉ raw_fields && :radius ∉ raw_fields
        radius = tsdata.gmrac / 0.7788
        diam = 2 * tsdata.radius
        radunits = tsdata.gmrunits
        gmrac = tsdata.gmrac
        gmrunits = radunits
    else
        radius = 1.0
        diam = 2.0
        gmrac = 0.7788
        radunits = "none"
        gmrunits = "none"
    end

    rac = :rac ∈ raw_fields ? tsdata.rac : :rdc ∈ raw_fields ? tsdata.rdc * 1.02 : 1.02
    rdc = :rdc ∈ raw_fields ? tsdata.rdc : :rac ∈ raw_fields ? tsdata.rac / 1.02 : 1.0

    normamps = :normamps ∈ raw_fields ? tsdata.normamps : :emergamps ∈ raw_fields ? tsdata.emergamps / 1.5 : 400.0
    emergamps = :emergamps ∈ raw_fields ? tsdata.emergamps : :normamps ∈ raw_fields ? tsdata.normamps * 1.5 : 600.0

    tsdata.diacable = tsdata.diacable * _convert_to_meters[radunits]
    tsdata.diains = tsdata.diains * _convert_to_meters[radunits]
    tsdata.diam = diam * _convert_to_meters[radunits]
    tsdata.diashield = tsdata.diashield * _convert_to_meters[radunits]
    tsdata.inslayer = tsdata.inslayer * _convert_to_meters[radunits]
    tsdata.tapelayer = tsdata.tapelayer * _convert_to_meters[radunits]
    tsdata.radunits = "m"

    tsdata.gmrac = gmrac * _convert_to_meters[gmrunits]
    tsdata.gmrunits = "m"

    tsdata.rac = rac / _convert_to_meters[tsdata.runits]
    tsdata.rdc = rdc / _convert_to_meters[tsdata.runits]
    tsdata.runits = "m"

    tsdata.normamps = normamps
    tsdata.emergamps = emergamps

    return tsdata
end

