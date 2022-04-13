"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssLoad
    raw_fields = collect(x.first for x in property_pairs)

    load = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :kw ∈ raw_fields && :pf ∈ raw_fields
        load.kvar = sign(load.pf) * load.kw * sqrt(1.0 / load.pf^2 - 1.0)
        load.kva = abs(load.kw) + load.kvar^2
    elseif :kw ∈ raw_fields && :kvar ∈ raw_fields
        load.kva = abs(load.kw) + load.kvar^2
        if load.kva > 0.0
            load.pf = load.kw / load.kva
            if load.kvar != 0.0
                load.pf *= sign(load.kw * load.kvar)
            end
        end
    elseif :kva ∈ raw_fields && :pf ∈ raw_fields
        load.kw = load.kva * abs(load.pf)
        load.kvar = sign(load.pf) * load.kw * sqrt(1.0 / load.pf^2 - 1.0)
    elseif :pf ∈ raw_fields && load.pf != 0.88
        load.kvar = sign(load.pf) * load.kw * sqrt(1.0 / load.pf^2 - 1.0)
        load.kva = abs(load.kw) + load.kvar^2
    end

    return load
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssStorage
    raw_fields = collect(x.first for x in property_pairs)

    strg = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if Symbol("%stored") ∈ raw_fields
        strg.kwhstored = strg[Symbol("%stored")] / 100 * strg.kwhrated
    elseif :kwhstored ∈ raw_fields
        strg[Symbol("%stored")] = strg.kwhstored / strg.kwhrated * 100
    end

    if :kwrated ∈ raw_fields && :kva ∉ raw_fields
        strg.kva = strg.kwrated
    elseif :kva ∈ raw_fields && :kwrated ∉ raw_fields
        strg.kwrated = strg.kva
    end

    if :pf ∈ raw_fields && :kw ∈ raw_fields
        strg.kvar = sign(strg.pf) * strg.kw * sqrt(1 / strg.pf^2 - 1)
    elseif :kw ∈ raw_fields && :kvar ∈ raw_fields
        strg.pf = strg.kw == 0.0 ? 1.0 : strg.kvar / strg.kw
    end

    return strg
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssGenerator
    raw_fields = collect(x.first for x in property_pairs)

    gen = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    gen.kva = :kva ∉ raw_fields ? gen.kw * 1.2 : gen.kva

    gen.pf = :pf ∉ raw_fields ? 0.80 : gen.pf

    if :kw ∈ raw_fields && :pf ∈ raw_fields
        gen.kvar = sign(gen.pf) * gen.kw * sqrt(1.0 / gen.pf^2 - 1.0)
        gen.kva = abs(gen.kw) + gen.kvar^2
    elseif :kw ∈ raw_fields && :kvar ∈ raw_fields
        gen.kva = abs(gen.kw) + gen.kvar^2
        if gen.kva > 0.0
            gen.pf = gen.kw / gen.kva
            if gen.kvar != 0.0
                gen.pf *= sign(gen.kw * gen.kvar)
            end
        end
    elseif :kva ∈ raw_fields && :pf ∈ raw_fields
        gen.kw = gen.kva * abs(gen.pf)
        gen.kvar = sign(gen.pf) * gen.kw * sqrt(1.0 / gen.pf^2 - 1.0)
    elseif :pf ∈ raw_fields && gen.pf != 0.80
        gen.kvar = sign(gen.pf) * gen.kw * sqrt(1.0 / gen.pf^2 - 1.0)
        gen.kva = abs(gen.kw) + gen.kvar^2
    end

    gen.maxkvar = :maxkvar ∉ raw_fields ? gen.kvar * 2.0 : gen.maxkvar
    gen.minkvar = :minkvar ∉ raw_fields ? -gen.maxkvar : gen.minkvar

    return gen
end


"""
"""
function create_dss_object(::Type{T}, property_pairs::Vector{Pair{String,String}}, dss::OpenDssDataModel, dss_raw::OpenDssRawDataModel)::T where T <: DssPvsystem
    raw_fields = collect(x.first for x in property_pairs)

    pv = _apply_property_pairs(T(), property_pairs, dss, dss_raw)

    if :kw ∈ raw_fields && :pf ∈ raw_fields
        pv.kvar = sign(pv.pf) * pv.kw * sqrt(1.0 / pv.pf^2 - 1.0)
        pv.kva = abs(pv.kw) + pv.kvar^2
    elseif :kw ∈ raw_fields && :kvar ∈ raw_fields
        pv.kva = abs(pv.kw) + pv.kvar^2
        if pv.kva > 0.0
            pv.pf = pv.kw / pv.kva
            if pv.kvar != 0.0
                pv.pf *= sign(pv.kw * pv.kvar)
            end
        end
    elseif :kva ∈ raw_fields && :pf ∈ raw_fields
        pv.kw = pv.kva * abs(pv.pf)
        pv.kvar = sign(pv.pf) * pv.kw * sqrt(1.0 / pv.pf^2 - 1.0)
    elseif :pf ∈ raw_fields && pv.pf != 0.88
        pv.kvar = sign(pv.pf) * pv.kw * sqrt(1.0 / pv.pf^2 - 1.0)
        pv.kva = abs(pv.kw) + pv.kvar^2
    end

    return pv
end
