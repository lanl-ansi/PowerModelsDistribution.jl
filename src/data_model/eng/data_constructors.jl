"""
"""
function create_eng_object(::Type{T}, dss_obj::DssLoadshape; import_all::Bool=false)::NTuple{2,T} where T <: EngTimeSeries
    eng_obj_p = T(;
        name = "$(dss_obj.name)_p",
        time = dss_obj.hour,
        offset = 0.0,
        replace = dss_obj.useactual,
        values = dss_obj.pmult,
        source_id = "loadshape.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )

    eng_obj_q = T(;
        name = "$(dss_obj.name)_q",
        time = dss_obj.hour,
        offset = 0.0,
        replace = dss_obj.useactual,
        values = dss_obj.qmult,
        source_id = "loadshape.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )

    return (eng_obj_p, eng_obj_q)
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssLinecode; import_all::Bool=false)::T where T <: EngLinecode
    nphases = dss_obj.nphases

    T(;
        name = dss_obj.name,
        rs = dss_obj.rmatrix,
        xs = dss_obj.xmatrix,
        b_fr = dss_obj.cmatrix ./ 2.0,
        b_to = dss_obj.cmatrix ./ 2.0,
        g_fr = fill(0.0, nphases, nphases),
        g_to = fill(0.0, nphases, nphases),
        cm_ub = fill(dss_obj.emergamps, nphases),
        source_id = "linecode.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssXfmrcode; import_all::Bool=false)::T where T <: EngXfmrcode
    nphases = dss_obj.phases
    nrw = dss_obj.windings

    T(;
        name = dss_obj.name,
        tm_nom = missing,
        tm_set = Vector{Vector{Float64}}([fill(tap, nphases) for tap in dss_obj.taps]),
        tm_lb = Vector{Vector{Float64}}(fill(fill(dss_obj.mintap, nphases), nrw)),
        tm_ub = Vector{Vector{Float64}}(fill(fill(dss_obj.maxtap, nphases), nrw)),
        tm_fix = Vector{Vector{Bool}}(fill(ones(Bool, nphases), nrw)),
        tm_step = Vector{Vector{Float64}}(fill(fill(1/32, nphases), nrw)),
        vm_nom = Vector{Float64}(dss_obj.kvs),
        sm_nom = Vector{Float64}(dss_obj.kvas),
        configurations = Vector{ConnConfig}(dss_obj.conns),
        rw = Vector{Float64}(dss_obj["%rs"] ./ 100),
        noloadloss = dss_obj["%noloadloss"] / 100,
        cmag = dss_obj["%imag"] / 100,
        xsc = nrw == 2 ? [dss_obj.xhl / 100] : [dss_obj.xhl, dss_obj.xht, dss_obj.xlt] ./ 100,
        sm_ub = dss_obj.emerghkva,
        source_id = "xfmrcode.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssRegcontrol; import_all::Bool=false)::T where T <: EngTransformerControls
    T(;
        windings = Int[dss_obj.winding],
        terminals = Vector{Int}[Int[dss_obj.ptphase]],
        vreg = Vector{Float64}[Float64[dss_obj.vreg]],
        band = Vector{Float64}[Float64[dss_obj.band]],
        ptratio = Vector{Float64}[Float64[dss_obj.ptratio]],
        ctprim = Vector{Float64}[Float64[dss_obj.ctprim]],
        r = Vector{Float64}[Float64[dss_obj.r]],
        x = Vector{Float64}[Float64[dss_obj.x]],
        dss = import_all ? Vector{DssRegcontrol}[DssRegcontrol[dss_obj]] : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssCapcontrol; import_all::Bool=false)::T where T <: EngShuntControls
    eng_obj = T(;
        type = CapControlType[dss_obj.type],
        elements = String[dss_obj.element],
        terminals = Int[dss_obj.terminal],
        onsetting = Float64[dss_obj.onsetting],
        offsetting = Float64[dss_obj.offsetting],
        voltoverride = Float64[dss_obj.voltoverride],
        ptratio = Float64[dss_obj.ptratio],
        ctratio = Float64[dss_obj.ctratio],
        vm_lb = Float64[dss_obj.vmin],
        vm_ub = Float64[dss_obj.vmax],
        dss = import_all ? DssCapcontrol[dss_obj] : missing,
    )
end
