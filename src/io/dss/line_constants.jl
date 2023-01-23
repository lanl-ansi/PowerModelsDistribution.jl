"Permeability of free space (N/A^2)"
const μ₀ = 4π*10^-7

"Permittivity of free space (F/m)"
const ε₀ = 8.8541878176e-12

"resistivity of copper tape shield (Ω-m)"
const ρₜₛ = 2.3718e-8


"gets line geometry data for line, including applying line spacing if specified"
function _get_geometry_data(data_dss::Dict{String,<:Any}, geometry_id::String)::Dict{String,Any}
    geom_obj = data_dss["linegeometry"][geometry_id]
    _apply_like!(geom_obj, data_dss, "linegeometry")

    geometry = _apply_ordered_properties(_create_linegeometry(geometry_id; _to_kwargs(geom_obj)...), geom_obj)

    if !isempty(get(geometry, "spacing", "")) && !isempty(get(get(data_dss,"linespacing",Dict()), string(geometry["spacing"]), Dict()))
        spacing = _get_spacing_data(data_dss, string(geometry["spacing"]))

        @assert geometry["nconds"] == spacing["nconds"] "Nconds on linegeometry.$(geometry_id) doesn't match nconds on linespacing.$(geometry["spacing"])"

        geometry["fx"] = spacing["fx"]
        geometry["fh"] = spacing["fh"]
    end

    return geometry
end


"get line spacing data for line or line geometry"
function _get_spacing_data(data_dss::Dict{String,<:Any}, spacing_id::String)
    spacing_obj = data_dss["linespacing"][spacing_id]
    _apply_like!(spacing_obj, data_dss, "linespacing")

    return _apply_ordered_properties(_create_linespacing(spacing_id; _to_kwargs(spacing_obj)...), spacing_obj)
end


"gets overhead wire data for line geometry"
function _get_wire_data(data_dss::Dict{String,<:Any}, wires::Vector{String})::Dict{String,Any}
    wiredata = Dict{String,Any}(id => get(get(data_dss,"wiredata",Dict()), id, Dict{String,Any}()) for id in filter(x->!isempty(x), wires))
    @assert !isempty(wiredata) && all(!isempty(wd) for (_,wd) in wiredata) "Some wiredata is missing, cannot continue"

    return Dict{String,Any}(
        id => _apply_ordered_properties(_create_wiredata(id; _to_kwargs(wd)...), wd) for (id,wd) in wiredata
    )
end


"gets concentric neutral cable data for line geometry"
function _get_cncable_data(data_dss::Dict{String,<:Any}, cncables::Vector{String})::Dict{String,Any}
    cncabledata = Dict{String,Any}(id => get(get(data_dss,"cndata",Dict()), id, Dict{String,Any}()) for id in filter(x->!isempty(x), cncables))
    @assert !isempty(cncabledata) && all(!isempty(cncd) for (_,cncd) in cncabledata) "Some cndata is missing, cannot continue"

    return Dict{String,Any}(
        id => _apply_ordered_properties(_create_cndata(id; _to_kwargs(cncd)...), cncd) for (id,cncd) in cncabledata
    )
end


"gets tape shielded cable data for line geometry"
function _get_tscable_data(data_dss::Dict{String,<:Any}, tscables::Vector{String})::Dict{String,Any}
    tscabledata = Dict{String,Any}(id => get(get(data_dss,"tsdata",Dict()), id, Dict{String,Any}()) for id in filter(x->!isempty(x), tscables))
    @assert !isempty(tscabledata) && all(!isempty(tscd) for (_,tscd) in tscabledata) "Some tsdata is missing, cannot continue"

    return Dict{String,Any}(
        id => _apply_ordered_properties(_create_tsdata(id; _to_kwargs(tscd)...), tscd) for (id,tscd) in tscabledata
    )
end


"""
    calculate_line_constants(data_dss::Dict{String,<:Any}, line_defaults::Dict{String,<:Any})::Tuple{Matrix{Complex},Matrix{Complex}}

Calculates line impedance and shunt admittance matrices for lines with line geometry, line spacing, wiredata, cncable, and/or tscable properties.
"""
function calculate_line_constants(data_dss::Dict{String,<:Any}, line_defaults::Dict{String,<:Any})::Tuple{Matrix{Complex},Matrix{Complex}}
    geometry = !isempty(line_defaults["geometry"]) ? _get_geometry_data(data_dss, string(line_defaults["geometry"])) : missing

    cncables = !ismissing(geometry) && !isempty(geometry["cncables"]) ? filter(x->!isempty(x), geometry["cncables"]) : !isempty(line_defaults["cncables"]) ? filter(x->!isempty(x), line_defaults["cncables"]) : missing
    ncncables = !ismissing(cncables) ? length(cncables) : missing
    cndata = !ismissing(cncables) ? _get_cncable_data(data_dss, cncables) : missing

    tscables = !ismissing(geometry) && !isempty(geometry["tscables"]) ? filter(x->!isempty(x), geometry["tscables"]) : !isempty(line_defaults["tscables"]) ? filter(x->!isempty(x), line_defaults["tscables"]) : missing
    ntscables = !ismissing(tscables) ? length(tscables) : missing
    tsdata = !ismissing(tscables) ? _get_tscable_data(data_dss, tscables) : missing

    spacing = !isempty(line_defaults["spacing"]) ? _get_spacing_data(data_dss, string(line_defaults["spacing"])) : missing

    wires = !ismissing(geometry) && !isempty(geometry["wires"]) ? filter(x->!isempty(x), geometry["wires"]) : !isempty(line_defaults["wires"]) ? filter(x->!isempty(x), line_defaults["wires"]) : missing
    nwires = !ismissing(wires) ? length(wires) : missing
    wiredata = !ismissing(wires) ? _get_wire_data(data_dss, wires) : missing

    @assert count(.!(ismissing.([wires, cncables, tscables]))) > 0 "No wire, tscable, or cncable data is defined on line.$(line_defaults["name"])"

    ω = 2π * get(line_defaults, "basefreq", get(get(data_dss, "options", Dict()), "defaultbasefreq", 60.0))
    ω₀ = 2π * get(get(data_dss, "options", Dict()), "defaultbasefreq", 60.0)

    x = !ismissing(geometry) ? geometry["fx"] : !ismissing(spacing) ? spacing["fx"] : missing
    y = !ismissing(geometry) ? geometry["fh"] : !ismissing(spacing) ? spacing["fh"] : missing

    gmr = !ismissing(wiredata) ? [wiredata[wire]["gmrac"] for wire in wires] : missing
    capradius = !ismissing(wiredata) ? [wiredata[wire]["capradius"] for wire in wires] : missing

    rac = !ismissing(wiredata) ? [wiredata[wire]["rac"] for wire in wires] : missing
    rdc = !ismissing(wiredata) ? [wiredata[wire]["rdc"] for wire in wires] : missing

    nphases = !ismissing(geometry) ? geometry["nphases"] : !ismissing(spacing) ? spacing["nphases"] : line_defaults["phases"]
    nconds = !ismissing(geometry) ? geometry["nconds"] : !ismissing(spacing) ? spacing["nconds"] : nphases

    rac = !ismissing(cndata) ? [cndata[cable]["rac"] for cable in cncables] : rac
    rdc = !ismissing(cndata) ? [cndata[cable]["rdc"] for cable in cncables] : rdc
    gmr = !ismissing(cndata) ? [cndata[cable]["gmrac"] for cable in cncables] : gmr

    rstrand = !ismissing(cndata) ? [cndata[cable]["rstrand"] for cable in cncables] : missing
    kstrand = !ismissing(cndata) ? [cndata[cable]["k"] for cable in cncables] : missing
    diacable = !ismissing(cndata) ? [cndata[cable]["diacable"] for cable in cncables] : missing
    diastrand = !ismissing(cndata) ? [cndata[cable]["diastrand"] for cable in cncables] : missing
    gmrstrand = !ismissing(cndata) ? [cndata[cable]["gmrstrand"] for cable in cncables] : missing
    epsr = !ismissing(cndata) ? [cndata[cable]["epsr"] for cable in cncables] : missing
    diains = !ismissing(cndata) ? [cndata[cable]["diains"] for cable in cncables] : missing
    inslayer = !ismissing(cndata) ? [cndata[cable]["inslayer"] for cable in cncables] : missing

    rac = !ismissing(tsdata) ? [tsdata[cable]["rac"] for cable in tscables] : rac
    rdc = !ismissing(tsdata) ? [tsdata[cable]["rdc"] for cable in tscables] : rdc
    gmr = !ismissing(tsdata) ? [tsdata[cable]["gmrac"] for cable in tscables] : gmr

    if count(.!(ismissing.([wires, cncables, tscables]))) > 1
        @assert sum(filter(x->!ismissing(x), [nwires, ncncables, ntscables])) == nconds "not enough wire/cable data for specified conductors on line.$(line_defaults["name"])"
        push!(rac, [wiredata[wire]["rac"] for wire in wires]...)
        push!(rdc, [wiredata[wire]["rdc"] for wire in wires]...)
        push!(gmr, [wiredata[wire]["gmrac"] for wire in wires]...)
    end

    diashield = !ismissing(tsdata) ? [tsdata[cable]["diashield"] for cable in tscables] : missing
    tapelayer = !ismissing(tsdata) ? [tsdata[cable]["tapelayer"] for cable in tscables] : missing
    tapelap = !ismissing(tsdata) ? [tsdata[cable]["tapelap"] for cable in tscables] : missing

    epsr = !ismissing(tsdata) ? [tsdata[cable]["epsr"] for cable in tscables] : epsr
    diains = !ismissing(tsdata) ? [tsdata[cable]["diains"] for cable in tscables] : diains
    inslayer = !ismissing(tsdata) ? [tsdata[cable]["inslayer"] for cable in tscables] : inslayer

    reduce = !ismissing(geometry) ? geometry["reduce"] : !ismissing(spacing) ? (line_defaults["phases"] == spacing["nconds"] ? false : true) : missing

    earth_model = !isempty(line_defaults["earthmodel"]) ? line_defaults["earthmodel"] : get(get(data_dss, "options", Dict()), "earthmodel", "deri")

    rho = line_defaults["rho"]

    Z, Y =  calculate_line_constants(
        x,
        y,
        ω,
        gmr,
        capradius,
        nconds,
        earth_model,
        rac,
        ω₀,
        rdc,
        rho,
        nphases,
        rstrand,
        kstrand,
        diacable,
        diastrand,
        gmrstrand,
        epsr,
        diains,
        inslayer,
        diashield,
        tapelayer,
        tapelap
    )

    if reduce
        Z, Y = _kron(Z, Y, nphases)
    end

    return Z, Y
end


"""
    calculate_line_constants(
        x::Vector{<:Real},
        y::Vector{<:Real},
        ω::Real,
        gmr::Vector{<:Real},
        ::Union{Missing,Vector{<:Real}},
        nconds::Int,
        earth_model::String,
        R_ac::Vector{<:Real},
        ω₀::Real,
        R_dc::Vector{<:Real},
        ρₑ::Real,
        nphases::Int,
        R_strand::Vector{<:Real},
        n_strand::Vector{<:Real},
        d_cable::Vector{<:Real},
        d_strand::Vector{<:Real},
        gmr_strand::Vector{<:Real},
        ε_ins::Vector{<:Real},
        d_ins::Vector{<:Real},
        t_ins::Vector{<:Real},
        ::Missing,
        ::Missing,
        ::Missing
    )::Tuple{Matrix{Complex},Matrix{Complex}}

Calculates the impedance and shunt admittance of concentric neutral cables

# References

- Nasser Tleis, Power Systems Modelling and Fault Analysis (Second Edition), Academic Press, 2019, ISBN 9780128151174.
- William H Kersting, Distribution System Modeling and Analysis (Forth Edition), CRC Press, 2018, ISBN 9781498772136.
- Andrea Ballanti, [Cable Modeling in OpenDSS](https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/Doc/TechNote%20CableModelling.pdf?format=raw), 2017.
"""
function calculate_line_constants(
    x::Vector{<:Real},
    y::Vector{<:Real},
    ω::Real,
    gmr::Vector{<:Real},
    ::Union{Missing,Vector{<:Real}},
    nconds::Int,
    earth_model::String,
    R_ac::Vector{<:Real},
    ω₀::Real,
    R_dc::Vector{<:Real},
    ρₑ::Real,
    nphases::Int,
    R_strand::Vector{<:Real},
    n_strand::Vector{<:Real},
    d_cable::Vector{<:Real},
    d_strand::Vector{<:Real},
    gmr_strand::Vector{<:Real},
    ε_ins::Vector{<:Real},
    d_ins::Vector{<:Real},
    t_ins::Vector{<:Real},
    ::Missing,
    ::Missing,
    ::Missing)::Tuple{Matrix{Complex},Matrix{Complex}}

    N = nconds + nphases

    Z = zeros(ComplexF64, N, N)
    Y = zeros(ComplexF64, nphases, nphases)

    for i in 1:nconds
        ii = i+nconds
        for j in 1:nconds
            jj = j+nconds

            Z_ije = Z_ije_ss = calc_earth_return_path_impedance(i, j, x, y, ρₑ, earth_model, ω, ω₀)
            D_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)

            if j <= i
                if i == j
                    Z_ic = calc_internal_impedance(R_ac[i], R_dc[i], earth_model, ω).re
                    Z_ig = 1im * ω*μ₀/2π*log(1/gmr[i])

                    Z[i,i] = Z_ic + Z_ig + Z_ije

                    if i <= nphases
                        r_cn = 1/2*(d_cable[i]-d_strand[i])
                        gmr_cn = (gmr_strand[i]*n_strand[i]*(r_cn^(n_strand[i])))^(1 / n_strand[i])

                        Z_ic_ss = R_strand[i] / n_strand[i]
                        Z_ig_ss = 1im * ω*μ₀/2π*log(1/gmr_cn)

                        Z[ii,ii] = Z_ic_ss + Z_ig_ss + Z_ije_ss
                    end
                else
                    Z_ijg = 1im * ω*μ₀/2π*log(1/D_ij)

                    Z[i,j] = Z[j,i] = Z_ijg + Z_ije

                    if i <= nphases
                        jj = j+nconds
                        D_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)

                        Z_ijg = 1im * ω*μ₀/2π*log(1 / D_ij)

                        Z[ii,jj] = Z[jj,ii] = Z_ijg + Z_ije
                    end
                end
            end

            if i <= nphases
                r_cn = 0.5 * (d_cable[i]-d_strand[i])
                if i == j
                    D_ij = r_cn
                else
                    D_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
                    D_ij = (D_ij^(n_strand[i]) - r_cn^(n_strand[i]))^(1/n_strand[i])
                end

                Z_ijg = 1im * ω*μ₀/2π*log(1/D_ij)

                Z[ii,j] = Z[j,ii] = Z_ijg + Z_ije
            end
        end

        if i <= nphases
            r_outer = d_ins[i]/2
            r_inner = r_outer - t_ins[i]
            Y[i,i] = 1im * 2π * ε₀ * ε_ins[i] * ω / log(r_outer / r_inner)
        end
    end

    Z = _kron(Z, nconds)

    C = Y .* 1e9 ./ ω

    return Z, C
end


"""
    calculate_line_constants(
        x::Vector{<:Real},
        y::Vector{<:Real},
        ω::Real,
        gmr::Vector{<:Real},
        ::Union{Missing,Vector{<:Real}},
        nconds::Int,
        earth_model::String,
        R_ac::Vector{<:Real},
        ω₀::Real,
        R_dc::Vector{<:Real},
        ρₑ::Real,
        nphases::Int,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ε_ins::Vector{<:Real},
        d_ins::Vector{<:Real},
        t_ins::Vector{<:Real},
        d_shield::Vector{<:Real},
        t_tape::Vector{<:Real},
        lap_tape::Vector{<:Real}
    )::Tuple{Matrix{Complex},Matrix{Complex}}

Calculates the impedance and shunt admittance of tape shielded cables

# References

- Nasser Tleis, Power Systems Modelling and Fault Analysis (Second Edition), Academic Press, 2019, ISBN 9780128151174.
- William H Kersting, Distribution System Modeling and Analysis (Forth Edition), CRC Press, 2018, ISBN 9781498772136.
- Andrea Ballanti, [Cable Modeling in OpenDSS](https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/Doc/TechNote%20CableModelling.pdf?format=raw), 2017.
"""
function calculate_line_constants(
    x::Vector{<:Real},
    y::Vector{<:Real},
    ω::Real,
    gmr::Vector{<:Real},
    ::Union{Missing,Vector{<:Real}},
    nconds::Int,
    earth_model::String,
    R_ac::Vector{<:Real},
    ω₀::Real,
    R_dc::Vector{<:Real},
    ρₑ::Real,
    nphases::Int,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ε_ins::Vector{<:Real},
    d_ins::Vector{<:Real},
    t_ins::Vector{<:Real},
    d_shield::Vector{<:Real},
    t_tape::Vector{<:Real},
    lap_tape::Vector{<:Real})::Tuple{Matrix{Complex},Matrix{Complex}}

    N = nconds + nphases
    Z = zeros(Complex, N, N)
    Y = zeros(Complex, nphases, nphases)

    for i in 1:nconds
        ii = i+nconds
        for j in 1:nconds
            jj = j+nconds

            Z_ije = Z_ije_ss = calc_earth_return_path_impedance(i, j, x, y, ρₑ, earth_model, ω, ω₀)
            D_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)

            if j <= i
                if i == j
                    Z_ic = calc_internal_impedance(R_ac[i], R_dc[i], earth_model, ω).re
                    Z_ig = 1im * ω*μ₀/2π*log(1/gmr[i])

                    Z[i,i] = Z_ic + Z_ig + Z_ije

                    if i <= nphases
                        gmr_ts = 1/2 * (d_shield[i] - t_tape[i])

                        Z_ic_ss = 0.3183 * ρₜₛ / (d_shield[i] * t_tape[i] * sqrt(50 / (100-lap_tape[i])))
                        Z_ig_ss = 1im * ω*μ₀/2π*log(1/gmr_ts)

                        Z[ii,ii] = Z_ic_ss + Z_ig_ss + Z_ije_ss
                    end
                else
                    Z_ijg = 1im * ω*μ₀/2π*log(1/D_ij)

                    Z[i,j] = Z[j,i] = Z_ijg + Z_ije

                    if i <= nphases
                        jj = j+nconds
                        D_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)

                        Z_ijg = 1im * ω*μ₀/2π*log(1 / D_ij)

                        Z[ii,jj] = Z[jj,ii] = Z_ijg + Z_ije
                    end
                end
            end

            if i <= nphases
                if i == j
                    D_ij = 1/2 * (d_shield[i] - t_tape[i])
                end

                Z_ijg = 1im * ω*μ₀/2π*log(1/D_ij)

                Z[ii,j] = Z[j,ii] = Z_ijg + Z_ije
            end
        end

        if i <= nphases
            r_outer = d_ins[i]/2
            r_inner = r_outer - t_ins[i]
            Y[i,i] = (0+1im) * 2π * ε₀ * ε_ins[i] * ω / log(r_outer / r_inner)
        end
    end

    Z = _kron(Z, nconds)

    C = Y .* 1e9 ./ ω

    return Z, C
end


"""
    calculate_line_constants(
        x::Vector{<:Real},
        y::Vector{<:Real},
        ω::Real,
        gmr::Vector{<:Real},
        r::Vector{<:Real},
        nconds::Int,
        earth_model::String,
        R_ac::Vector{<:Real},
        ω₀::Real,
        R_dc::Vector{<:Real},
        ρₑ::Real,
        ::Int,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing,
        ::Missing
    )::Tuple{Matrix{Complex},Matrix{Complex}}

Calculates impedance and shunt admittance for overhead lines.

# References

- Nasser Tleis, Power Systems Modelling and Fault Analysis (Second Edition), Academic Press, 2019, ISBN 9780128151174.
- William H Kersting, Distribution System Modeling and Analysis (Forth Edition), CRC Press, 2018, ISBN 9781498772136.
"""
function calculate_line_constants(
    x::Vector{<:Real},
    y::Vector{<:Real},
    ω::Real,
    gmr::Vector{<:Real},
    r::Vector{<:Real},
    nconds::Int,
    earth_model::String,
    R_ac::Vector{<:Real},
    ω₀::Real,
    R_dc::Vector{<:Real},
    ρₑ::Real,
    ::Int,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing,
    ::Missing)::Tuple{Matrix{Complex},Matrix{Complex}}

    Z = zeros(ComplexF64, nconds, nconds)
    P = zeros(ComplexF64, nconds, nconds)

    for i in 1:nconds
        for j in 1:i
            if i == j
                Z_ic = calc_internal_impedance(R_ac[i], R_dc[i], earth_model, ω).re
                # Tleis Eq 3.11
                Z_ig = 1im * ω*μ₀/2π*log(1/gmr[i])
                Z_ie = calc_earth_return_path_impedance(i,j, x, y, ρₑ, earth_model, ω, ω₀)

                # Tleis Eq 3.5b
                Z[i,i] = Z_ic + Z_ig + Z_ie

                # Kersting Eq 5.9
                P[i,i] = 1im * -1/(2π * ε₀ * ω) * log(2 * y[i] / r[i])
            else
                d_ij = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)

                # Tleis Eq 3.12b
                X_ijg = (μ₀*ω/2π)*log(1/d_ij)
                Z[i,j] = Z[j,i] = 1im * X_ijg + calc_earth_return_path_impedance(i,j, x, y, ρₑ, earth_model, ω, ω₀)

                # Kersting Pg 82
                S_ij = sqrt((x[i] - x[j])^2 + (y[i] + y[j])^2)

                # Kersting Eq 5.10
                P[i,j] = P[j,i] = 1im * -1/(2π * ε₀ * ω) * log(S_ij/d_ij)
            end
        end
    end

    C = LinearAlgebra.pinv(P) .* 1e9 ./ ω

    return Z, C
end


"""
    calc_earth_return_path_impedance(i::Int, j::Int, x::Vector{<:Real}, y::Vector{<:Real}, ρₑ::Real, earth_model::String, ω::Real, ω₀::Real)

Calculates the earth return path impedance

# References

- Nasser Tleis, Power Systems Modelling and Fault Analysis (Second Edition), Academic Press, 2019, ISBN 9780128151174.
- A. Deri, G. Tevan, A. Semlyen, A. Castanheira, The Complex Ground Return Plane – A Simplified Model for Homogeneous and Multi-Layer Earth Return, IEEE Transactions on Power Apparatus and Systems, Vol. PAS-100, No. 8, pp. 3686-3693, August 1981.
"""
function calc_earth_return_path_impedance(i::Int, j::Int, x::Vector{<:Real}, y::Vector{<:Real}, ρₑ::Real, earth_model::String, ω::Real, ω₀::Real)
    # Tleis Eq 3.7c
    δ = 503.292*sqrt(ρₑ/(ω₀/2π))

    # Tleis Eq 3.15
    D_erc = 1.309125*δ

    if earth_model ∈ ["simplecarson", "carson"]
        Z_ij = (μ₀*ω/2π) * (π/4 + 1im * log(D_erc))
    elseif earth_model == "fullcarson"
        # Tleis Pg 111
        b1 = 1/(3*sqrt(2))
        b2 = 1/16
        b3 = b1/(3*5)
        b4 = b2/(4*6)
        b5 = -b3/(5*7)

        c2 = 1.3659315
        c4 = c2 + 1/4 + 1/6

        d2 = π/4*b2
        d4 = π/4*b4


        if i == j
            # Tleis Eq 3.14a
            D_ij = 2 * y[i]
            # Tleis Eq 3.14b
            θ_ij = 0
        else
            # Tleis Eq 3.14a
            D_ij = sqrt((x[i]-x[j])^2 + (y[i]+y[j])^2)
            # Tleis Eq 3.14b
            θ_ij = acos((y[i]+y[j])/D_ij)
        end

        # Tleis Eq 3.14a
        m_ij = sqrt(2) * D_ij / δ

        # Tleis Eq 3.13a
        R_ije = (π/8 - b1*m_ij*cos(θ_ij) + b2*m_ij^2*(log(exp(c2)/m_ij)*cos(2*θ_ij) + θ_ij*sin(2*θ_ij)) + b3*m_ij^3*cos(3*θ_ij) - d4*m_ij^4*cos(4*θ_ij) - b5*m_ij^5*cos(5*θ_ij))

        # Tleis Eq 3.13b
        X_ije = (1/2*log(1.85138/m_ij) + b1*m_ij*cos(θ_ij) - d2*m_ij^2*cos(2*θ_ij) + b3*m_ij^3*cos(3*θ_ij) - b4*m_ij^4*(log(exp(c4)/m_ij)*cos(4*θ_ij)+θ_ij*sin(4*θ_ij)) + b5*m_ij^5*cos(5*θ_ij))

        # Tleis Eq 3.17b
        X_ije += 1/2*log(D_ij)

        Z_ij = ω/2π*μ₀ * 2 * (R_ije + 1im * X_ije)
    elseif earth_model == "deri"
        # Deri Eq 18
        p = 1 / sqrt(1im * ω * μ₀ / ρₑ)

        if i == j
            # Deri Eq 3
            Z_ij = 1im * ω*μ₀/2π * log(2*(y[i] + p))
        else
            # Deri Eq 4
            Z_ij = 1im * ω*μ₀/2π * log(sqrt((y[i] + y[j] + 2*p)^2 + (x[i]-x[j])^2))
        end
    else
        error("earth model $(earth_model) not recognized")
    end

    return Z_ij
end


"""
    calc_internal_impedance(R_ac::Real, R_dc::Real, earth_model::String, ω::Real)

Calculates the internal impedance of the conductor

# References

- Nasser Tleis, Power Systems Modelling and Fault Analysis (Second Edition), Academic Press, 2019, ISBN 9780128151174.
"""
function calc_internal_impedance(R_ac::Real, R_dc::Real, earth_model::String, ω::Real)
    if occursin("carson", earth_model)
        L_i = μ₀/8π  # internal inductance
        Z_int = R_ac + 1im * ω*L_i
    elseif earth_model == "deri"
        δ = sqrt(R_dc*(ω/2π)*μ₀)
        m = (1+1im)*sqrt(ω/2π*μ₀/R_dc)
        Z_int = (1+1im) * (abs(m) > 35 ? (1+1im) : SpecialFunctions.besselj0(m) / SpecialFunctions.besselj1(m)) * δ / 2
    else
        error("earth model $(earth_model) not recognized")
    end

    return Z_int
end


"""
    _kron(Z::Matrix{T}, Y::Matrix{T}, nconds::Int)::Tuple{Matrix{T}, Matrix{T}} where T <: Complex

Kron reduces impedance and shunt admittance matrices down to size `(nconds, nconds)`
"""
function _kron(Z::Matrix{T}, Y::Matrix{T}, nconds::Int)::Tuple{Matrix{T}, Matrix{T}} where T <: Complex
    Zout = deepcopy(Z)
    Yout = deepcopy(Y)

    nrow, ncol = size(Z)
    @assert nrow == ncol "cannot kron reduce non-square matrix"

    _Z = zeros(ComplexF64, nconds, nconds)
    _Y = zeros(ComplexF64, nconds, nconds)
    N = nrow
    while N > nconds
        _Z = zeros(ComplexF64, N-1, N-1)
        _Y = zeros(ComplexF64, N-1, N-1)

        for i in 1:N-1
            for j in 1:N-1
                _Z[i,j] = Zout[i,j] - Zout[i,N]*Zout[j,N]/Zout[N,N]
                _Y[i,j] = Yout[i,j]
            end
        end

        Zout = _Z
        Yout = _Y
        N -= 1
    end

    return Zout, Yout
end


"""
    _kron(Z::Matrix{T}, nconds::Int)::Matrix{T} where T <: Complex

Kron reduces impedance matrix down to size `(nconds, nconds)`
"""
function _kron(Z::Matrix{T}, nconds::Int)::Matrix{T} where T <: Complex
    Zout = deepcopy(Z)

    nrow, ncol = size(Z)
    @assert nrow == ncol "cannot kron reduce non-square matrix"

    _Z = zeros(T, nconds, nconds)
    N = nrow
    while N > nconds
        _Z = zeros(T, N-1, N-1)
        for i in 1:N-1
            for j in 1:N-1
                _Z[i,j] = Zout[i,j] - Zout[i,N]*Zout[j,N]/Zout[N,N]
            end
        end
        Zout = _Z
        N -= 1
    end

    return Zout
end
