"""
    _solve_mn_mc_opb(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)

solve test mn mc problem
"""
function _solve_mn_mc_opb(data::Union{DistributionModel,String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, _build_mc_mn_opb; ref_extensions=[ref_add_connected_components!], multinetwork=true, kwargs...)
end


"""
    _build_mc_mn_opb(pm::AbstractUnbalancedPowerModel)

Constructor for Optimal Power Balance
"""
function _build_mc_mn_opb(pm::AbstractUnbalancedPowerModel)
    for (n, network) in nws(pm)
        variable_mc_generator_power(pm; nw=n)

        for i in ids(pm, :components, nw=n)
            constraint_mc_network_power_balance(pm, i; nw=n)
        end
    end

    objective_mc_min_fuel_cost(pm)
end

"""
    ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})

Ref-extension for opb problem type to add connected components ref
"""
function ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::MathematicalModel)
    apply_pmd!(_ref_add_connected_components!, ref, data)
end


"adds connected components for opb problem type"
function _ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::MathematicalModel{NetworkModel})
    component_sets = calc_connected_components(data)
    ref[:components] = Dict(i => c for (i,c) in enumerate(sort(collect(component_sets); by = length)))
end

"""
    ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})

Ref-extension for opb problem type to add connected components ref.
This method is used when InfrastructureModels passes the raw mathematical data dict
to ref extensions.
"""
function ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    _ref_add_connected_components!(ref, data)
end

function _pmd_ref(ref::Dict{Symbol,<:Any})
    if haskey(ref, :nw)
        return ref
    elseif haskey(ref, :it) && haskey(ref[:it], pmd_it_sym)
        return ref[:it][pmd_it_sym]
    elseif haskey(ref, :it) && haskey(ref[:it], :pmd)
        return ref[:it][:pmd]
    else
        error("could not locate PMD ref; keys(ref)=$(collect(keys(ref)))")
    end
end

"adds connected components for opb problem type from raw mathematical data"
function _ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    pmd_ref = _pmd_ref(ref)

    if get(data, "multinetwork", false) == true || haskey(data, "nw")
        @assert haskey(pmd_ref, :nw) "expected PMD ref to have :nw; keys(pmd_ref)=$(collect(keys(pmd_ref)))"

        for (n, nw_data) in data["nw"]
            nw_id = parse(Int, n)

            if !haskey(nw_data, "bus")
                @warn "skipping connected-components ref for network without bus table" n=n keys=collect(keys(nw_data))
                continue
            end

            nw_ref = get!(pmd_ref[:nw], nw_id, Dict{Symbol,Any}())

            component_sets = _calc_connected_components_math(nw_data)

            nw_ref[:components] = Dict(
                i => c for (i, c) in enumerate(sort(collect(component_sets); by=length))
            )
        end
    else
        component_sets = _calc_connected_components_math(data)

        pmd_ref[:components] = Dict(
            i => c for (i, c) in enumerate(sort(collect(component_sets); by=length))
        )
    end

    return nothing
end