
#=
""
function apply_func_array(data::Dict{String,Any}, key::String, func)
    if haskey(data, key)
        if isa(data[key], Array)
            data[key] = [func(v) for v in data[key]]
        else
            data[key] = func(data[key])
        end
    end
end


"Transforms network data into per-unit"
function make_per_unit(data::Dict{String,Any})
    if haskey(data, "multiphase") && data["multiphase"]
        if !haskey(data, "per_unit") || data["per_unit"] == false
            data["per_unit"] = true
            mva_base = data["baseMVA"]
            if data["multinetwork"] == true
                for (i,nw_data) in data["nw"]
                    _make_per_unit(nw_data, mva_base)
                end
            else
                 _make_per_unit(data, mva_base)
            end
        end
    else
        error(LOGGER, "this version of make_per_unit requires multiphase data")
    end
end


function _make_per_unit(data::Dict{String,Any}, mva_base::Real)
    rescale      = x -> x/mva_base
    rescale_dual = x -> x*mva_base

    if haskey(data, "bus")
        for (i, bus) in data["bus"]
            apply_func_array(bus, "va", deg2rad)

            apply_func_array(bus, "lam_kcl_r", rescale_dual)
            apply_func_array(bus, "lam_kcl_i", rescale_dual)
        end
    end

    if haskey(data, "load")
        for (i, load) in data["load"]
            apply_func_array(load, "pd", rescale)
            apply_func_array(load, "qd", rescale)
        end
    end

    if haskey(data, "shunt")
        for (i, shunt) in data["shunt"]
            apply_func_array(shunt, "gs", rescale)
            apply_func_array(shunt, "bs", rescale)
        end
    end

    if haskey(data, "branch")
        for (i, branch) in data["branch"]
            apply_func_array(branch, "rate_a", rescale)
            apply_func_array(branch, "rate_b", rescale)
            apply_func_array(branch, "rate_c", rescale)

            apply_func_array(branch, "angmax", deg2rad)
            apply_func_array(branch, "angmin", deg2rad)

            apply_func_array(branch, "mu_sm_fr", rescale_dual)
            apply_func_array(branch, "mu_sm_to", rescale_dual)
        end
    end

    if haskey(data, "gen")
        for (i, gen) in data["gen"]
            apply_func_array(gen, "pg", rescale)
            apply_func_array(gen, "qg", rescale)

            apply_func_array(gen, "pmax", rescale)
            apply_func_array(gen, "pmin", rescale)

            apply_func_array(gen, "qmax", rescale)
            apply_func_array(gen, "qmin", rescale)

            if "model" in keys(gen) && "cost" in keys(gen)
                if gen["model"] == 1
                    for i in 1:2:length(gen["cost"])
                        gen["cost"][i] = gen["cost"][i]/mva_base
                    end
                elseif gen["model"] == 2
                    degree = length(gen["cost"])
                    for (i, item) in enumerate(gen["cost"])
                        gen["cost"][i] = item*mva_base^(degree-i)
                    end
                else
                    warn(LOGGER, "Skipping generator cost model of type $(gen["model"]) in per unit transformation")
                end
            end
        end
    end

end


"Transforms network data into mixed-units (inverse of per-unit)"
function make_mixed_units(data::Dict{String,Any})
    if haskey(data, "multiphase") && data["multiphase"]
        if haskey(data, "per_unit") && data["per_unit"] == true
            data["per_unit"] = false
            mva_base = data["baseMVA"]
            if data["multinetwork"]
                for (i,nw_data) in data["nw"]
                    _make_mixed_units(nw_data, mva_base)
                end
            else
                 _make_mixed_units(data, mva_base)
            end
        end
    else
        error(LOGGER, "this version of make_mixed_units requires multiphase data")
    end
end


function _make_mixed_units(data::Dict{String,Any}, mva_base::Real)
    rescale      = x -> x*mva_base
    rescale_dual = x -> x/mva_base

    if haskey(data, "bus")
        for (i, bus) in data["bus"]
            apply_func_array(bus, "va", rad2deg)

            apply_func_array(bus, "lam_kcl_r", rescale_dual)
            apply_func_array(bus, "lam_kcl_i", rescale_dual)
        end
    end

    if haskey(data, "load")
        for (i, load) in data["load"]
            apply_func_array(load, "pd", rescale)
            apply_func_array(load, "qd", rescale)
        end
    end

    if haskey(data, "shunt")
        for (i, shunt) in data["shunt"]
            apply_func_array(shunt, "gs", rescale)
            apply_func_array(shunt, "bs", rescale)
        end
    end

    if haskey(data, "branch")
        for (i, branch) in data["branch"]
            apply_func_array(branch, "rate_a", rescale)
            apply_func_array(branch, "rate_b", rescale)
            apply_func_array(branch, "rate_c", rescale)

            apply_func_array(branch, "angmax", rad2deg)
            apply_func_array(branch, "angmin", rad2deg)

            apply_func_array(branch, "mu_sm_fr", rescale_dual)
            apply_func_array(branch, "mu_sm_to", rescale_dual)
        end
    end

    if haskey(data, "gen")
        for (i, gen) in data["gen"]
            apply_func_array(gen, "pg", rescale)
            apply_func_array(gen, "qg", rescale)

            apply_func_array(gen, "pmax", rescale)
            apply_func_array(gen, "pmin", rescale)

            apply_func_array(gen, "qmax", rescale)
            apply_func_array(gen, "qmin", rescale)

            if "model" in keys(gen) && "cost" in keys(gen)
                if gen["model"] == 1
                    for i in 1:2:length(gen["cost"])
                        gen["cost"][i] = gen["cost"][i]*mva_base
                    end
                elseif gen["model"] == 2
                    degree = length(gen["cost"])
                    for (i, item) in enumerate(gen["cost"])
                        gen["cost"][i] = item/mva_base^(degree-i)
                    end
                else
                    warn(LOGGER, "Skipping generator cost model of type $(gen["model"]) in mixed units transformation")
                end
            end
        end
    end

end
=#