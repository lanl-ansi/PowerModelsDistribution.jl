""
function shift_phase_angles!(data::Dict{String,Any})
    assert(haskey(data, "per_unit") && data["per_unit"])

    n = data["phases"]
    angshift = MultiPhaseVector([2*pi/n*(h-1) for h in 1:n])

    if ~haskey(data, "phases_offset") || ~data["phases_offset"]
        for k in keys(data["bus"])
            data["bus"][k]["va"] += angshift
        end
        data["phases_offset"] = true
    elseif data["phases_offset"]
        for k in keys(data["bus"])
            data["bus"][k]["va"] -= angshift
        end
        data["phases_offset"] = false
    end
end
