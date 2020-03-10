""
function _map_math2eng!(data_math)
    @assert get(data_math, "data_model", "mathematical") == "mathematical" "Cannot map data to engineering model: provided data is not a mathematical model"
    @assert haskey(data_math, "map") "Cannot map data to engineering model: no mapping from mathematical to engineering data model is provided"

    data_eng = Dict{<:Any,<:Any}()

    map_keys = sort(keys(data_math["map"]); rev=true)
    for map in map_keys
        # TODO
    end

end


# MAP SOLUTION UP
""
function solution_math2eng(sol_math::Dict, data_math::Dict; make_si::Bool=true, make_deg::Bool=true)
    sol_eng = Dict{String, Any}()

    if make_deg || make_si
        sol_math = solution_make_si(sol_math, data_math; mult_vbase=make_si, mult_sbase=make_si, convert_rad2deg=make_deg)
    end

    map_keys = sort(collect(keys(data_math["map"])); rev=true)
    for map_key in map_keys
        map = data_math["map"][map_key]
        umap_type = map[:unmap_function]

        if     umap_type==:_map_math2eng_sourcebus!
        elseif umap_type==:_map_math2eng_transformer!
        elseif umap_type==:_map_math2eng_root!
        else
            comp_type_eng,  name  = split(map[:from], ".")
            comp_type_math, index = split(map[:to], ".")
            if !haskey(sol_eng, comp_type_eng)
                sol_eng[comp_type_eng] = Dict{String, Any}()
            end
            if haskey(sol_math, comp_type_math)
                sol_eng[comp_type_eng][name] = sol_math[comp_type_math][index]
            else
                #TODO add empty dicts if math object has no solution object?
                sol_eng[comp_type_eng][name] = Dict{String, Any}()
            end
        end
    end

    return sol_eng
end
