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
function solution_math2eng(solution_math::Dict, data_math::Dict; make_si::Bool=true, make_deg::Bool=true)
    solution_eng = Dict{String, Any}()

    if make_deg || make_si
        solution_math = solution_make_si(solution_math, data_math; mult_vbase=make_si, mult_sbase=make_si, convert_rad2deg=make_deg)
    end

    map_keys = sort(collect(keys(data_math["map"])); rev=true)
    for map_key in map_keys
        map = data_math["map"][map_key]
        umap_type = map[:unmap_function]

        if     umap_type==:_map_math2eng_sourcebus!
            # nothing to do for the solution
        elseif umap_type==:_map_math2eng_root!
            # nothing to do for the solution
        elseif umap_type==:_map_math2eng_transformer!
            _,  name  = split(map[:from_id], ".")
            trans_2wa_ids = [index for (comp_type, index) in split.(map[:to_id], ".", limit=2) if comp_type=="transformer"]

            if !haskey(solution_eng, "transformer")
                solution_eng["transformer"] = Dict{String, Any}()
            end

            trans = Dict()
            prop_map = Dict("pf"=>"p", "qf"=>"q")
            for (prop_from, prop_to) in prop_map
                trans[prop_to] = [get(solution_math["transformer"][id], prop_from, NaN) for id in trans_2wa_ids]
            end

            solution_eng["transformer"][name] = trans

        else
            comp_type_eng,  name  = split(map[:from], ".")
            comp_type_math, index = split(map[:to], ".")
            if !haskey(solution_eng, comp_type_eng)
                solution_eng[comp_type_eng] = Dict{String, Any}()
            end
            if haskey(solution_math, comp_type_math)
                solution_eng[comp_type_eng][name] = solution_math[comp_type_math][index]
            else
                #TODO add empty dicts if math object has no solution object?
                solution_eng[comp_type_eng][name] = Dict{String, Any}()
            end
        end
    end

    return solution_eng
end
