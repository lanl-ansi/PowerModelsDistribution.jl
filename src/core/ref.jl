""
function calc_tp_voltage_product_bounds(pm::GenericPowerModel, buspairs; nw::Int=pm.cnw)
    wr_min = Dict([(bp, -Inf) for bp in buspairs])
    wr_max = Dict([(bp,  Inf) for bp in buspairs])
    wi_min = Dict([(bp, -Inf) for bp in buspairs])
    wi_max = Dict([(bp,  Inf) for bp in buspairs])

    for (i, j, c, d) in buspairs
        if i == j
            bus = ref(pm, nw, :bus)[i]
            vm_fr_max = bus["vmax"][c]
            vm_to_max = bus["vmax"][d]
            vm_fr_min = bus["vmin"][c]
            vm_to_min = bus["vmin"][d]
        else
            buspair = ref(pm, nw, :buspairs)[(i, j)]
            vm_fr_max = buspair["vm_fr_max"][c]
            vm_to_max = buspair["vm_to_max"][d]
            vm_fr_min = buspair["vm_fr_min"][c]
            vm_to_min = buspair["vm_to_min"][d]
        end

        wr_max[(i, j, c, d)] =  vm_fr_max * vm_to_max
        wr_min[(i, j, c, d)] = -vm_fr_max * vm_to_max
        wi_max[(i, j, c, d)] =  vm_fr_max * vm_to_max
        wi_min[(i, j, c, d)] = -vm_fr_max * vm_to_max
    end

    return wr_min, wr_max, wi_min, wi_max
end

function find_ref_buses(pm::GenericPowerModel, nw)
    buses = ref(pm, nw, :bus)
    return [b for (b,bus) in buses if bus["bus_type"]==3]
    # return [bus for (b,bus) in buses ]
end

"Adds arcs for TPPM transformers; for dclines and branches this is done in PMs"
function add_arcs_trans!(pm::GenericPowerModel)
    if haskey(ref(pm), :trans)
    # dirty fix add arcs_from/to_trans and bus_arcs_trans
        pm.ref[:nw][0][:arcs_from_trans] = [(i, trans["f_bus"], trans["t_bus"]) for (i,trans) in ref(pm, :trans)]
        pm.ref[:nw][0][:arcs_to_trans] = [(i, trans["t_bus"], trans["f_bus"]) for (i,trans) in ref(pm, :trans)]
        pm.ref[:nw][0][:arcs_trans] = [pm.ref[:nw][0][:arcs_from_trans]..., pm.ref[:nw][0][:arcs_to_trans]...]
        pm.ref[:nw][0][:bus_arcs_trans] = Dict{Int64, Array{Any, 1}}()
        for i in ids(pm, :bus)
            pm.ref[:nw][0][:bus_arcs_trans][i] = [e for e in pm.ref[:nw][0][:arcs_trans] if e[2]==i]
        end
    end
end

function calc_tp_trans_Tvi(pm::GenericPowerModel, i::Int; nw=pm.cnw)
    trans = ref(pm, nw, :trans,  i)
    # transformation matrices
    # Tv and Ti will be compositions of these
    #Tminn = [diagm(0=>ones(Float64, 3)) -ones(Float64, 3)]  # substract neutral from other conductors
    #Tnoz = [diagm(0=>ones(Float64, 3)) zeros(Float64, 3)]   # remove the neutral
    Tbr = [0 0 1; 1 0 0; 0 1 0]                             # barrel roll
    Tdelt  = [1 -1 0; 0 1 -1; -1 0 1]                       # delta transform
    # connection transformers
    perm = trans["conn"][1:3]
    #TODO allow this finally or not?
    @assert(!(perm in ["213", "321", "132"]))
    polarity = trans["conn"][4]
    dyz = trans["conn"][5]
    grounded = length(trans["conn"])>5 && trans["conn"][6]=='n'
    # Tw will contain transformations related to permutation and polarity
    perm_to_trans = Dict(
        "123"=>diagm(0=>ones(Float64, 3)),
        "312"=>Tbr,
        "231"=>Tbr*Tbr
    )
    Tw = perm_to_trans[perm]
    Tw = (polarity=='+') ? Tw : -Tw
    #Tw = diagm(0=>ones(Float64, 3))
    vmult = 1.0 # compensate for change in LN
    if dyz=='y'
        Tv_fr = Tw
        Tv_im = diagm(0=>ones(Float64, 3))
        Ti_fr = Tw
        Ti_im = diagm(0=>ones(Float64, 3))
        # if !grounded
        #     # if not grounded, phase currents should sum to zero
        #     Ti_fr = [Ti_fr; ones(1,3)]
        #     Ti_im = [Ti_im; zeros(1,3)]
        # end
    elseif dyz=='d'
        Tv_fr = Tdelt*Tw
        Tv_im = diagm(0=>ones(Float64, 3))
        Ti_fr = Tw
        Ti_im = Tdelt'
        vmult = sqrt(3)
    elseif dyz=='z'
        #TODO zif-zag here
    else
        error(LOGGER, "Illegal connection specification $dyz.")
    end
    #Tv_fr = diagm(0=>ones(Float64, 3))*Tminn
    Cv_to = trans["vnom_kv"][1]*vmult/trans["vnom_kv"][2]
    # make equations dimensionless
    bkv_fr = ref(pm, nw, :bus, trans["f_bus"], "base_kv")
    bkv_to = ref(pm, nw, :bus, trans["t_bus"], "base_kv")
    Cv_to = Cv_to*bkv_to/bkv_fr
    return (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to)
end
