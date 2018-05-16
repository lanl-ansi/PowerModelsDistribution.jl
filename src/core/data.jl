# Three-phase specific branch functions


""
function calc_branch_t(branch::Dict{String,Any})
    tap_ratio = branch["tap"]
    angle_shift = branch["shift"]

    tr = map(*, tap_ratio, map(cos, angle_shift))
    ti = map(*, tap_ratio, map(sin, angle_shift))

    trout = Array{Float64}(length(tr))
    for i in range(1, length(tr))
        trout[i] = tr[i]
    end

    tiout = Array{Float64}(length(ti))
    for i in range(1, length(ti))
        tiout[i] = ti[i]
    end


    return trout, tiout
end


""
function calc_branch_y(branch::Dict{String,Any})
    r = branch["br_r"]
    x = branch["br_x"]

    ym = map(+, map(*, r, r), map(*, x, x))

    g = map(PMs._div_zero, r, ym)
    b = map(-, map(PMs._div_zero, x, ym))

    gout = Matrix{Float64}(size(g)...)
    for i in range(1, size(g)[1])
        for j in range(1, size(g)[2])
            gout[i, j] = g[i, j]
        end
    end

    bout = Matrix{Float64}(size(b)...)
    for i in range(1, size(b)[1])
        for j in range(1, size(b)[2])
            bout[i, j] = b[i, j]
        end
    end

    return gout, bout
end