"wraps angles in degrees to 180"
function _wrap_to_180(degrees)
    return degrees - 360*floor.((degrees .+ 180)/360)
end


"wraps angles in radians to pi"
function _wrap_to_pi(radians)
    return radians - 2*pi*floor.((radians .+ pi)/(2*pi))
end


"rolls a 1d array left or right by idx"
function _roll(array::Array{T, 1}, idx::Int; right=true) where T <: Number
    out = Array{T}(undef, size(array))
    pos = idx % length(out)

    if right
        out[1+pos:end] = array[1:end-pos]
        out[1:pos] = array[end-(pos-1):end]
    else
        out[1:end-pos] = array[1+pos:end]
        out[end-(pos-1):end] = array[1:pos]
    end

    return out
end
