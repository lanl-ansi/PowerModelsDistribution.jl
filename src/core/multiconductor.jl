# https://stackoverflow.com/questions/39039553/lower-triangular-matrix-in-julia
function vec2utri(v::Vector{T}) where T
    d = length(v)
    n = Int((sqrt(8d+1)+1)/2)
    n*(n-1)/2 == d || error("vec2utri: length of vector is not triangular")
    [ i<j ? v[Int((j-1)*(j-2)/2)+i] : 0 for i=1:n, j=1:n ]
end

function vec2ltri(v::Vector{T}) where T
    vec2utri(v)'
end

function mat2utrivec(m::Matrix{T}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i < j]
end

function mat2ltrivec(m::Matrix{T}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[j,i] for i=1:n, j=1:n if i < j]
end

function make_hermitian_matrix_variable(diag, lowertrianglereal, lowertriangleimag)
    matrixreal = diagm(0 =>   diag) + vec2ltri(lowertrianglereal) + vec2utri(lowertrianglereal)
    matriximag = diagm(0 => 0*diag) + vec2ltri(lowertriangleimag) - vec2utri(lowertriangleimag)
    #TODO if not multiplied with 0, array is not a JuMP type
    return matrixreal, matriximag
end

function make_full_matrix_variable(diag, lowertriangle, uppertriangle)
    matrix = diagm(0 => diag) + vec2ltri(lowertriangle) + vec2utri(uppertriangle)
    return matrix
end
