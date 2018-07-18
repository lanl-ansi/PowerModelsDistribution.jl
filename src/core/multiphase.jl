# https://stackoverflow.com/questions/39039553/lower-triangular-matrix-in-julia
function vec2ltri{T}(v::Vector{T})
    d = length(v)
    n = Int((sqrt(8d+1)+1)/2)
    n*(n-1)/2 == d || error("vec2ltri: length of vector is not triangular")
    [ i>j ? v[Int((2n-j)*(j-1)/2)+i-j] : 0 for i=1:n, j=1:n ]
end

# https://stackoverflow.com/questions/39039553/lower-triangular-matrix-in-julia
function vec2utri{T}(v::Vector{T})
    d = length(v)
    n = Int((sqrt(8d+1)+1)/2)
    n*(n-1)/2 == d || error("vec2utri: length of vector is not triangular")
    [ i<j ? v[Int((j-1)*(j-2)/2)+i] : 0 for i=1:n, j=1:n ]
end

function mat2utrivec{T}(m::Matrix{T})
    assert(size(m,1) == size(m,2))
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i < j]
end


function mat2ltrivec{T}(m::Matrix{T})
    assert(size(m,1) == size(m,2))
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i > j]
end

function make_hermitian_matrix_variable(diag, lowertrianglereal, lowertriangleimag)
    matrixreal = diagm(  diag) + vec2ltri(lowertrianglereal) + vec2utri(lowertrianglereal)
    matriximag = diagm(0*diag) + vec2ltri(lowertriangleimag) - vec2utri(lowertriangleimag)
    #TODO if not multiplied with 0, array is not a JuMP type
    return matrixreal, matriximag
end

function make_full_matrix_variable(diag, lowertriangle, uppertriangle)
    matrix = diagm(diag) + vec2ltri(lowertriangle) + vec2utri(uppertriangle)
    return matrix
end
