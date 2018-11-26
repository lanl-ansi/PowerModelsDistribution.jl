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
    #TODO clean up
    matrixreal = []
    if length(diag) == 3
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]
        ]
    elseif length(diag) == 4
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2]    lowertrianglereal[4];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3]    lowertrianglereal[5];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]                 lowertrianglereal[6];
        lowertrianglereal[4]    lowertrianglereal[5]    lowertrianglereal[6]    diag[4]
        ]
    elseif length(diag) == 5
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2]    lowertrianglereal[4]    lowertrianglereal[7];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3]    lowertrianglereal[5]    lowertrianglereal[8];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]                 lowertrianglereal[6]    lowertrianglereal[9];
        lowertrianglereal[4]    lowertrianglereal[5]    lowertrianglereal[6]    diag[4]                 lowertrianglereal[10];
        lowertrianglereal[7]    lowertrianglereal[8]    lowertrianglereal[9]    lowertrianglereal[10]    diag[5]
        ]
    end

    matriximag = []
    if length(diag) == 3
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2];
        lowertriangleimag[1]    0                       -lowertriangleimag[3];
        lowertriangleimag[2]    lowertriangleimag[3]    0
        ]
    elseif length(diag) == 4
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2]   -lowertriangleimag[4];
        lowertriangleimag[1]    0                       -lowertriangleimag[3]   -lowertriangleimag[5];
        lowertriangleimag[2]    lowertriangleimag[3]    0                       -lowertriangleimag[6];
        lowertriangleimag[4]    lowertriangleimag[5]    lowertriangleimag[6]    0
        ]
    elseif length(diag) == 5
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2]   -lowertriangleimag[4]   -lowertriangleimag[7];
        lowertriangleimag[1]    0                       -lowertriangleimag[3]   -lowertriangleimag[5]   -lowertriangleimag[8];
        lowertriangleimag[2]    lowertriangleimag[3]    0                       -lowertriangleimag[6]   -lowertriangleimag[9];
        lowertriangleimag[4]    lowertriangleimag[5]    lowertriangleimag[6]    0                       -lowertriangleimag[10];
        lowertriangleimag[7]    lowertriangleimag[8]    lowertriangleimag[9]    lowertriangleimag[10]    0
        ]
    end
    return matrixreal, matriximag
end

function make_full_matrix_variable(diag, lowertriangle, uppertriangle)
    #TODO clean up
    matrix = []
    if length(diag) == 3
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2];
        lowertriangle[1]    diag[2]             uppertriangle[3];
        lowertriangle[2]    lowertriangle[3]    diag[3]
        ]
    elseif length(diag) == 4
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2]    uppertriangle[4];
        lowertriangle[1]    diag[2]             uppertriangle[3]    uppertriangle[5];
        lowertriangle[2]    lowertriangle[3]    diag[3]             uppertriangle[6];
        lowertriangle[4]    lowertriangle[5]    lowertriangle[6]    diag[4]
        ]
    elseif length(diag) == 5
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2]    uppertriangle[4]    uppertriangle[7];
        lowertriangle[1]    diag[2]             uppertriangle[3]    uppertriangle[5]    uppertriangle[8];
        lowertriangle[2]    lowertriangle[3]    diag[3]             uppertriangle[6]    uppertriangle[9];
        lowertriangle[4]    lowertriangle[5]    lowertriangle[6]    diag[4]             uppertriangle[10]
        lowertriangle[7]    lowertriangle[8]    lowertriangle[9]    lowertriangle[10]    diag[5]
        ]
    end
    # matrix = diagm(0 => diag) + vec2ltri(lowertriangle) + vec2utri(uppertriangle)
    return matrix
end
