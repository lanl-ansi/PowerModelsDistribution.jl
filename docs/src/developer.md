# Developer Documentation

## Variable and parameter naming scheme

- suffix `_mx` for matrix-based formulation

### Voltage variables (Hermitian matrix)

Defining voltage product $\mathbf{W}_{i} = \mathbf{V}_i \cdot (\mathbf{V}_i)^H$, with $\mathbf{V}_i =  \begin{pmatrix} V_{i,a}\\ V_{i,b} \\ V_{i,c} \end{pmatrix} \in \mathbb{C}^{3\times1}$ and $\mathbf{W}_{i} \in \mathbb{C}^{3\times3}$ where
$\mathbf{W}_{i} = \mathbf{W}^{\text{re}}_{i} + j\mathbf{W}^{\text{re}}_{i} $
= `W` = `Wr` $ + j\cdot$ `Wi`:
- `W` (capitalized, short for VV):   complex matrix of voltage vector outer product (V$^2$)
- `Wr` (capitalized, short for VVr): real part of matrix of voltage vector outer product (V$^2$)
- `Wi` (capitalized, short for VVi): imaginary part of matrix of voltage vector outer product (V$^2$)

The diagonal of $\mathbf{W}_{i}$ is $diag(\mathbf{W}_{i}) = \begin{pmatrix} V_{i,a} V^*_{i,a}\\ V_{i,b} V^*_{i,b} \\ V_{i,c} V^*_{i,c} \end{pmatrix}$ = `w`
- `w` (lowercase, short for vv): diagonal vector of voltage vector outer product (V$^2$)
