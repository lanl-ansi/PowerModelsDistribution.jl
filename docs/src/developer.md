# Developer Documentation

## Variable and parameter naming scheme

### Voltage variables (Hermitian matrix)

Defining voltage product $\mathbf{W}_{i} = \mathbf{V}_i \cdot (\mathbf{V}_i)^H$, with $\mathbf{V}_i =  \begin{pmatrix} v_{i,a}\\ v_{i,b} \\ v_{i,c} \end{pmatrix}$ and
$\mathbf{W}_{i} = \mathbf{W}^{\text{re}}_{i} + j\mathbf{W}^{\text{re}}_{i} $
= `W` = `Wr` $ + j\cdot$ `Wi`:
- `W` (short for VV):   complex matrix of voltage vector outer product (V$^2$)
- `Wr` (short for VVr): real part of matrix of voltage vector outer product (V$^2$)
- `Wi` (short for VVi): imaginary part of matrix of voltage vector outer product (V$^2$)

The diagonal of $\mathbf{W}_{i}$ is $diag(\mathbf{W}_{i}) = \begin{pmatrix} v_{i,a}v^*_{i,a}\\ v_{i,b} v^*_{i,b} \\ v_{i,c} v^*_{i,c} \end{pmatrix}$ = `w`
- `w` (short for vv): diagonal vector of voltage vector outer product (V$^2$)
